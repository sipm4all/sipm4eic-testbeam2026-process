#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <boost/program_options.hpp>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace po = boost::program_options;

int main(int argc, char **argv)
{
  std::string input, output;
  double spatial = 3.2, time = 5.;
  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help")
    ("input,i", po::value<std::string>(&input)->required(), "input ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("spatial-distance", po::value<double>(&spatial)->default_value(3.2), "maximum XY distance [mm]")
    ("time-distance", po::value<double>(&time)->default_value(5.), "maximum time difference [native units]");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) { std::cout << options << std::endl; return 0; }
    po::notify(vm);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl << options << std::endl;
    return 2;
  }
  if (spatial <= 0. || time < 0.) {
    std::cerr << "ERROR: clustering distances must be non-negative (spatial > 0)" << std::endl;
    return 2;
  }

  auto fin = TFile::Open(input.c_str(), "READ");
  if (!fin || fin->IsZombie()) { std::cerr << "ERROR: cannot open input: " << input << std::endl; return 1; }
  auto frames = dynamic_cast<TTree *>(fin->Get("frames"));
  auto cherenkov = dynamic_cast<TTree *>(fin->Get("cherenkov"));
  if (!frames || !cherenkov) {
    std::cerr << "ERROR: input must contain 'frames' and 'cherenkov' trees" << std::endl;
    return 1;
  }
  for (const char *name : {"spill", "nhits", "x", "y", "time"}) {
    TTree *tree = (std::string(name) == "spill") ? frames : cherenkov;
    if (!tree->GetBranch(name)) {
      std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
      return 1;
    }
  }
  if (frames->GetEntries() != cherenkov->GetEntries()) {
    std::cerr << "ERROR: frames/cherenkov entry-count mismatch" << std::endl;
    return 1;
  }

  auto fout = TFile::Open(output.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) { std::cerr << "ERROR: cannot create output: " << output << std::endl; return 1; }
  // Copy every input tree byte-for-byte at the entry level before adding clusters.
  TIter keys(fin->GetListOfKeys());
  while (auto key = dynamic_cast<TKey *>(keys())) {
    if (std::string(key->GetClassName()) != "TTree") continue;
    // clusters is recomputed below; do not copy an old version and create a
    // second ROOT key with the same name.
    if (std::string(key->GetName()) == "clusters") continue;
    auto tree = dynamic_cast<TTree *>(key->ReadObj());
    auto clone = tree ? tree->CloneTree(-1, "fast") : nullptr;
    if (!clone || !clone->Write()) {
      std::cerr << "ERROR: failed to clone tree '" << key->GetName() << "'" << std::endl;
      return 1;
    }
  }

  fout->cd();
  TTree clusters("clusters", "space-time hit clusters, one entry per frame");
  constexpr std::size_t max_clusters = 65535;
  uint16_t nclusters = 0;
  uint8_t size[max_clusters];
  float x[max_clusters], y[max_clusters], cluster_time[max_clusters];
  clusters.Branch("nclusters", &nclusters, "nclusters/s");
  clusters.Branch("size", size, "size[nclusters]/b");
  clusters.Branch("x", x, "x[nclusters]/F");
  clusters.Branch("y", y, "y[nclusters]/F");
  clusters.Branch("time", cluster_time, "time[nclusters]/F");

  TTreeReader reader(cherenkov);
  TTreeReaderValue<UShort_t> nhits(reader, "nhits");
  TTreeReaderArray<Float_t> hit_x(reader, "x");
  TTreeReaderArray<Float_t> hit_y(reader, "y");
  TTreeReaderArray<Float_t> hit_time(reader, "time");
  const Long64_t entries = frames->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (!reader.Next()) { std::cerr << "ERROR: failed to read cherenkov entry " << entry << std::endl; return 1; }
    const std::size_t n = *nhits;
    if (hit_x.GetSize() < n || hit_y.GetSize() < n || hit_time.GetSize() < n) {
      std::cerr << "ERROR: invalid hit arrays at frame " << entry << std::endl; return 1;
    }
    std::vector<bool> visited(n, false);
    nclusters = 0;
    for (std::size_t seed = 0; seed < n; ++seed) {
      if (visited[seed]) continue;
      if (nclusters == max_clusters) { std::cerr << "ERROR: too many clusters at frame " << entry << std::endl; return 1; }
      std::vector<std::size_t> pending{seed};
      visited[seed] = true;
      double sx = 0., sy = 0., st = 0.;
      std::size_t count = 0;
      while (!pending.empty()) {
        const auto current = pending.back(); pending.pop_back();
        sx += hit_x[current]; sy += hit_y[current]; st += hit_time[current]; ++count;
        for (std::size_t other = 0; other < n; ++other) {
          if (visited[other]) continue;
          const double dx = hit_x[current] - hit_x[other];
          const double dy = hit_y[current] - hit_y[other];
          if (std::hypot(dx, dy) <= spatial &&
              std::abs(hit_time[current] - hit_time[other]) <= time) {
            visited[other] = true;
            pending.push_back(other);
          }
        }
      }
      if (count > std::numeric_limits<uint8_t>::max()) {
        std::cerr << "ERROR: cluster size " << count
                  << " exceeds uint8_t capacity at frame " << entry << std::endl;
        return 1;
      }
      size[nclusters] = static_cast<uint8_t>(count);
      x[nclusters] = static_cast<float>(sx / count);
      y[nclusters] = static_cast<float>(sy / count);
      cluster_time[nclusters] = static_cast<float>(st / count);
      ++nclusters;
    }
    clusters.Fill();
  }
  fout->Close();
  fin->Close();
  std::cout << "frames processed: " << entries << std::endl
            << "spatial distance [mm]: " << spatial << std::endl
            << "time distance [native units]: " << time << std::endl
            << "output: " << output << std::endl;
  return 0;
}
