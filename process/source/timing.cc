#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <boost/program_options.hpp>

#include "TimingEstimator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <string>

namespace {

constexpr int maxframes = 65536;

int
timing_do_channel(int device, int fifo, int column, int pixel)
{
  static const int eo2do[32] = {
    22, 20, 18, 16, 24, 26, 28, 30,
    25, 27, 29, 31, 23, 21, 19, 17,
    9,  11, 13, 15, 7,  5,  3,  1,
    6,  4,  2,  0,  8,  10, 12, 14
  };

  if (device != 200 || fifo < 0 || fifo > 7)
    return -1;

  int eoch = pixel + 4 * column;
  if (eoch < 0 || eoch >= 32)
    return -1;
  return eo2do[eoch];
}

struct timing_hit_t {
  int slot;
  double time;
};

bool
has_branch(TTree *tree, const std::string &name)
{
  if (tree && tree->GetBranch(name.c_str()))
    return true;

  std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
  return false;
}

bool
check_output_branch(TTree *tree, const std::string &name)
{
  if (!tree->GetBranch(name.c_str()))
    return true;

  std::cerr << "ERROR: output timing branch already exists: " << name << std::endl;
  return false;
}

bool
copy_tree_if_present(TFile *fin, TFile *fout, const char *name)
{
  auto tree = (TTree *)fin->Get(name);
  if (!tree)
    return true;

  fout->cd();
  auto clone = tree->CloneTree(-1, "fast");
  if (!clone) {
    std::cerr << "ERROR: failed to copy tree '" << name << "'" << std::endl;
    return false;
  }
  return true;
}

}

bool
timing(const std::string &filename,
       const std::string &outfilename)
{
  auto fin = TFile::Open(filename.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << filename << std::endl;
    return false;
  }

  auto tin = (TTree *)fin->Get("frames");
  if (!tin) {
    std::cerr << "ERROR: could not find 'frames' tree in input file" << std::endl;
    fin->Close();
    return false;
  }

  for (auto name : {
         "nframes",
         "ntiminghits",
         "timing_frame_start",
         "timing_frame_nhits",
         "timing_device",
         "timing_fifo",
         "timing_type",
         "timing_column",
         "timing_pixel",
         "timing_time"
       }) {
    if (!has_branch(tin, name)) {
      fin->Close();
      return false;
    }
  }

  for (auto name : {
         "timing_valid",
         "T0",
         "sigma0",
         "T1",
         "sigma1",
         "T",
         "sigmaT"
       }) {
    if (!check_output_branch(tin, name)) {
      fin->Close();
      return false;
    }
  }

  TTreeReader reader(tin);
  TTreeReaderValue<int> nframes(reader, "nframes");
  TTreeReaderValue<int> ntiminghits(reader, "ntiminghits");
  TTreeReaderArray<int> timing_frame_start(reader, "timing_frame_start");
  TTreeReaderArray<int> timing_frame_nhits(reader, "timing_frame_nhits");
  TTreeReaderArray<int> timing_device(reader, "timing_device");
  TTreeReaderArray<int> timing_fifo(reader, "timing_fifo");
  TTreeReaderArray<int> timing_type(reader, "timing_type");
  TTreeReaderArray<int> timing_column(reader, "timing_column");
  TTreeReaderArray<int> timing_pixel(reader, "timing_pixel");
  TTreeReaderArray<double> timing_time(reader, "timing_time");

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    fin->Close();
    return false;
  }

  fout->cd();
  auto tout = tin->CloneTree(-1, "fast");
  if (!tout) {
    std::cerr << "ERROR: could not clone 'frames' tree" << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  auto timing_valid = std::make_unique<int[]>(maxframes);
  auto T0 = std::make_unique<double[]>(maxframes);
  auto sigma0 = std::make_unique<double[]>(maxframes);
  auto T1 = std::make_unique<double[]>(maxframes);
  auto sigma1 = std::make_unique<double[]>(maxframes);
  auto T = std::make_unique<double[]>(maxframes);
  auto sigmaT = std::make_unique<double[]>(maxframes);

  int nframes_out = 0;
  tout->SetBranchAddress("nframes", &nframes_out);
  auto b_timing_valid = tout->Branch("timing_valid", timing_valid.get(), "timing_valid[nframes]/I");
  auto b_T0 = tout->Branch("T0", T0.get(), "T0[nframes]/D");
  auto b_sigma0 = tout->Branch("sigma0", sigma0.get(), "sigma0[nframes]/D");
  auto b_T1 = tout->Branch("T1", T1.get(), "T1[nframes]/D");
  auto b_sigma1 = tout->Branch("sigma1", sigma1.get(), "sigma1[nframes]/D");
  auto b_T = tout->Branch("T", T.get(), "T[nframes]/D");
  auto b_sigmaT = tout->Branch("sigmaT", sigmaT.get(), "sigmaT[nframes]/D");

  TimingEstimator estimator;
  auto nan = std::numeric_limits<double>::quiet_NaN();

  Long64_t nspills = 0;
  Long64_t nframes_total = 0;
  Long64_t nframes_timed = 0;
  Long64_t nframes_incomplete = 0;
  Long64_t nduplicate_hits = 0;
  Long64_t nignored_hits = 0;

  auto entries = tin->GetEntries();
  while (reader.Next()) {
    ++nspills;
    nframes_out = *nframes;

    if (*nframes < 0 || *nframes > maxframes) {
      std::cerr << "ERROR: invalid nframes=" << *nframes
                << " at spill entry " << (nspills - 1) << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    for (int iframe = 0; iframe < *nframes; ++iframe) {
      timing_valid[iframe] = 0;
      T0[iframe] = nan;
      sigma0[iframe] = nan;
      T1[iframe] = nan;
      sigma1[iframe] = nan;
      T[iframe] = nan;
      sigmaT[iframe] = nan;

      ++nframes_total;

      int first = timing_frame_start[iframe];
      int nhits = timing_frame_nhits[iframe];
      if (*ntiminghits < 0 || first < 0 || nhits < 0 ||
          first + nhits > *ntiminghits) {
        std::cerr << "ERROR: invalid timing frame indexing"
                  << " spill_entry=" << (nspills - 1)
                  << " frame=" << iframe
                  << " first=" << first
                  << " nhits=" << nhits
                  << " total=" << *ntiminghits
                  << std::endl;
        fout->Close();
        fin->Close();
        return false;
      }

      std::array<double, 64> times{};
      std::array<bool, 64> found{};
      std::vector<timing_hit_t> selected;
      selected.reserve(nhits);

      for (int ihit = first; ihit < first + nhits; ++ihit) {
        if (timing_type[ihit] != 1) {
          ++nignored_hits;
          continue;
        }

        int doch = timing_do_channel(timing_device[ihit],
                                     timing_fifo[ihit],
                                     timing_column[ihit],
                                     timing_pixel[ihit]);
        if (doch < 0) {
          ++nignored_hits;
          continue;
        }

        int offset = timing_fifo[ihit] < 4 ? 0 : 32;
        selected.push_back({offset + doch, timing_time[ihit]});
      }

      std::sort(selected.begin(), selected.end(),
                [](const timing_hit_t &a, const timing_hit_t &b) {
                  if (a.slot != b.slot)
                    return a.slot < b.slot;
                  return a.time < b.time;
                });

      for (const auto &hit : selected) {
        if (!found[hit.slot]) {
          times[hit.slot] = hit.time;
          found[hit.slot] = true;
        } else {
          ++nduplicate_hits;
        }
      }

      bool complete = true;
      for (bool ok : found)
        complete = complete && ok;

      if (!complete) {
        ++nframes_incomplete;
        continue;
      }

      auto result = estimator.estimate(times);
      timing_valid[iframe] = 1;
      T0[iframe] = result.T0;
      sigma0[iframe] = result.sigma0;
      T1[iframe] = result.T1;
      sigma1[iframe] = result.sigma1;
      T[iframe] = result.T;
      sigmaT[iframe] = result.sigmaT;
      ++nframes_timed;
    }

    b_timing_valid->Fill();
    b_T0->Fill();
    b_sigma0->Fill();
    b_T1->Fill();
    b_sigma1->Fill();
    b_T->Fill();
    b_sigmaT->Fill();
  }

  if (tout->GetEntries() != entries) {
    std::cerr << "ERROR: output entry count mismatch: input=" << entries
              << " output=" << tout->GetEntries() << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  for (auto branch : {b_timing_valid, b_T0, b_sigma0, b_T1, b_sigma1, b_T, b_sigmaT}) {
    if (branch->GetEntries() != entries) {
      std::cerr << "ERROR: timing branch entry count mismatch: branch="
                << branch->GetName()
                << " entries=" << branch->GetEntries()
                << " expected=" << entries
                << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }
  }

  if (!copy_tree_if_present(fin, fout, "spill_participation")) {
    fout->Close();
    fin->Close();
    return false;
  }

  if (fout->Write() <= 0) {
    std::cerr << "ERROR: failed writing output file: " << outfilename << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  fout->Close();
  fin->Close();

  std::cout << "input spill entries:        " << entries << std::endl;
  std::cout << "output spill entries:       " << entries << std::endl;
  std::cout << "frames processed:           " << nframes_total << std::endl;
  std::cout << "frames with timing estimate:" << nframes_timed << std::endl;
  std::cout << "incomplete timing frames:   " << nframes_incomplete << std::endl;
  std::cout << "duplicate timing hits:      " << nduplicate_hits << std::endl;
  std::cout << "ignored timing hits:        " << nignored_hits << std::endl;

  return true;
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string output;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input triggered ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file with timing branches");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) {
      std::cout << options << std::endl;
      return 0;
    }
    po::notify(vm);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    std::cerr << options << std::endl;
    return 1;
  }

  return timing(input, output) ? 0 : 1;
}
