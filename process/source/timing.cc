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
  auto ttiming_in = (TTree *)fin->Get("timing");
  if (!tin || !ttiming_in) {
    std::cerr << "ERROR: input must contain 'frames' and 'timing' trees" << std::endl;
    fin->Close();
    return false;
  }

  for (auto name : {"spill"}) {
    if (!has_branch(tin, name)) {
      fin->Close();
      return false;
    }
  }
  for (auto name : {"nhits", "device", "fifo", "column", "pixel", "time"}) {
    if (!has_branch(ttiming_in, name)) {
      fin->Close();
      return false;
    }
  }
  for (auto name : {"timing_valid", "T0", "sigma0", "T1", "sigma1", "T", "sigmaT"}) {
    if (!check_output_branch(ttiming_in, name)) {
      fin->Close();
      return false;
    }
  }

  auto ttrigger_in = (TTree *)fin->Get("trigger");
  auto tcherenkov_in = (TTree *)fin->Get("cherenkov");
  auto entries = tin->GetEntries();
  if (ttiming_in->GetEntries() != entries ||
      (ttrigger_in && ttrigger_in->GetEntries() != entries) ||
      (tcherenkov_in && tcherenkov_in->GetEntries() != entries)) {
    std::cerr << "ERROR: frame-aligned tree entry-count mismatch"
              << " frames=" << entries
              << " timing=" << ttiming_in->GetEntries();
    if (ttrigger_in)
      std::cerr << " trigger=" << ttrigger_in->GetEntries();
    if (tcherenkov_in)
      std::cerr << " cherenkov=" << tcherenkov_in->GetEntries();
    std::cerr << std::endl;
    fin->Close();
    return false;
  }

  TTreeReader timing_reader(ttiming_in);
  TTreeReaderValue<UShort_t> timing_nhits(timing_reader, "nhits");
  TTreeReaderArray<UChar_t> timing_device(timing_reader, "device");
  TTreeReaderArray<UChar_t> timing_fifo(timing_reader, "fifo");
  TTreeReaderArray<UChar_t> timing_column(timing_reader, "column");
  TTreeReaderArray<UChar_t> timing_pixel(timing_reader, "pixel");
  TTreeReaderArray<Float_t> timing_time(timing_reader, "time");

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    fin->Close();
    return false;
  }

  fout->cd();
  auto tout = tin->CloneTree(0);
  auto ttrigger_out = ttrigger_in ? ttrigger_in->CloneTree(0) : nullptr;
  auto ttiming_out = ttiming_in->CloneTree(0);
  auto tcherenkov_out = tcherenkov_in ? tcherenkov_in->CloneTree(0) : nullptr;
  if (!tout || !ttiming_out || (ttrigger_in && !ttrigger_out) ||
      (tcherenkov_in && !tcherenkov_out)) {
    std::cerr << "ERROR: failed to create output frame trees" << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  bool timing_valid = false;
  float T0 = 0.;
  float sigma0 = 0.;
  float T1 = 0.;
  float sigma1 = 0.;
  float T = 0.;
  float sigmaT = 0.;
  auto b_timing_valid = ttiming_out->Branch("timing_valid", &timing_valid, "timing_valid/O");
  auto b_T0 = ttiming_out->Branch("T0", &T0, "T0/F");
  auto b_sigma0 = ttiming_out->Branch("sigma0", &sigma0, "sigma0/F");
  auto b_T1 = ttiming_out->Branch("T1", &T1, "T1/F");
  auto b_sigma1 = ttiming_out->Branch("sigma1", &sigma1, "sigma1/F");
  auto b_T = ttiming_out->Branch("T", &T, "T/F");
  auto b_sigmaT = ttiming_out->Branch("sigmaT", &sigmaT, "sigmaT/F");

  TimingEstimator estimator;
  auto nan = std::numeric_limits<float>::quiet_NaN();
  Long64_t nframes_timed = 0;
  Long64_t nframes_incomplete = 0;
  Long64_t nduplicate_hits = 0;
  Long64_t nignored_hits = 0;

  for (Long64_t iframe = 0; iframe < entries; ++iframe) {
    if (tin->GetEntry(iframe) <= 0 || !timing_reader.Next()) {
      std::cerr << "ERROR: failed to read synchronized frame entry " << iframe << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    timing_valid = false;
    T0 = nan;
    sigma0 = nan;
    T1 = nan;
    sigma1 = nan;
    T = nan;
    sigmaT = nan;

    const int nhits = static_cast<int>(*timing_nhits);
    if (timing_device.GetSize() < nhits || timing_fifo.GetSize() < nhits ||
        timing_column.GetSize() < nhits || timing_pixel.GetSize() < nhits ||
        timing_time.GetSize() < nhits) {
      std::cerr << "ERROR: invalid timing hit array at frame entry " << iframe
                << " nhits=" << nhits << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    std::array<double, 64> times{};
    std::array<bool, 64> found{};
    std::vector<timing_hit_t> selected;
    selected.reserve(nhits);

    for (int ihit = 0; ihit < nhits; ++ihit) {
      int doch = timing_do_channel(timing_device[ihit], timing_fifo[ihit],
                                   timing_column[ihit], timing_pixel[ihit]);
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
    } else {
      auto result = estimator.estimate(times);
      timing_valid = true;
      T0 = static_cast<float>(result.T0);
      sigma0 = static_cast<float>(result.sigma0);
      T1 = static_cast<float>(result.T1);
      sigma1 = static_cast<float>(result.sigma1);
      T = static_cast<float>(result.T);
      sigmaT = static_cast<float>(result.sigmaT);
      ++nframes_timed;
    }

    if (ttrigger_in)
      ttrigger_in->GetEntry(iframe);
    if (tcherenkov_in)
      tcherenkov_in->GetEntry(iframe);
    tout->Fill();
    ttiming_out->Fill();
    if (ttrigger_out)
      ttrigger_out->Fill();
    if (tcherenkov_out)
      tcherenkov_out->Fill();
  }

  for (auto tree : {tout, ttiming_out, ttrigger_out, tcherenkov_out}) {
    if (tree && tree->GetEntries() != entries) {
      std::cerr << "ERROR: output tree entry-count mismatch: " << tree->GetName()
                << "=" << tree->GetEntries() << " expected=" << entries << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }
  }
  for (auto branch : {b_timing_valid, b_T0, b_sigma0, b_T1, b_sigma1, b_T, b_sigmaT}) {
    if (branch->GetEntries() != entries) {
      std::cerr << "ERROR: timing branch entry-count mismatch: " << branch->GetName()
                << "=" << branch->GetEntries() << " expected=" << entries << std::endl;
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
  std::cout << "input frame entries:        " << entries << std::endl;
  std::cout << "output frame entries:       " << entries << std::endl;
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
