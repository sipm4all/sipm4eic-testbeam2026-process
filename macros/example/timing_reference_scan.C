#include "../lib/trigger_reader.h"

#include <TFile.h>
#include <TH2D.h>
#include <TString.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

namespace {

constexpr int nchannels = 32;

const int eo2do[nchannels] = {22, 20, 18, 16, 24, 26, 28, 30,
                              25, 27, 29, 31, 23, 21, 19, 17,
                              9,  11, 13, 15, 7,  5,  3,  1,
                              6,  4,  2,  0,  8,  10, 12, 14};

int timing_detector(const hit_t &hit)
{
  if (hit.device != 200)
    return -1;
  if (hit.fifo >= 0 && hit.fifo <= 3)
    return 0;
  if (hit.fifo >= 4 && hit.fifo <= 7)
    return 1;
  return -1;
}

int electronics_channel(const hit_t &hit)
{
  int channel = hit.pixel + 4 * hit.column;
  if (channel < 0 || channel >= nchannels)
    return -1;
  return channel;
}

double do_distance(int a, int b)
{
  int ax = a % 4;
  int ay = a / 4;
  int bx = b % 4;
  int by = b / 4;
  int dx = ax - bx;
  int dy = ay - by;
  return std::sqrt(dx * dx + dy * dy);
}

double reflected_distance(int a, int b, double pitch, double thickness)
{
  int ax = a % 4;
  int ay = a / 4;
  int bx = b % 4;
  int by = b / 4;
  double dx = pitch * (ax - bx);
  double dy = pitch * (ay - by);
  return 2. * std::sqrt(dx * dx + dy * dy + thickness * thickness);
}

bool load_timing0_offsets(const char *filename, double offset[nchannels])
{
  bool loaded[nchannels] = {};

  for (int i = 0; i < nchannels; ++i)
    offset[i] = 0.;

  std::ifstream in(filename);
  if (!in) {
    std::cerr << "ERROR: could not open calibration file: " << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(in, line)) {
    auto comment = line.find('#');
    if (comment != std::string::npos)
      line.erase(comment);

    std::istringstream iss(line);
    std::string first;
    if (!(iss >> first))
      continue;
    if (first[0] == '[')
      continue;

    int device = std::atoi(first.c_str());
    int fifo = 0;
    int column = 0;
    int pixel = 0;
    double value = 0.;
    if (!(iss >> fifo >> column >> pixel >> value))
      continue;

    if (device != 200 || fifo < 0 || fifo > 3)
      continue;

    int eoch = pixel + 4 * column;
    if (eoch < 0 || eoch >= nchannels)
      continue;

    offset[eoch] = value;
    loaded[eoch] = true;
  }

  bool ok = true;
  for (int eoch = 0; eoch < nchannels; ++eoch) {
    if (!loaded[eoch]) {
      std::cerr << "ERROR: missing TIMING0 calibration for EOCH " << eoch << std::endl;
      ok = false;
    }
  }

  return ok;
}

} // namespace

void timing_reference_scan(const char *filename,
                           const char *calibfilename = "timing_offsets_from_cherenkov860.conf",
                           const char *outfilename = "timing_reference_scan.root",
                           int nbins = 2048,
                           double range = 32.,
                           double pitch = 1.,
                           double thickness = 2.)
{
  double offset[nchannels] = {};
  if (!load_timing0_offsets(calibfilename, offset))
    return;

  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto fout = TFile::Open(outfilename, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return;
  }

  TH2D *hAll[nchannels] = {};
  TH2D *hFirst[nchannels] = {};
  TH2D *hFirstDistance[nchannels] = {};
  TH2D *hFirstReflectedDistance[nchannels] = {};

  for (int ref = 0; ref < nchannels; ++ref) {
    hAll[ref] = new TH2D(TString::Format("hDeltaVsChannel_ref%02d", ref),
                         TString::Format("TIMING0;DO channel;time - DO %02d [clock]", ref),
                         nchannels, -0.5, nchannels - 0.5,
                         nbins, -range, range);

    hFirst[ref] = new TH2D(TString::Format("hDeltaVsChannel_first_ref%02d", ref),
                           TString::Format("TIMING0, DO %02d is first;DO channel;time - DO %02d [clock]",
                                           ref, ref),
                           nchannels, -0.5, nchannels - 0.5,
                           nbins, -range, range);

    hFirstDistance[ref] = new TH2D(TString::Format("hDeltaVsDistance_first_ref%02d", ref),
                                   TString::Format("TIMING0, DO %02d is first;distance from DO %02d;time - DO %02d [clock]",
                                                   ref, ref, ref),
                                   64, -0.05, 8.0 - 0.05,
                                   nbins, -range, range);

    hFirstReflectedDistance[ref] =
      new TH2D(TString::Format("hDeltaVsReflectedDistance_first_ref%02d", ref),
               TString::Format("TIMING0, DO %02d is first;reflected path length from DO %02d [cm];time - DO %02d [clock]",
                               ref, ref, ref),
               128, 0., 20.,
               nbins, -range, range);
  }

  Long64_t nframes = 0;
  Long64_t nframes_with_timing0 = 0;
  Long64_t nfills_all = 0;
  Long64_t nfills_first = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++nframes;

      bool have[nchannels] = {};
      double time[nchannels] = {};

      for (const auto &hit : reader.timing_hits()) {
        if (timing_detector(hit) != 0)
          continue;

        int eoch = electronics_channel(hit);
        if (eoch < 0)
          continue;

        int doch = eo2do[eoch];
        double calibrated_time = hit.time - offset[eoch];
        if (!have[doch] || calibrated_time < time[doch]) {
          have[doch] = true;
          time[doch] = calibrated_time;
        }
      }

      int first = -1;
      double first_time = std::numeric_limits<double>::max();
      for (int doch = 0; doch < nchannels; ++doch) {
        if (!have[doch])
          continue;
        if (time[doch] < first_time) {
          first = doch;
          first_time = time[doch];
        }
      }

      if (first < 0)
        continue;
      ++nframes_with_timing0;

      for (int ref = 0; ref < nchannels; ++ref) {
        if (!have[ref])
          continue;
        for (int doch = 0; doch < nchannels; ++doch) {
          if (!have[doch])
            continue;

          double deltat = time[doch] - time[ref];
          hAll[ref]->Fill(doch, deltat);
          ++nfills_all;

          if (first == ref) {
            hFirst[ref]->Fill(doch, deltat);
            hFirstDistance[ref]->Fill(do_distance(doch, ref), deltat);
            hFirstReflectedDistance[ref]->Fill(reflected_distance(doch, ref, pitch, thickness),
                                               deltat);
            ++nfills_first;
          }
        }
      }
    }
  }

  fout->cd();
  for (int ref = 0; ref < nchannels; ++ref) {
    hAll[ref]->Write();
    hFirst[ref]->Write();
    hFirstDistance[ref]->Write();
    hFirstReflectedDistance[ref]->Write();
  }
  fout->Close();

  std::cout << "frames processed:            " << nframes << std::endl;
  std::cout << "frames with TIMING0:         " << nframes_with_timing0 << std::endl;
  std::cout << "all-reference fills:         " << nfills_all << std::endl;
  std::cout << "first-reference fills:       " << nfills_first << std::endl;
  std::cout << "reflected path pitch [cm]:   " << pitch << std::endl;
  std::cout << "reflected path z [cm]:       " << thickness << std::endl;
  std::cout << "ROOT output:                 " << outfilename << std::endl;
}
