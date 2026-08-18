#include "../lib/trigger_reader.h"

#include <TFile.h>
#include <TH1D.h>
#include <TString.h>

#include <iostream>

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

} // namespace

void timing_pairs(const char *filename,
                  const char *outfilename = "timing_pairs.root",
                  int nbins = 2048,
                  double range = 32.)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto fout = TFile::Open(outfilename, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return;
  }

  TH1D *hDelta[2][nchannels][nchannels] = {};
  for (int timing = 0; timing < 2; ++timing) {
    for (int a = 0; a < nchannels; ++a) {
      for (int b = a + 1; b < nchannels; ++b) {
        auto name = TString::Format("deltat_%02d_%02d_timing%d", a, b, timing);
        auto title = TString::Format("TIMING%d DO %02d - DO %02d;#Deltat;entries",
                                     timing, a, b);
        hDelta[timing][a][b] = new TH1D(name.Data(), title.Data(), nbins, -range, range);
      }
    }
  }

  Long64_t nframes = 0;
  Long64_t nframes_with_timing[2] = {};
  Long64_t nfills[2] = {};

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++nframes;

      bool have[2][nchannels] = {};
      double time[2][nchannels] = {};

      for (const auto &hit : reader.timing_hits()) {
        int timing = timing_detector(hit);
        if (timing < 0)
          continue;

        int eoch = electronics_channel(hit);
        if (eoch < 0)
          continue;

        int doch = eo2do[eoch];
        if (!have[timing][doch] || hit.time < time[timing][doch]) {
          have[timing][doch] = true;
          time[timing][doch] = hit.time;
        }
      }

      for (int timing = 0; timing < 2; ++timing) {
        bool any = false;
        for (int a = 0; a < nchannels; ++a)
          any = any || have[timing][a];
        if (any)
          ++nframes_with_timing[timing];

        for (int a = 0; a < nchannels; ++a) {
          if (!have[timing][a])
            continue;
          for (int b = a + 1; b < nchannels; ++b) {
            if (!have[timing][b])
              continue;
            hDelta[timing][a][b]->Fill(time[timing][a] - time[timing][b]);
            ++nfills[timing];
          }
        }
      }
    }
  }

  for (int timing = 0; timing < 2; ++timing) {
    for (int a = 0; a < nchannels; ++a) {
      for (int b = a + 1; b < nchannels; ++b)
        hDelta[timing][a][b]->Write();
    }
  }
  fout->Close();

  std::cout << "frames processed:        " << nframes << std::endl;
  std::cout << "frames with TIMING0:     " << nframes_with_timing[0] << std::endl;
  std::cout << "frames with TIMING1:     " << nframes_with_timing[1] << std::endl;
  std::cout << "TIMING0 pair fills:      " << nfills[0] << std::endl;
  std::cout << "TIMING1 pair fills:      " << nfills[1] << std::endl;
  std::cout << "ROOT output:             " << outfilename << std::endl;
}
