#include "../lib/trigger_reader.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

constexpr int nchannels = 32;

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

int timing_channel(const hit_t &hit)
{
  int channel = hit.pixel + 4 * hit.column;
  if (channel < 0 || channel >= nchannels)
    return -1;
  return channel;
}

double median(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  int n = values.size();
  if (n % 2 == 1)
    return values[n / 2];
  return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

bool robust_mean(const std::vector<double> &times,
                 int min_channels,
                 double outlier_window,
                 double &mean,
                 int &nused)
{
  nused = 0;
  mean = 0.;

  if ((int)times.size() < min_channels)
    return false;

  auto center = median(times);
  for (auto time : times) {
    if (std::fabs(time - center) > outlier_window)
      continue;
    mean += time;
    ++nused;
  }

  if (nused < min_channels)
    return false;

  mean /= nused;
  return true;
}

} // namespace

void timing(const char *filename,
            const char *outfilename = "timing.root",
            int min_channels = 16,
            double outlier_window = 2.0,
            double residual_range = 10.0,
            double delta_range = 20.0)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  auto fout = TFile::Open(outfilename, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return;
  }

  auto hDeltaTiming0 = new TH2D("hDeltaTiming0", "", nchannels, 0., nchannels,
                                400, -residual_range, residual_range);
  auto hDeltaTiming1 = new TH2D("hDeltaTiming1", "", nchannels, 0., nchannels,
                                400, -residual_range, residual_range);
  auto hDelta = new TH1D("hDelta", "", 400, -delta_range, delta_range);

  int nframes = 0;
  int nwith0 = 0;
  int nwith1 = 0;
  int naccepted = 0;
  int nreject0 = 0;
  int nreject1 = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++nframes;

      std::vector<double> channel_time[2];
      bool have_channel[2][nchannels] = {};
      double earliest[2][nchannels] = {};

      for (const auto &hit : reader.timing_hits()) {
        int det = timing_detector(hit);
        if (det < 0)
          continue;
        int channel = timing_channel(hit);
        if (channel < 0)
          continue;

        if (!have_channel[det][channel] || hit.time < earliest[det][channel]) {
          earliest[det][channel] = hit.time;
          have_channel[det][channel] = true;
        }
      }

      for (int det = 0; det < 2; ++det) {
        for (int channel = 0; channel < nchannels; ++channel) {
          if (have_channel[det][channel])
            channel_time[det].push_back(earliest[det][channel]);
        }
      }

      if (!channel_time[0].empty()) ++nwith0;
      if (!channel_time[1].empty()) ++nwith1;

      double mean0 = 0.;
      double mean1 = 0.;
      int nused0 = 0;
      int nused1 = 0;
      bool ok0 = robust_mean(channel_time[0], min_channels, outlier_window, mean0, nused0);
      bool ok1 = robust_mean(channel_time[1], min_channels, outlier_window, mean1, nused1);

      if (!ok0) ++nreject0;
      if (!ok1) ++nreject1;
      if (!ok0 || !ok1)
        continue;

      ++naccepted;
      hDelta->Fill(mean0 - mean1);

      for (int channel = 0; channel < nchannels; ++channel) {
        if (have_channel[0][channel] && std::fabs(earliest[0][channel] - mean0) <= outlier_window)
          hDeltaTiming0->Fill(channel, earliest[0][channel] - mean0);
        if (have_channel[1][channel] && std::fabs(earliest[1][channel] - mean1) <= outlier_window)
          hDeltaTiming1->Fill(channel, earliest[1][channel] - mean1);
      }
    }
  }

  std::cout << "frames processed:               " << nframes << std::endl;
  std::cout << "frames with TIMING0 hits:        " << nwith0 << std::endl;
  std::cout << "frames with TIMING1 hits:        " << nwith1 << std::endl;
  std::cout << "frames accepted:                 " << naccepted << std::endl;
  std::cout << "frames rejected by TIMING0:      " << nreject0 << std::endl;
  std::cout << "frames rejected by TIMING1:      " << nreject1 << std::endl;

  hDeltaTiming0->Write();
  hDeltaTiming1->Write();
  hDelta->Write();
  fout->Close();
}
