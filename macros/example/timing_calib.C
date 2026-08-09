#include "../lib/trigger_reader.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

constexpr int nchannels = 32;
constexpr int ntiming = 64;
constexpr int reference_channel = 0;

struct frame_info_t {
  bool have[2][nchannels] = {};
  double time[2][nchannels] = {};
};

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

int global_channel(int det, int channel)
{
  return channel + nchannels * det;
}

double median(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  int n = values.size();
  if (n % 2 == 1)
    return values[n / 2];
  return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

bool selected_channels(const frame_info_t &frame,
                       int det,
                       const double *offset,
                       int min_channels,
                       double outlier_window,
                       std::vector<int> &channels,
                       double &mean)
{
  channels.clear();
  mean = 0.;

  std::vector<double> times;
  for (int channel = 0; channel < nchannels; ++channel) {
    if (!frame.have[det][channel])
      continue;
    int gch = global_channel(det, channel);
    times.push_back(frame.time[det][channel] - offset[gch]);
  }

  if ((int)times.size() < min_channels)
    return false;

  auto center = median(times);
  for (int channel = 0; channel < nchannels; ++channel) {
    if (!frame.have[det][channel])
      continue;
    int gch = global_channel(det, channel);
    auto time = frame.time[det][channel] - offset[gch];
    if (std::fabs(time - center) > outlier_window)
      continue;
    channels.push_back(channel);
    mean += time;
  }

  if ((int)channels.size() < min_channels)
    return false;

  mean /= channels.size();
  return true;
}

void write_channel_calibration(const char *filename, const double *offset)
{
  std::ofstream out(filename);
  if (!out) {
    std::cerr << "ERROR: could not create calibration output: " << filename << std::endl;
    return;
  }

  out << "[CHANNEL]\n";
  out << "# device fifo column pixel offset\n";
  out << "# TIMING channel offsets from timing_calib.C\n";
  out << "# calibrated_time = raw_time - offset\n";
  out << std::setprecision(12);

  for (int det = 0; det < 2; ++det) {
    int first_fifo = det == 0 ? 0 : 4;
    int last_fifo = det == 0 ? 3 : 7;
    for (int channel = 0; channel < nchannels; ++channel) {
      int column = channel / 4;
      int pixel = channel % 4;
      int gch = global_channel(det, channel);
      for (int fifo = first_fifo; fifo <= last_fifo; ++fifo) {
        out << 200 << ' '
            << fifo << ' '
            << column << ' '
            << pixel << ' '
            << offset[gch] << '\n';
      }
    }
  }
}

} // namespace

void timing_calib(const char *filename,
                  const char *outfilename = "timing_calib.root",
                  const char *calibfilename = "timing_channel_offsets.conf",
                  int min_channels = 16,
                  double outlier_window = 2.0,
                  double delta_range = 20.0,
                  double offset_range = 20.0)
{
  (void)offset_range;

  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  std::vector<frame_info_t> frames;
  int nframes = 0;
  int nwith0 = 0;
  int nwith1 = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++nframes;
      frame_info_t frame;

      for (const auto &hit : reader.timing_hits()) {
        int det = timing_detector(hit);
        if (det < 0)
          continue;
        int channel = timing_channel(hit);
        if (channel < 0)
          continue;

        if (!frame.have[det][channel] || hit.time < frame.time[det][channel]) {
          frame.time[det][channel] = hit.time;
          frame.have[det][channel] = true;
        }
      }

      bool have0 = false;
      bool have1 = false;
      for (int channel = 0; channel < nchannels; ++channel) {
        have0 = have0 || frame.have[0][channel];
        have1 = have1 || frame.have[1][channel];
      }
      if (have0) ++nwith0;
      if (have1) ++nwith1;

      frames.push_back(frame);
    }
  }

  double zero_offsets[ntiming] = {};
  double offset[ntiming] = {};
  double residual_sum[ntiming] = {};
  int residual_count[ntiming] = {};
  std::vector<int> selected0;
  std::vector<int> selected1;
  std::vector<frame_info_t> fit_frames;
  fit_frames.reserve(frames.size());

  int nfit = 0;
  for (const auto &frame : frames) {
    double mean0 = 0.;
    double mean1 = 0.;
    bool ok0 = selected_channels(frame, 0, zero_offsets, min_channels,
                                 outlier_window, selected0, mean0);
    auto channels0 = selected0;
    bool ok1 = selected_channels(frame, 1, zero_offsets, min_channels,
                                 outlier_window, selected1, mean1);
    auto channels1 = selected1;
    if (!ok0 || !ok1)
      continue;

    for (auto channel : channels0) {
      int gch = global_channel(0, channel);
      residual_sum[gch] += frame.time[0][channel] - mean0;
      ++residual_count[gch];
    }
    for (auto channel : channels1) {
      int gch = global_channel(1, channel);
      residual_sum[gch] += frame.time[1][channel] - mean1;
      ++residual_count[gch];
    }

    fit_frames.push_back(frame);
    ++nfit;
  }

  if (nfit == 0) {
    std::cerr << "ERROR: no frames survived timing calibration selection" << std::endl;
    return;
  }

  for (int gch = 0; gch < ntiming; ++gch) {
    if (residual_count[gch] == 0) {
      std::cerr << "WARNING: timing channel " << gch
                << " was never selected; offset kept at 0" << std::endl;
      continue;
    }
    offset[gch] = residual_sum[gch] / residual_count[gch];
  }

  double ref = offset[reference_channel];
  for (int gch = 0; gch < ntiming; ++gch)
    offset[gch] -= ref;
  offset[reference_channel] = 0.;

  double delta_sum = 0.;
  int delta_count = 0;
  for (const auto &frame : fit_frames) {
    double mean0 = 0.;
    double mean1 = 0.;
    bool ok0 = selected_channels(frame, 0, offset, min_channels,
                                 outlier_window, selected0, mean0);
    bool ok1 = selected_channels(frame, 1, offset, min_channels,
                                 outlier_window, selected1, mean1);
    if (!ok0 || !ok1)
      continue;
    delta_sum += mean0 - mean1;
    ++delta_count;
  }

  if (delta_count > 0) {
    double timing1_shift = -delta_sum / delta_count;
    for (int channel = 0; channel < nchannels; ++channel)
      offset[global_channel(1, channel)] += timing1_shift;
  }

  auto fout = TFile::Open(outfilename, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return;
  }

  auto hDeltaBefore = new TH1D("hDeltaBefore", "", 400, -delta_range, delta_range);
  auto hDeltaAfter = new TH1D("hDeltaAfter", "", 400, -delta_range, delta_range);
  auto hOffset = new TH1D("hOffset", "", ntiming, 0., ntiming);
  auto hDeltaTiming0Before = new TH2D("hDeltaTiming0Before", "", nchannels, 0., nchannels,
                                      400, -delta_range, delta_range);
  auto hDeltaTiming1Before = new TH2D("hDeltaTiming1Before", "", nchannels, 0., nchannels,
                                      400, -delta_range, delta_range);
  auto hDeltaTiming0After = new TH2D("hDeltaTiming0After", "", nchannels, 0., nchannels,
                                     400, -delta_range, delta_range);
  auto hDeltaTiming1After = new TH2D("hDeltaTiming1After", "", nchannels, 0., nchannels,
                                     400, -delta_range, delta_range);

  for (int gch = 0; gch < ntiming; ++gch)
    hOffset->SetBinContent(gch + 1, offset[gch]);

  int ndiag = 0;
  for (const auto &frame : fit_frames) {
    double mean0_before = 0.;
    double mean1_before = 0.;
    double mean0_after = 0.;
    double mean1_after = 0.;
    bool ok0_before = selected_channels(frame, 0, zero_offsets, min_channels,
                                        outlier_window, selected0, mean0_before);
    bool ok1_before = selected_channels(frame, 1, zero_offsets, min_channels,
                                        outlier_window, selected1, mean1_before);
    bool ok0_after = selected_channels(frame, 0, offset, min_channels,
                                       outlier_window, selected0, mean0_after);
    auto channels0_after = selected0;
    bool ok1_after = selected_channels(frame, 1, offset, min_channels,
                                       outlier_window, selected1, mean1_after);
    auto channels1_after = selected1;

    if (ok0_before && ok1_before)
      hDeltaBefore->Fill(mean0_before - mean1_before);
    if (ok0_after && ok1_after) {
      hDeltaAfter->Fill(mean0_after - mean1_after);
      ++ndiag;

      for (auto channel : channels0_after) {
        int gch = global_channel(0, channel);
        hDeltaTiming0After->Fill(channel, frame.time[0][channel] - offset[gch] - mean0_after);
      }
      for (auto channel : channels1_after) {
        int gch = global_channel(1, channel);
        hDeltaTiming1After->Fill(channel, frame.time[1][channel] - offset[gch] - mean1_after);
      }
    }

    double dummy_mean = 0.;
    if (ok0_before) {
      selected_channels(frame, 0, zero_offsets, min_channels,
                        outlier_window, selected0, dummy_mean);
      for (auto channel : selected0)
        hDeltaTiming0Before->Fill(channel, frame.time[0][channel] - mean0_before);
    }
    if (ok1_before) {
      selected_channels(frame, 1, zero_offsets, min_channels,
                        outlier_window, selected1, dummy_mean);
      for (auto channel : selected1)
        hDeltaTiming1Before->Fill(channel, frame.time[1][channel] - mean1_before);
    }
  }

  hDeltaBefore->Write();
  hDeltaAfter->Write();
  hOffset->Write();
  hDeltaTiming0Before->Write();
  hDeltaTiming1Before->Write();
  hDeltaTiming0After->Write();
  hDeltaTiming1After->Write();
  fout->Close();

  write_channel_calibration(calibfilename, offset);

  std::cout << "frames processed:               " << nframes << std::endl;
  std::cout << "frames with TIMING0 hits:        " << nwith0 << std::endl;
  std::cout << "frames with TIMING1 hits:        " << nwith1 << std::endl;
  std::cout << "frames used for calibration:     " << nfit << std::endl;
  std::cout << "frames used for diagnostics:     " << ndiag << std::endl;
  std::cout << "reference channel fixed:         fifo=0 column=0 pixel=0 offset=0" << std::endl;
  std::cout << "ROOT output:                     " << outfilename << std::endl;
  std::cout << "calibration output:              " << calibfilename << std::endl;
}
