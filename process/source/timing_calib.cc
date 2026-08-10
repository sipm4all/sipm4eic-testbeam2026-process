#include "../../macros/lib/trigger_reader.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMinuit.h>

#include <boost/program_options.hpp>

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
const int eo2do[nchannels] = {22, 20, 18, 16, 24, 26, 28, 30, 25, 27, 29, 31, 23, 21, 19, 17, 9, 11, 13, 15, 7, 5, 3, 1, 6, 4, 2, 0, 8, 10, 12, 14};

struct frame_info_t {
  bool have[2][nchannels] = {};
  double time[2][nchannels] = {};
};

struct fit_event_t {
  std::vector<int> channel[2];
  std::vector<double> time[2];
};

struct shape_event_t {
  double delta = 0.;
  double spread0 = 0.;
  double spread1 = 0.;
  double slope_x0 = 0.;
  double slope_x1 = 0.;
  double slope_y0 = 0.;
  double slope_y1 = 0.;
  double left_right0 = 0.;
  double left_right1 = 0.;
  double even_odd0 = 0.;
  double even_odd1 = 0.;
  double first_x0 = 0.;
  double first_x1 = 0.;
  double first_y0 = 0.;
  double first_y1 = 0.;
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


const std::vector<fit_event_t> *g_min_events = nullptr;
int g_min_channels = 16;
double g_outlier_window = 2.0;
double g_delta_weight = 1.0;
double g_residual_weight = 1.0;

int par_to_global(int ipar)
{
  return ipar < reference_channel ? ipar : ipar + 1;
}

void normalize_reference(double *offset)
{
  double ref = offset[reference_channel];
  for (int gch = 0; gch < ntiming; ++gch)
    offset[gch] -= ref;
  offset[reference_channel] = 0.;
}

int accumulate_residual_offsets(const std::vector<frame_info_t> &frames,
                                double *offset,
                                double damping,
                                std::vector<frame_info_t> *fit_frames = nullptr)
{
  double residual_sum[ntiming] = {};
  int residual_count[ntiming] = {};
  std::vector<int> selected0;
  std::vector<int> selected1;

  if (fit_frames)
    fit_frames->clear();

  int nfit = 0;
  for (const auto &frame : frames) {
    double mean0 = 0.;
    double mean1 = 0.;
    bool ok0 = selected_channels(frame, 0, offset, g_min_channels,
                                 g_outlier_window, selected0, mean0);
    auto channels0 = selected0;
    bool ok1 = selected_channels(frame, 1, offset, g_min_channels,
                                 g_outlier_window, selected1, mean1);
    auto channels1 = selected1;
    if (!ok0 || !ok1)
      continue;

    for (auto channel : channels0) {
      int gch = global_channel(0, channel);
      double time = frame.time[0][channel] - offset[gch];
      double reference = channels0.size() > 1 ?
        (mean0 * channels0.size() - time) / (channels0.size() - 1) : mean0;
      residual_sum[gch] += time - reference;
      ++residual_count[gch];
    }
    for (auto channel : channels1) {
      int gch = global_channel(1, channel);
      double time = frame.time[1][channel] - offset[gch];
      double reference = channels1.size() > 1 ?
        (mean1 * channels1.size() - time) / (channels1.size() - 1) : mean1;
      residual_sum[gch] += time - reference;
      ++residual_count[gch];
    }

    if (fit_frames)
      fit_frames->push_back(frame);
    ++nfit;
  }

  for (int gch = 0; gch < ntiming; ++gch) {
    if (residual_count[gch] == 0) {
      std::cerr << "WARNING: timing channel " << gch
                << " was never selected; offset unchanged" << std::endl;
      continue;
    }
    offset[gch] += damping * residual_sum[gch] / residual_count[gch];
  }

  normalize_reference(offset);
  return nfit;
}

int shift_timing1_to_zero_delta(const std::vector<frame_info_t> &frames,
                                double *offset)
{
  std::vector<int> selected0;
  std::vector<int> selected1;
  double delta_sum = 0.;
  int delta_count = 0;

  for (const auto &frame : frames) {
    double mean0 = 0.;
    double mean1 = 0.;
    bool ok0 = selected_channels(frame, 0, offset, g_min_channels,
                                 g_outlier_window, selected0, mean0);
    bool ok1 = selected_channels(frame, 1, offset, g_min_channels,
                                 g_outlier_window, selected1, mean1);
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

  normalize_reference(offset);
  return delta_count;
}

double timing_objective(const double *offset)
{
  double chi2 = 0.;
  int nterms = 0;

  for (const auto &event : *g_min_events) {
    if (event.channel[0].empty() || event.channel[1].empty())
      continue;

    double mean[2] = {};
    for (int det = 0; det < 2; ++det) {
      for (size_t i = 0; i < event.channel[det].size(); ++i)
        mean[det] += event.time[det][i] - offset[global_channel(det, event.channel[det][i])];
      mean[det] /= event.channel[det].size();
    }

    for (int det = 0; det < 2; ++det) {
      for (size_t i = 0; i < event.channel[det].size(); ++i) {
        double time = event.time[det][i] - offset[global_channel(det, event.channel[det][i])];
        double residual = time - mean[det];
        chi2 += g_residual_weight * residual * residual;
        ++nterms;
      }
    }

    double delta = mean[0] - mean[1];
    chi2 += g_delta_weight * delta * delta;
    ++nterms;
  }

  if (nterms == 0)
    return 1.e30;
  return chi2 / nterms;
}

void timing_fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  double offset[ntiming] = {};
  offset[reference_channel] = 0.;
  for (int ipar = 0; ipar < ntiming - 1; ++ipar)
    offset[par_to_global(ipar)] = par[ipar];
  f = timing_objective(offset);
}


std::vector<fit_event_t> make_fit_events(const std::vector<frame_info_t> &frames,
                                         const double *offset)
{
  std::vector<fit_event_t> events;
  events.reserve(frames.size());
  std::vector<int> selected;

  for (const auto &frame : frames) {
    fit_event_t event;
    bool ok = true;
    for (int det = 0; det < 2; ++det) {
      double mean = 0.;
      if (!selected_channels(frame, det, offset, g_min_channels,
                             g_outlier_window, selected, mean)) {
        ok = false;
        break;
      }
      for (auto channel : selected) {
        event.channel[det].push_back(channel);
        event.time[det].push_back(frame.time[det][channel]);
      }
    }
    if (ok)
      events.push_back(event);
  }

  return events;
}

bool run_full_minimization(const std::vector<fit_event_t> &events,
                           double *offset,
                           double step,
                           int max_calls)
{
  g_min_events = &events;

  TMinuit minuit(ntiming - 1);
  minuit.SetFCN(timing_fcn);

  Double_t arglist[2] = {-1., 0.};
  Int_t ierflg = 0;
  minuit.mnexcm("SET PRINT", arglist, 1, ierflg);

  for (int ipar = 0; ipar < ntiming - 1; ++ipar) {
    int gch = par_to_global(ipar);
    minuit.DefineParameter(ipar, Form("offset_%02d", gch), offset[gch], step, 0., 0.);
  }

  arglist[0] = max_calls;
  arglist[1] = 0.01;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

  for (int ipar = 0; ipar < ntiming - 1; ++ipar) {
    int gch = par_to_global(ipar);
    Double_t value = 0.;
    Double_t error = 0.;
    minuit.GetParameter(ipar, value, error);
    offset[gch] = value;
  }
  offset[reference_channel] = 0.;
  normalize_reference(offset);

  return ierflg == 0;
}



int first_selected_channel(const frame_info_t &frame,
                           int det,
                           const double *offset,
                           const std::vector<int> &channels)
{
  int first_channel = -1;
  double first_time = 0.;
  for (auto channel : channels) {
    int gch = global_channel(det, channel);
    double time = frame.time[det][channel] - offset[gch];
    if (first_channel < 0 || time < first_time) {
      first_channel = channel;
      first_time = time;
    }
  }
  return first_channel;
}

double first_selected_time(const frame_info_t &frame,
                           int det,
                           const double *offset,
                           const std::vector<int> &channels)
{
  double first = 0.;
  bool found = false;
  for (auto channel : channels) {
    int gch = global_channel(det, channel);
    double time = frame.time[det][channel] - offset[gch];
    if (!found || time < first) {
      first = time;
      found = true;
    }
  }
  return first;
}

bool event_shape(const frame_info_t &frame,
                 int det,
                 const double *offset,
                 const std::vector<int> &channels,
                 double mean,
                 double &spread,
                 double &slope_x,
                 double &slope_y,
                 double &left_right,
                 double &even_odd)
{
  if (channels.empty())
    return false;

  spread = 0.;
  slope_x = 0.;
  slope_y = 0.;
  left_right = 0.;
  even_odd = 0.;

  double sx = 0.;
  double sy = 0.;
  double sxx = 0.;
  double syy = 0.;
  double sxt = 0.;
  double syt = 0.;
  double left = 0.;
  double right = 0.;
  double even = 0.;
  double odd = 0.;
  int nleft = 0;
  int nright = 0;
  int neven = 0;
  int nodd = 0;

  for (auto channel : channels) {
    int gch = global_channel(det, channel);
    int detector_channel = eo2do[channel];
    double x = detector_channel % 4;
    double y = detector_channel / 4;
    double time = frame.time[det][channel] - offset[gch];
    double residual = time - mean;
    spread += residual * residual;
    sx += x;
    sy += y;
    sxx += x * x;
    syy += y * y;
    sxt += x * residual;
    syt += y * residual;

    if (x < 2) {
      left += residual;
      ++nleft;
    } else {
      right += residual;
      ++nright;
    }

    if (detector_channel % 2 == 0) {
      even += residual;
      ++neven;
    } else {
      odd += residual;
      ++nodd;
    }
  }

  spread = std::sqrt(spread / channels.size());
  double denx = channels.size() * sxx - sx * sx;
  double deny = channels.size() * syy - sy * sy;
  slope_x = denx != 0. ? channels.size() * sxt / denx : 0.;
  slope_y = deny != 0. ? channels.size() * syt / deny : 0.;
  left_right = nleft > 0 && nright > 0 ? left / nleft - right / nright : 0.;
  even_odd = neven > 0 && nodd > 0 ? even / neven - odd / nodd : 0.;
  return true;
}
double shape_value(const shape_event_t &event, int index)
{
  switch (index) {
  case 0: return 1.;
  case 1: return event.spread0;
  case 2: return event.spread1;
  case 3: return event.slope_x0;
  case 4: return event.slope_x1;
  case 5: return event.slope_y0;
  case 6: return event.slope_y1;
  case 7: return event.left_right0;
  case 8: return event.left_right1;
  case 9: return event.even_odd0;
  case 10: return event.even_odd1;
  case 11: return event.first_x0;
  case 12: return event.first_x1;
  case 13: return event.first_y0;
  case 14: return event.first_y1;
  default: return 0.;
  }
}

std::vector<double> fit_shape_model(const std::vector<shape_event_t> &events,
                                    int nvariables)
{
  std::vector<std::vector<double>> a(nvariables, std::vector<double>(nvariables + 1, 0.));

  for (const auto &event : events) {
    for (int i = 0; i < nvariables; ++i) {
      double xi = shape_value(event, i);
      for (int j = 0; j < nvariables; ++j)
        a[i][j] += xi * shape_value(event, j);
      a[i][nvariables] += xi * event.delta;
    }
  }

  for (int i = 0; i < nvariables; ++i) {
    int pivot = i;
    for (int k = i + 1; k < nvariables; ++k) {
      if (std::fabs(a[k][i]) > std::fabs(a[pivot][i]))
        pivot = k;
    }
    std::swap(a[i], a[pivot]);
    if (std::fabs(a[i][i]) < 1.e-20)
      continue;
    double div = a[i][i];
    for (int j = i; j <= nvariables; ++j)
      a[i][j] /= div;
    for (int k = 0; k < nvariables; ++k) {
      if (k == i)
        continue;
      double factor = a[k][i];
      for (int j = i; j <= nvariables; ++j)
        a[k][j] -= factor * a[i][j];
    }
  }

  std::vector<double> beta(nvariables, 0.);
  for (int i = 0; i < nvariables; ++i)
    beta[i] = a[i][nvariables];
  return beta;
}

double eval_shape_model(const shape_event_t &event,
                        const std::vector<double> &beta)
{
  double value = 0.;
  for (size_t i = 0; i < beta.size(); ++i)
    value += beta[i] * shape_value(event, i);
  return value;
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

bool timing_calib(const std::string &filename,
                  const std::string &outfilename,
                  const std::string &calibfilename,
                  int min_channels,
                  double outlier_window,
                  double delta_range,
                  double offset_range,
                  int pre_iterations,
                  double iteration_damping,
                  double minimizer_step,
                  int minimizer_calls,
                  double delta_weight,
                  double residual_weight,
                  int max_frames)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return false;

  std::vector<frame_info_t> frames;
  int nframes = 0;
  int nwith_trigger = 0;
  int nwithout_trigger = 0;
  int nwith0 = 0;
  int nwith1 = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++nframes;
      if (reader.trigger_hits().empty()) {
        ++nwithout_trigger;
        continue;
      }
      ++nwith_trigger;

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
      if (max_frames > 0 && (int)frames.size() >= max_frames)
        break;
    }
    if (max_frames > 0 && (int)frames.size() >= max_frames)
      break;
  }

  g_min_channels = min_channels;
  g_outlier_window = outlier_window;
  g_delta_weight = delta_weight;
  g_residual_weight = residual_weight;

  double zero_offsets[ntiming] = {};
  double offset[ntiming] = {};
  std::vector<frame_info_t> fit_frames;
  fit_frames.reserve(frames.size());

  int nfit = 0;
  for (int iter = 0; iter < pre_iterations; ++iter) {
    nfit = accumulate_residual_offsets(frames, offset, iteration_damping,
                                       iter == pre_iterations - 1 ? &fit_frames : nullptr);
    if (nfit == 0) {
      std::cerr << "ERROR: no frames survived timing calibration selection" << std::endl;
      return false;
    }
    shift_timing1_to_zero_delta(fit_frames.empty() ? frames : fit_frames, offset);
  }

  if (fit_frames.empty()) {
    nfit = accumulate_residual_offsets(frames, offset, iteration_damping, &fit_frames);
    if (nfit == 0) {
      std::cerr << "ERROR: no frames survived timing calibration selection" << std::endl;
      return false;
    }
  }

  auto fit_events = make_fit_events(fit_frames, offset);
  g_min_events = &fit_events;
  double objective_before = timing_objective(offset);
  bool minuit_ok = true;
  if (minimizer_calls > 0)
    minuit_ok = run_full_minimization(fit_events, offset, minimizer_step, minimizer_calls);
  double objective_after = timing_objective(offset);

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return false;
  }

  auto hDeltaBefore = new TH1D("hDeltaBefore", "", 400, -delta_range, delta_range);
  auto hDeltaAfter = new TH1D("hDeltaAfter", "", 400, -delta_range, delta_range);
  auto hDeltaFirstBefore = new TH1D("hDeltaFirstBefore", "", 400, -delta_range, delta_range);
  auto hDeltaFirstAfter = new TH1D("hDeltaFirstAfter", "", 400, -delta_range, delta_range);
  auto hDeltaShapeCorrected = new TH1D("hDeltaShapeCorrected", "", 400, -delta_range, delta_range);
  auto hDeltaShapePositionCorrected = new TH1D("hDeltaShapePositionCorrected", "", 400, -delta_range, delta_range);
  auto hDeltaVsSpread0 = new TH2D("hDeltaVsSpread0", "", 200, 0., delta_range, 400, -delta_range, delta_range);
  auto hDeltaVsSpread1 = new TH2D("hDeltaVsSpread1", "", 200, 0., delta_range, 400, -delta_range, delta_range);
  auto hTiming0RelativeVsSpread0 = new TH2D("hTiming0RelativeVsSpread0", "", 200, 0., delta_range, 400, -delta_range, delta_range);
  auto hTiming1RelativeVsSpread1 = new TH2D("hTiming1RelativeVsSpread1", "", 200, 0., delta_range, 400, -delta_range, delta_range);
  auto hDeltaVsSlopeX0 = new TH2D("hDeltaVsSlopeX0", "", 200, -0.2, 0.2, 400, -delta_range, delta_range);
  auto hDeltaVsSlopeX1 = new TH2D("hDeltaVsSlopeX1", "", 200, -0.2, 0.2, 400, -delta_range, delta_range);
  auto hDeltaVsSlopeY0 = new TH2D("hDeltaVsSlopeY0", "", 200, -0.2, 0.2, 400, -delta_range, delta_range);
  auto hDeltaVsSlopeY1 = new TH2D("hDeltaVsSlopeY1", "", 200, -0.2, 0.2, 400, -delta_range, delta_range);
  auto hDeltaVsLeftRight0 = new TH2D("hDeltaVsLeftRight0", "", 200, -delta_range, delta_range, 400, -delta_range, delta_range);
  auto hDeltaVsLeftRight1 = new TH2D("hDeltaVsLeftRight1", "", 200, -delta_range, delta_range, 400, -delta_range, delta_range);
  auto hDeltaVsFirstX0 = new TH2D("hDeltaVsFirstX0", "", 4, 0., 4., 400, -delta_range, delta_range);
  auto hDeltaVsFirstX1 = new TH2D("hDeltaVsFirstX1", "", 4, 0., 4., 400, -delta_range, delta_range);
  auto hDeltaVsFirstY0 = new TH2D("hDeltaVsFirstY0", "", 8, 0., 8., 400, -delta_range, delta_range);
  auto hDeltaVsFirstY1 = new TH2D("hDeltaVsFirstY1", "", 8, 0., 8., 400, -delta_range, delta_range);
  auto hFirstPosition0 = new TH2D("hFirstPosition0", "", 4, 0., 4., 8, 0., 8.);
  auto hFirstPosition1 = new TH2D("hFirstPosition1", "", 4, 0., 4., 8, 0., 8.);
  auto hFirstDetectorChannel0Vs1 = new TH2D("hFirstDetectorChannel0Vs1", "", nchannels, 0., nchannels, nchannels, 0., nchannels);
  auto hFirstElectronicsChannel0Vs1 = new TH2D("hFirstElectronicsChannel0Vs1", "", nchannels, 0., nchannels, nchannels, 0., nchannels);
  auto hFirstPositionDistance = new TH1D("hFirstPositionDistance", "", 100, 0., 9.);
  auto hDeltaVsFirstPositionDistance = new TH2D("hDeltaVsFirstPositionDistance", "", 100, 0., 9., 400, -delta_range, delta_range);
  auto hFirstPositionDxDy = new TH2D("hFirstPositionDxDy", "", 15, -7.5, 7.5, 15, -7.5, 7.5);
  auto hDeltaSameFirstDetectorChannel = new TH2D("hDeltaSameFirstDetectorChannel", "", nchannels, 0., nchannels,
                                                400, -delta_range, delta_range);
  std::array<TH1D *, nchannels> hDeltaSameFirstDetectorChannel1D{};
  for (int channel = 0; channel < nchannels; ++channel) {
    std::ostringstream name;
    name << "hDeltaSameFirstDetectorChannel";
    if (channel < 10)
      name << "0";
    name << channel;
    hDeltaSameFirstDetectorChannel1D[channel] = new TH1D(name.str().c_str(), "", 400, -delta_range, delta_range);
  }
  auto hTiming0SpreadBefore = new TH1D("hTiming0SpreadBefore", "", 400, 0., delta_range);
  auto hTiming1SpreadBefore = new TH1D("hTiming1SpreadBefore", "", 400, 0., delta_range);
  auto hTiming0SpreadAfter = new TH1D("hTiming0SpreadAfter", "", 400, 0., delta_range);
  auto hTiming1SpreadAfter = new TH1D("hTiming1SpreadAfter", "", 400, 0., delta_range);
  auto hExpectedDeltaFromSpreadBefore = new TH1D("hExpectedDeltaFromSpreadBefore", "", 400, 0., delta_range);
  auto hExpectedDeltaFromSpreadAfter = new TH1D("hExpectedDeltaFromSpreadAfter", "", 400, 0., delta_range);
  auto hOffset = new TH1D("hOffset", "", ntiming, 0., ntiming);
  auto hOffsetValue = new TH1D("hOffsetValue", "", 400, -offset_range, offset_range);
  auto hDeltaTiming0Before = new TH2D("hDeltaTiming0Before", "", nchannels, 0., nchannels,
                                      400, -delta_range, delta_range);
  auto hDeltaTiming1Before = new TH2D("hDeltaTiming1Before", "", nchannels, 0., nchannels,
                                      400, -delta_range, delta_range);
  auto hDeltaTiming0After = new TH2D("hDeltaTiming0After", "", nchannels, 0., nchannels,
                                     400, -delta_range, delta_range);
  auto hDeltaTiming1After = new TH2D("hDeltaTiming1After", "", nchannels, 0., nchannels,
                                     400, -delta_range, delta_range);

  for (int gch = 0; gch < ntiming; ++gch) {
    hOffset->SetBinContent(gch + 1, offset[gch]);
    hOffsetValue->Fill(offset[gch]);
  }

  std::vector<int> selected0;
  std::vector<int> selected1;

  auto selected_spread = [](const frame_info_t &frame,
                            int det,
                            const double *offset,
                            const std::vector<int> &channels,
                            double mean) {
    double s2 = 0.;
    for (auto channel : channels) {
      int gch = global_channel(det, channel);
      double residual = frame.time[det][channel] - offset[gch] - mean;
      s2 += residual * residual;
    }
    return channels.empty() ? 0. : std::sqrt(s2 / channels.size());
  };

  std::vector<shape_event_t> shape_events;
  shape_events.reserve(fit_frames.size());

  int ndiag = 0;
  for (const auto &frame : fit_frames) {
    double mean0_before = 0.;
    double mean1_before = 0.;
    double mean0_after = 0.;
    double mean1_after = 0.;
    bool ok0_before = selected_channels(frame, 0, zero_offsets, min_channels,
                                        outlier_window, selected0, mean0_before);
    auto channels0_before = selected0;
    bool ok1_before = selected_channels(frame, 1, zero_offsets, min_channels,
                                        outlier_window, selected1, mean1_before);
    auto channels1_before = selected1;
    bool ok0_after = selected_channels(frame, 0, offset, min_channels,
                                       outlier_window, selected0, mean0_after);
    auto channels0_after = selected0;
    bool ok1_after = selected_channels(frame, 1, offset, min_channels,
                                       outlier_window, selected1, mean1_after);
    auto channels1_after = selected1;

    if (ok0_before && ok1_before) {
      hDeltaBefore->Fill(mean0_before - mean1_before);
      hDeltaFirstBefore->Fill(first_selected_time(frame, 0, zero_offsets, channels0_before) -
                              first_selected_time(frame, 1, zero_offsets, channels1_before));
      double spread0 = selected_spread(frame, 0, zero_offsets, channels0_before, mean0_before);
      double spread1 = selected_spread(frame, 1, zero_offsets, channels1_before, mean1_before);
      hTiming0SpreadBefore->Fill(spread0);
      hTiming1SpreadBefore->Fill(spread1);
      hExpectedDeltaFromSpreadBefore->Fill(std::sqrt(spread0 * spread0 / channels0_before.size() +
                                                     spread1 * spread1 / channels1_before.size()));
    }
    if (ok0_after && ok1_after) {
      hDeltaAfter->Fill(mean0_after - mean1_after);
      hDeltaFirstAfter->Fill(first_selected_time(frame, 0, offset, channels0_after) -
                             first_selected_time(frame, 1, offset, channels1_after));
      double spread0 = selected_spread(frame, 0, offset, channels0_after, mean0_after);
      double spread1 = selected_spread(frame, 1, offset, channels1_after, mean1_after);
      hTiming0SpreadAfter->Fill(spread0);
      hTiming1SpreadAfter->Fill(spread1);
      hExpectedDeltaFromSpreadAfter->Fill(std::sqrt(spread0 * spread0 / channels0_after.size() +
                                                    spread1 * spread1 / channels1_after.size()));

      shape_event_t shape;
      shape.delta = mean0_after - mean1_after;
      shape.spread0 = spread0;
      shape.spread1 = spread1;
      event_shape(frame, 0, offset, channels0_after, mean0_after,
                  shape.spread0, shape.slope_x0, shape.slope_y0, shape.left_right0, shape.even_odd0);
      event_shape(frame, 1, offset, channels1_after, mean1_after,
                  shape.spread1, shape.slope_x1, shape.slope_y1, shape.left_right1, shape.even_odd1);
      int first0 = first_selected_channel(frame, 0, offset, channels0_after);
      int first1 = first_selected_channel(frame, 1, offset, channels1_after);
      int first_detector0 = -1;
      int first_detector1 = -1;
      if (first0 >= 0) {
        first_detector0 = eo2do[first0];
        shape.first_x0 = first_detector0 % 4;
        shape.first_y0 = first_detector0 / 4;
      }
      if (first1 >= 0) {
        first_detector1 = eo2do[first1];
        shape.first_x1 = first_detector1 % 4;
        shape.first_y1 = first_detector1 / 4;
      }
      shape_events.push_back(shape);
      hDeltaVsSpread0->Fill(shape.spread0, shape.delta);
      hDeltaVsSpread1->Fill(shape.spread1, shape.delta);
      double event_reference = 0.5 * (mean0_after + mean1_after);
      hTiming0RelativeVsSpread0->Fill(shape.spread0, mean0_after - event_reference);
      hTiming1RelativeVsSpread1->Fill(shape.spread1, mean1_after - event_reference);
      hDeltaVsSlopeX0->Fill(shape.slope_x0, shape.delta);
      hDeltaVsSlopeX1->Fill(shape.slope_x1, shape.delta);
      hDeltaVsSlopeY0->Fill(shape.slope_y0, shape.delta);
      hDeltaVsSlopeY1->Fill(shape.slope_y1, shape.delta);
      hDeltaVsLeftRight0->Fill(shape.left_right0, shape.delta);
      hDeltaVsLeftRight1->Fill(shape.left_right1, shape.delta);
      hDeltaVsFirstX0->Fill(shape.first_x0, shape.delta);
      hDeltaVsFirstX1->Fill(shape.first_x1, shape.delta);
      hDeltaVsFirstY0->Fill(shape.first_y0, shape.delta);
      hDeltaVsFirstY1->Fill(shape.first_y1, shape.delta);
      hFirstPosition0->Fill(shape.first_x0, shape.first_y0);
      hFirstPosition1->Fill(shape.first_x1, shape.first_y1);
      if (first_detector0 >= 0 && first_detector1 >= 0) {
        hFirstDetectorChannel0Vs1->Fill(first_detector0, first_detector1);
        double dx = shape.first_x0 - shape.first_x1;
        double dy = shape.first_y0 - shape.first_y1;
        double distance = std::sqrt(dx * dx + dy * dy);
        hFirstPositionDistance->Fill(distance);
        hDeltaVsFirstPositionDistance->Fill(distance, shape.delta);
        hFirstPositionDxDy->Fill(dx, dy);
        if (first_detector0 == first_detector1) {
          hDeltaSameFirstDetectorChannel->Fill(first_detector0, shape.delta);
          hDeltaSameFirstDetectorChannel1D[first_detector0]->Fill(shape.delta);
        }
      }
      if (first0 >= 0 && first1 >= 0)
        hFirstElectronicsChannel0Vs1->Fill(first0, first1);
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

  auto shape_beta = fit_shape_model(shape_events, 11);
  auto shape_position_beta = fit_shape_model(shape_events, 15);
  for (const auto &shape : shape_events) {
    hDeltaShapeCorrected->Fill(shape.delta - eval_shape_model(shape, shape_beta));
    hDeltaShapePositionCorrected->Fill(shape.delta - eval_shape_model(shape, shape_position_beta));
  }

  double timing0_spread_before = hTiming0SpreadBefore->GetMean();
  double timing0_spread_after = hTiming0SpreadAfter->GetMean();
  double timing1_spread_before = hTiming1SpreadBefore->GetMean();
  double timing1_spread_after = hTiming1SpreadAfter->GetMean();
  double expected_delta_before = hExpectedDeltaFromSpreadBefore->GetMean();
  double expected_delta_after = hExpectedDeltaFromSpreadAfter->GetMean();
  double observed_delta_rms_before = hDeltaBefore->GetRMS();
  double observed_delta_rms_after = hDeltaAfter->GetRMS();
  double observed_delta_first_rms_before = hDeltaFirstBefore->GetRMS();
  double observed_delta_first_rms_after = hDeltaFirstAfter->GetRMS();
  double shape_corrected_delta_rms = hDeltaShapeCorrected->GetRMS();
  double shape_position_corrected_delta_rms = hDeltaShapePositionCorrected->GetRMS();

  hDeltaBefore->Write();
  hDeltaAfter->Write();
  hDeltaFirstBefore->Write();
  hDeltaFirstAfter->Write();
  hDeltaShapeCorrected->Write();
  hDeltaShapePositionCorrected->Write();
  hDeltaVsSpread0->Write();
  hDeltaVsSpread1->Write();
  hTiming0RelativeVsSpread0->Write();
  hTiming1RelativeVsSpread1->Write();
  hDeltaVsSlopeX0->Write();
  hDeltaVsSlopeX1->Write();
  hDeltaVsSlopeY0->Write();
  hDeltaVsSlopeY1->Write();
  hDeltaVsLeftRight0->Write();
  hDeltaVsLeftRight1->Write();
  hDeltaVsFirstX0->Write();
  hDeltaVsFirstX1->Write();
  hDeltaVsFirstY0->Write();
  hDeltaVsFirstY1->Write();
  hFirstPosition0->Write();
  hFirstPosition1->Write();
  hFirstDetectorChannel0Vs1->Write();
  hFirstElectronicsChannel0Vs1->Write();
  hFirstPositionDistance->Write();
  hDeltaVsFirstPositionDistance->Write();
  hFirstPositionDxDy->Write();
  hDeltaSameFirstDetectorChannel->Write();
  for (auto hist : hDeltaSameFirstDetectorChannel1D)
    hist->Write();
  hTiming0SpreadBefore->Write();
  hTiming1SpreadBefore->Write();
  hTiming0SpreadAfter->Write();
  hTiming1SpreadAfter->Write();
  hExpectedDeltaFromSpreadBefore->Write();
  hExpectedDeltaFromSpreadAfter->Write();
  hOffset->Write();
  hOffsetValue->Write();
  hDeltaTiming0Before->Write();
  hDeltaTiming1Before->Write();
  hDeltaTiming0After->Write();
  hDeltaTiming1After->Write();
  fout->Close();

  write_channel_calibration(calibfilename.c_str(), offset);

  std::cout << "frames processed:               " << nframes << std::endl;
  std::cout << "frames with trigger hits:        " << nwith_trigger << std::endl;
  std::cout << "frames without trigger hits:     " << nwithout_trigger << std::endl;
  std::cout << "frames with TIMING0 hits:        " << nwith0 << std::endl;
  std::cout << "frames with TIMING1 hits:        " << nwith1 << std::endl;
  std::cout << "frames retained:                 " << frames.size() << std::endl;
  std::cout << "frames used for calibration:     " << nfit << std::endl;
  std::cout << "events used by Minuit:           " << fit_events.size() << std::endl;
  std::cout << "pre-calibration iterations:      " << pre_iterations << std::endl;
  std::cout << "iteration damping:               " << iteration_damping << std::endl;
  std::cout << "Minuit max calls:                " << minimizer_calls << std::endl;
  std::cout << "Minuit status:                   " << (minuit_ok ? "OK" : "WARNING") << std::endl;
  std::cout << "delta weight:                    " << g_delta_weight << std::endl;
  std::cout << "residual weight:                 " << g_residual_weight << std::endl;
  std::cout << "objective before Minuit:         " << objective_before << std::endl;
  std::cout << "objective after Minuit:          " << objective_after << std::endl;
  std::cout << "TIMING0 spread before/after:     " << timing0_spread_before
            << " / " << timing0_spread_after << std::endl;
  std::cout << "TIMING1 spread before/after:     " << timing1_spread_before
            << " / " << timing1_spread_after << std::endl;
  std::cout << "expected delta spread before/after: " << expected_delta_before
            << " / " << expected_delta_after << std::endl;
  std::cout << "observed delta RMS before/after: " << observed_delta_rms_before
            << " / " << observed_delta_rms_after << std::endl;
  std::cout << "first-hit delta RMS before/after: " << observed_delta_first_rms_before
            << " / " << observed_delta_first_rms_after << std::endl;
  std::cout << "shape-corrected delta RMS:       " << shape_corrected_delta_rms << std::endl;
  std::cout << "shape+first-position delta RMS:  " << shape_position_corrected_delta_rms << std::endl;
  std::cout << "frames used for diagnostics:     " << ndiag << std::endl;
  std::cout << "reference channel fixed:         fifo=0 column=0 pixel=0 offset=0" << std::endl;
  std::cout << "ROOT output:                     " << outfilename << std::endl;
  std::cout << "calibration output:              " << calibfilename << std::endl;
  return true;
}

int main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string output = "timing_calib.root";
  std::string calibration_output = "timing_channel_offsets.conf";
  int min_channels = 16;
  double outlier_window = 2.0;
  double delta_range = 20.0;
  double offset_range = 20.0;
  int pre_iterations = 3;
  double iteration_damping = 0.5;
  double minimizer_step = 0.01;
  int minimizer_calls = 0;
  double delta_weight = 1.0;
  double residual_weight = 1.0;
  int max_frames = 0;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input triggered-frame ROOT file")
    ("output,o", po::value<std::string>(&output)->default_value(output), "output diagnostic ROOT file")
    ("calibration-output,c", po::value<std::string>(&calibration_output)->default_value(calibration_output), "output [CHANNEL] calibration snippet")
    ("min-channels", po::value<int>(&min_channels)->default_value(min_channels), "minimum selected channels per timing scintillator")
    ("outlier-window", po::value<double>(&outlier_window)->default_value(outlier_window), "robust mean outlier window")
    ("delta-range", po::value<double>(&delta_range)->default_value(delta_range), "delta/residual histogram half range")
    ("offset-range", po::value<double>(&offset_range)->default_value(offset_range), "offset-value histogram half range")
    ("pre-iterations", po::value<int>(&pre_iterations)->default_value(pre_iterations), "iterative residual pre-calibration passes")
    ("iteration-damping", po::value<double>(&iteration_damping)->default_value(iteration_damping), "damping factor for each iterative residual update")
    ("minimizer-step", po::value<double>(&minimizer_step)->default_value(minimizer_step), "initial TMinuit parameter step")
    ("minimizer-calls", po::value<int>(&minimizer_calls)->default_value(minimizer_calls), "maximum TMinuit MIGRAD calls, 0 skips Minuit")
    ("delta-weight", po::value<double>(&delta_weight)->default_value(delta_weight), "weight for TIMING0_mean - TIMING1_mean term")
    ("residual-weight", po::value<double>(&residual_weight)->default_value(residual_weight), "weight for intra-scintillator channel residual terms")
    ("max-frames", po::value<int>(&max_frames)->default_value(max_frames), "maximum frames to read, 0 means all frames")
    ;

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

  if (min_channels <= 0 || min_channels > nchannels) {
    std::cerr << "ERROR: --min-channels must be in [1," << nchannels << "]" << std::endl;
    return 1;
  }
  if (pre_iterations < 0) {
    std::cerr << "ERROR: --pre-iterations must be non-negative" << std::endl;
    return 1;
  }
  if (iteration_damping <= 0. || iteration_damping > 1.) {
    std::cerr << "ERROR: --iteration-damping must be in (0,1]" << std::endl;
    return 1;
  }
  if (minimizer_step <= 0.) {
    std::cerr << "ERROR: --minimizer-step must be positive" << std::endl;
    return 1;
  }
  if (minimizer_calls < 0) {
    std::cerr << "ERROR: --minimizer-calls must be non-negative" << std::endl;
    return 1;
  }
  if (delta_weight < 0. || residual_weight < 0. || delta_weight + residual_weight <= 0.) {
    std::cerr << "ERROR: --delta-weight and --residual-weight must be non-negative, and at least one must be positive" << std::endl;
    return 1;
  }
  if (max_frames < 0) {
    std::cerr << "ERROR: --max-frames must be non-negative" << std::endl;
    return 1;
  }

  return timing_calib(input,
                      output,
                      calibration_output,
                      min_channels,
                      outlier_window,
                      delta_range,
                      offset_range,
                      pre_iterations,
                      iteration_damping,
                      minimizer_step,
                      minimizer_calls,
                      delta_weight,
                      residual_weight,
                      max_frames) ? 0 : 1;
}
