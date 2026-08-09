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

struct frame_info_t {
  bool have[2][nchannels] = {};
  double time[2][nchannels] = {};
};

struct fit_event_t {
  std::vector<int> channel[2];
  std::vector<double> time[2];
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
      residual_sum[gch] += frame.time[0][channel] - mean0;
      ++residual_count[gch];
    }
    for (auto channel : channels1) {
      int gch = global_channel(1, channel);
      residual_sum[gch] += frame.time[1][channel] - mean1;
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
    offset[gch] += residual_sum[gch] / residual_count[gch];
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
    nfit = accumulate_residual_offsets(frames, offset,
                                       iter == pre_iterations - 1 ? &fit_frames : nullptr);
    if (nfit == 0) {
      std::cerr << "ERROR: no frames survived timing calibration selection" << std::endl;
      return false;
    }
    shift_timing1_to_zero_delta(fit_frames.empty() ? frames : fit_frames, offset);
  }

  if (fit_frames.empty()) {
    nfit = accumulate_residual_offsets(frames, offset, &fit_frames);
    if (nfit == 0) {
      std::cerr << "ERROR: no frames survived timing calibration selection" << std::endl;
      return false;
    }
  }

  auto fit_events = make_fit_events(fit_frames, offset);
  g_min_events = &fit_events;
  double objective_before = timing_objective(offset);
  bool minuit_ok = run_full_minimization(fit_events, offset, minimizer_step, minimizer_calls);
  double objective_after = timing_objective(offset);

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return false;
  }

  auto hDeltaBefore = new TH1D("hDeltaBefore", "", 400, -delta_range, delta_range);
  auto hDeltaAfter = new TH1D("hDeltaAfter", "", 400, -delta_range, delta_range);
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
  hOffsetValue->Write();
  hDeltaTiming0Before->Write();
  hDeltaTiming1Before->Write();
  hDeltaTiming0After->Write();
  hDeltaTiming1After->Write();
  fout->Close();

  write_channel_calibration(calibfilename.c_str(), offset);

  std::cout << "frames processed:               " << nframes << std::endl;
  std::cout << "frames with TIMING0 hits:        " << nwith0 << std::endl;
  std::cout << "frames with TIMING1 hits:        " << nwith1 << std::endl;
  std::cout << "frames retained:                 " << frames.size() << std::endl;
  std::cout << "frames used for calibration:     " << nfit << std::endl;
  std::cout << "events used by Minuit:           " << fit_events.size() << std::endl;
  std::cout << "pre-calibration iterations:      " << pre_iterations << std::endl;
  std::cout << "Minuit max calls:                " << minimizer_calls << std::endl;
  std::cout << "Minuit status:                   " << (minuit_ok ? "OK" : "WARNING") << std::endl;
  std::cout << "delta weight:                    " << g_delta_weight << std::endl;
  std::cout << "residual weight:                 " << g_residual_weight << std::endl;
  std::cout << "objective before Minuit:         " << objective_before << std::endl;
  std::cout << "objective after Minuit:          " << objective_after << std::endl;
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
  int pre_iterations = 5;
  double minimizer_step = 0.01;
  int minimizer_calls = 5000;
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
    ("minimizer-step", po::value<double>(&minimizer_step)->default_value(minimizer_step), "initial TMinuit parameter step")
    ("minimizer-calls", po::value<int>(&minimizer_calls)->default_value(minimizer_calls), "maximum TMinuit MIGRAD calls")
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
  if (minimizer_step <= 0.) {
    std::cerr << "ERROR: --minimizer-step must be positive" << std::endl;
    return 1;
  }
  if (minimizer_calls <= 0) {
    std::cerr << "ERROR: --minimizer-calls must be positive" << std::endl;
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
                      minimizer_step,
                      minimizer_calls,
                      delta_weight,
                      residual_weight,
                      max_frames) ? 0 : 1;
}
