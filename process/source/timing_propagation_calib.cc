#include "../../macros/lib/trigger_reader.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMinuit.h>

#include <boost/program_options.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int nchannels = 32;
constexpr int pitch_parameter = 0;
constexpr int depth_parameter = 1;
constexpr int slope_parameter = 2;
constexpr int power_parameter = 3;
constexpr int nparameters = 4;

const int eo2do[nchannels] = {22, 20, 18, 16, 24, 26, 28, 30,
                              25, 27, 29, 31, 23, 21, 19, 17,
                              9,  11, 13, 15, 7,  5,  3,  1,
                              6,  4,  2,  0,  8,  10, 12, 14};

struct event_t {
  Long64_t frame = -1;
  std::vector<int> channel;
  std::vector<double> time;
};

struct event_fit_t {
  bool ok = false;
  double x = 0.;
  double y = 0.;
  double t0 = 0.;
  double sum2 = 0.;
  int n = 0;
};

struct fit_context_t {
  int timing = 0;
  const double *offset = nullptr;
  int min_channels = 8;
  int grid_steps = 9;
  std::vector<event_t> events;
  long long residuals = 0;
};

fit_context_t *g_context = nullptr;

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

double channel_x(int channel, double pitch)
{
  return pitch * (channel % 4);
}

double channel_y(int channel, double pitch)
{
  return pitch * (channel / 4);
}

double propagation_path(double x, double y, int channel, double pitch, double depth)
{
  double dx = channel_x(channel, pitch) - x;
  double dy = channel_y(channel, pitch) - y;
  return std::sqrt(dx * dx + dy * dy + depth * depth);
}

double propagation(double x,
                   double y,
                   int channel,
                   double pitch,
                   double depth,
                   double slope,
                   double power)
{
  return slope * std::pow(propagation_path(x, y, channel, pitch, depth), power);
}

double event_objective(const event_t &event,
                       const double offset[nchannels],
                       double pitch,
                       double depth,
                       double slope,
                       double power,
                       double x,
                       double y,
                       double &t0)
{
  double sum = 0.;
  int n = 0;
  for (size_t i = 0; i < event.channel.size(); ++i) {
    int channel = event.channel[i];
    sum += event.time[i] - offset[channel] -
           propagation(x, y, channel, pitch, depth, slope, power);
    ++n;
  }

  if (n == 0) {
    t0 = 0.;
    return std::numeric_limits<double>::max();
  }

  t0 = sum / n;

  double sum2 = 0.;
  for (size_t i = 0; i < event.channel.size(); ++i) {
    int channel = event.channel[i];
    double expected = t0 + offset[channel] +
                      propagation(x, y, channel, pitch, depth, slope, power);
    double residual = event.time[i] - expected;
    sum2 += residual * residual;
  }

  return sum2;
}

event_fit_t fit_event_position(const event_t &event,
                               const double offset[nchannels],
                               double pitch,
                               double depth,
                               double slope,
                               double power,
                               int min_channels,
                               int grid_steps)
{
  event_fit_t fit;
  fit.n = event.channel.size();
  if ((int)event.channel.size() < min_channels)
    return fit;

  double xmin = -0.5 * pitch;
  double xmax = 3.5 * pitch;
  double ymin = -0.5 * pitch;
  double ymax = 7.5 * pitch;

  double best = std::numeric_limits<double>::max();
  double best_x = 0.;
  double best_y = 0.;
  double best_t0 = 0.;

  auto scan = [&](double cx, double cy, double half_range) {
    int nstep = std::max(3, grid_steps);
    for (int ix = 0; ix < nstep; ++ix) {
      double x = cx - half_range + 2. * half_range * ix / std::max(1, nstep - 1);
      x = std::clamp(x, xmin, xmax);
      for (int iy = 0; iy < nstep; ++iy) {
        double y = cy - half_range + 2. * half_range * iy / std::max(1, nstep - 1);
        y = std::clamp(y, ymin, ymax);
        double t0 = 0.;
        double value = event_objective(event, offset, pitch, depth, slope, power, x, y, t0);
        if (value < best) {
          best = value;
          best_x = x;
          best_y = y;
          best_t0 = t0;
        }
      }
    }
  };

  scan(0.5 * (xmin + xmax), 0.5 * (ymin + ymax), std::max(xmax - xmin, ymax - ymin));
  scan(best_x, best_y, 1.5 * pitch);
  scan(best_x, best_y, 0.35 * pitch);
  scan(best_x, best_y, 0.08 * pitch);

  fit.ok = true;
  fit.x = best_x;
  fit.y = best_y;
  fit.t0 = best_t0;
  fit.sum2 = best;
  return fit;
}

bool load_offsets(const std::string &filename, double offset[2][nchannels])
{
  bool loaded[2][nchannels] = {};
  for (int timing = 0; timing < 2; ++timing)
    for (int channel = 0; channel < nchannels; ++channel)
      offset[timing][channel] = 0.;

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

    if (device != 200)
      continue;

    int timing = -1;
    if (fifo >= 0 && fifo <= 3)
      timing = 0;
    if (fifo >= 4 && fifo <= 7)
      timing = 1;
    if (timing < 0)
      continue;

    int eoch = pixel + 4 * column;
    if (eoch < 0 || eoch >= nchannels)
      continue;

    int doch = eo2do[eoch];
    offset[timing][doch] = value;
    loaded[timing][doch] = true;
  }

  bool ok = true;
  for (int timing = 0; timing < 2; ++timing) {
    for (int channel = 0; channel < nchannels; ++channel) {
      if (!loaded[timing][channel]) {
        std::cerr << "ERROR: missing initial calibration for TIMING" << timing
                  << " DO " << channel << std::endl;
        ok = false;
      }
    }
  }

  return ok;
}

void fcn(int &, double *, double &value, double *par, int)
{
  const auto &events = g_context->events;

  double pitch = par[pitch_parameter];
  double depth = par[depth_parameter];
  double slope = par[slope_parameter];
  double power = par[power_parameter];

  if (pitch <= 0. || depth <= 0. || power <= 0.) {
    value = 1.e30;
    return;
  }

  double sum2 = 0.;
  long long n = 0;

  for (const auto &event : events) {
    auto fit = fit_event_position(event,
                                  g_context->offset,
                                  pitch,
                                  depth,
                                  slope,
                                  power,
                                  g_context->min_channels,
                                  g_context->grid_steps);
    if (!fit.ok)
      continue;

    sum2 += fit.sum2;
    n += fit.n;
  }

  value = n > 0 ? sum2 / n : 1.e30;
}

std::vector<event_t> collect_events(const std::string &filename,
                                    int timing,
                                    int min_channels,
                                    int max_events,
                                    bool require_trigger)
{
  trigger_reader_t reader;
  if (!reader.open(filename))
    return {};

  std::vector<event_t> events;
  Long64_t frame = 0;
  while (reader.next_spill()) {
    while (reader.next_frame()) {
      if (require_trigger && reader.trigger_hits().empty()) {
        ++frame;
        continue;
      }

      bool have[nchannels] = {};
      double time[nchannels] = {};

      for (const auto &hit : reader.timing_hits()) {
        if (timing_detector(hit) != timing)
          continue;

        int eoch = electronics_channel(hit);
        if (eoch < 0)
          continue;

        int doch = eo2do[eoch];
        if (!have[doch] || hit.time < time[doch]) {
          have[doch] = true;
          time[doch] = hit.time;
        }
      }

      event_t event;
      event.frame = frame++;
      for (int channel = 0; channel < nchannels; ++channel) {
        if (!have[channel])
          continue;
        event.channel.push_back(channel);
        event.time.push_back(time[channel]);
      }

      if ((int)event.channel.size() < min_channels)
        continue;

      events.push_back(std::move(event));
      if (max_events > 0 && (int)events.size() >= max_events)
        return events;
    }
  }

  return events;
}

bool fit_timing(int timing,
                const std::vector<event_t> &events,
                const double initial_offset[nchannels],
                double initial_pitch,
                double initial_depth,
                double initial_slope,
                double initial_power,
                int max_calls,
                int print_level,
                int min_channels,
                int grid_steps,
                double result[nparameters])
{
  fit_context_t context;
  context.timing = timing;
  context.offset = initial_offset;
  context.min_channels = min_channels;
  context.grid_steps = grid_steps;
  context.events = events;
  for (const auto &event : events)
    context.residuals += (long long)event.channel.size();
  g_context = &context;

  TMinuit minuit(nparameters);
  minuit.SetFCN(fcn);
  minuit.SetPrintLevel(print_level);

  int ierr = 0;
  minuit.mnparm(pitch_parameter, "pitch_cm", initial_pitch, 0.01, 0.05, 2.0, ierr);
  minuit.mnparm(depth_parameter, "depth_cm", initial_depth, 0.01, 0.01, 10.0, ierr);
  minuit.mnparm(slope_parameter, "slope_clock_per_cm", initial_slope, 0.001, -5.0, 5.0, ierr);
  minuit.mnparm(power_parameter, "power", initial_power, 0.01, 0.1, 5.0, ierr);

  double arglist[2] = {(double)max_calls, 0.1};
  minuit.mnexcm("MIGRAD", arglist, 2, ierr);
  if (ierr != 0)
    std::cerr << "WARNING: MIGRAD returned status " << ierr
              << " for TIMING" << timing << std::endl;

  for (int ipar = 0; ipar < nparameters; ++ipar) {
    double error = 0.;
    minuit.GetParameter(ipar, result[ipar], error);
  }

  std::cout << "TIMING" << timing
            << " events=" << events.size()
            << " residuals=" << context.residuals
            << " pitch=" << result[pitch_parameter]
            << " depth=" << result[depth_parameter]
            << " slope=" << result[slope_parameter]
            << " power=" << result[power_parameter]
            << std::endl;

  return true;
}

void write_outputs(const std::string &outroot,
                   const std::string &outconf,
                   const double initial_offset[2][nchannels],
                   double initial_pitch,
                   double initial_depth,
                   double initial_slope,
                   double initial_power,
                   const double result[2][nparameters],
                   const std::vector<event_t> events[2])
{
  TFile fout(outroot.c_str(), "RECREATE");

  for (int timing = 0; timing < 2; ++timing) {
    auto hOffset = new TH1D(TString::Format("hOffset_timing%d", timing),
                           ";DO channel;fixed offset [clock]",
                           nchannels, -0.5, nchannels - 0.5);
    auto hResidualInitial = new TH1D(TString::Format("hResidualInitial_timing%d", timing),
                                    ";measured - model [clock];entries",
                                    2048, -32., 32.);
    auto hResidualFit = new TH1D(TString::Format("hResidualFit_timing%d", timing),
                                ";measured - model [clock];entries",
                                2048, -32., 32.);
    auto hResidualVsPath = new TH2D(TString::Format("hResidualVsPath_timing%d", timing),
                                   ";direct path [cm];residual [clock]",
                                   200, 0., 5., 1024, -8., 8.);
    auto hYvsX = new TH2D(TString::Format("hYvsX_timing%d", timing),
                          ";fitted x [cm];fitted y [cm]",
                          120, -0.5 * result[timing][pitch_parameter],
                          3.5 * result[timing][pitch_parameter],
                          200, -0.5 * result[timing][pitch_parameter],
                          7.5 * result[timing][pitch_parameter]);

    for (int channel = 0; channel < nchannels; ++channel) {
      hOffset->SetBinContent(channel + 1, initial_offset[timing][channel]);
    }

    double initial_par[nparameters] = {};
    initial_par[pitch_parameter] = initial_pitch;
    initial_par[depth_parameter] = initial_depth;
    initial_par[slope_parameter] = initial_slope;
    initial_par[power_parameter] = initial_power;

    auto fill_event_residuals = [&](const event_t &event, const double *par,
                                    TH1D *hist, bool fill_diagnostics) {
      auto fit = fit_event_position(event,
                                    initial_offset[timing],
                                    par[pitch_parameter],
                                    par[depth_parameter],
                                    par[slope_parameter],
                                    par[power_parameter],
                                    2,
                                    11);
      if (!fit.ok)
        return;
      if (fill_diagnostics)
        hYvsX->Fill(fit.x, fit.y);

      for (size_t i = 0; i < event.channel.size(); ++i) {
        int channel = event.channel[i];
        double path = propagation_path(fit.x, fit.y, channel,
                                     par[pitch_parameter],
                                     par[depth_parameter]);
        double expected = fit.t0 + initial_offset[timing][channel] +
                          par[slope_parameter] * std::pow(path, par[power_parameter]);
        double residual = event.time[i] - expected;
        hist->Fill(residual);
        if (fill_diagnostics)
          hResidualVsPath->Fill(path, residual);
      }
    };

    for (const auto &event : events[timing]) {
      fill_event_residuals(event, initial_par, hResidualInitial, false);
      fill_event_residuals(event, result[timing], hResidualFit, true);
    }

    hOffset->Write();
    hResidualInitial->Write();
    hResidualFit->Write();
    hResidualVsPath->Write();
    hYvsX->Write();
  }

  auto hDeltaTimingEventT0 =
    new TH1D("hDeltaTimingEventT0",
             ";fitted t0 TIMING0 - fitted t0 TIMING1 [clock];events",
             2048, -32., 32.);
  auto hX0vsX1 =
    new TH2D("hX0vsX1",
             ";TIMING0 fitted x [cm];TIMING1 fitted x [cm]",
             120, -0.5 * result[0][pitch_parameter],
             3.5 * result[0][pitch_parameter],
             120, -0.5 * result[1][pitch_parameter],
             3.5 * result[1][pitch_parameter]);
  auto hY0vsY1 =
    new TH2D("hY0vsY1",
             ";TIMING0 fitted y [cm];TIMING1 fitted y [cm]",
             200, -0.5 * result[0][pitch_parameter],
             7.5 * result[0][pitch_parameter],
             200, -0.5 * result[1][pitch_parameter],
             7.5 * result[1][pitch_parameter]);
  auto hPositionDistance =
    new TH1D("hPositionDistance",
             ";distance between fitted TIMING0 and TIMING1 positions [cm];events",
             1000, 0., 10.);

  size_t i0 = 0;
  size_t i1 = 0;
  while (i0 < events[0].size() && i1 < events[1].size()) {
    if (events[0][i0].frame < events[1][i1].frame) {
      ++i0;
      continue;
    }
    if (events[1][i1].frame < events[0][i0].frame) {
      ++i1;
      continue;
    }

    auto fit0 = fit_event_position(events[0][i0],
                                   initial_offset[0],
                                   result[0][pitch_parameter],
                                   result[0][depth_parameter],
                                   result[0][slope_parameter],
                                   result[0][power_parameter],
                                   2,
                                   11);
    auto fit1 = fit_event_position(events[1][i1],
                                   initial_offset[1],
                                   result[1][pitch_parameter],
                                   result[1][depth_parameter],
                                   result[1][slope_parameter],
                                   result[1][power_parameter],
                                   2,
                                   11);

    if (fit0.ok && fit1.ok) {
      hDeltaTimingEventT0->Fill(fit0.t0 - fit1.t0);
      hX0vsX1->Fill(fit0.x, fit1.x);
      hY0vsY1->Fill(fit0.y, fit1.y);
      double dx = fit0.x - fit1.x;
      double dy = fit0.y - fit1.y;
      hPositionDistance->Fill(std::sqrt(dx * dx + dy * dy));
    }

    ++i0;
    ++i1;
  }

  hDeltaTimingEventT0->Write();
  hX0vsX1->Write();
  hY0vsY1->Write();
  hPositionDistance->Write();

  fout.Close();

  std::ofstream conf(outconf);
  conf << "# TIMING propagation model fit\n";
  conf << "# Channel offsets were fixed to the input calibration.\n";
  conf << "# timing pitch_cm depth_cm slope_clock_per_cm power\n";
  for (int timing = 0; timing < 2; ++timing)
    conf << "TIMING" << timing << " "
         << std::setprecision(12)
         << result[timing][pitch_parameter] << " "
         << result[timing][depth_parameter] << " "
         << result[timing][slope_parameter] << " "
         << result[timing][power_parameter] << "\n";
}

} // namespace

int main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string calibration;
  std::string output;
  std::string calibration_output;
  int max_events = 0;
  int max_calls = 5000;
  int print_level = -1;
  int min_channels = 8;
  int grid_steps = 9;
  bool require_trigger = false;
  double pitch = 0.37;
  double depth = 1.;
  double slope = 0.05;
  double power = 0.5;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input", po::value<std::string>(&input)->required(), "triggered-frame input ROOT file")
    ("calibration", po::value<std::string>(&calibration)->required(), "initial TIMING [CHANNEL] calibration file")
    ("output", po::value<std::string>(&output)->default_value("timing_propagation_calib.root"), "output ROOT file")
    ("calibration-output", po::value<std::string>(&calibration_output)->default_value("timing_propagation_calib.conf"), "output propagation model text file")
    ("max-events", po::value<int>(&max_events)->default_value(max_events), "maximum frames per TIMING scintillator, 0 means all")
    ("max-calls", po::value<int>(&max_calls)->default_value(max_calls), "maximum MIGRAD calls")
    ("print-level", po::value<int>(&print_level)->default_value(print_level), "TMinuit print level")
    ("min-channels", po::value<int>(&min_channels)->default_value(min_channels), "minimum channels required in each event fit")
    ("grid-steps", po::value<int>(&grid_steps)->default_value(grid_steps), "grid points per axis for nested event-position scans")
    ("require-trigger", po::bool_switch(&require_trigger), "use only frames containing at least one trigger hit")
    ("pitch", po::value<double>(&pitch)->default_value(pitch), "initial channel pitch in cm")
    ("depth", po::value<double>(&depth)->default_value(depth), "initial effective emission depth in cm")
    ("thickness", po::value<double>(&depth), "deprecated alias for --depth")
    ("slope", po::value<double>(&slope)->default_value(slope), "initial light propagation slope in clock/cm^power")
    ("power", po::value<double>(&power)->default_value(power), "initial light propagation distance power");

  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) {
      std::cout << options << std::endl;
      return 0;
    }
    po::notify(vm);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << "\n\n" << options << std::endl;
    return 1;
  }

  double initial_offset[2][nchannels] = {};
  if (!load_offsets(calibration, initial_offset))
    return 1;

  std::vector<event_t> events[2];
  for (int timing = 0; timing < 2; ++timing) {
    events[timing] = collect_events(input, timing, min_channels, max_events, require_trigger);
    if (events[timing].empty()) {
      std::cerr << "ERROR: no TIMING" << timing << " events found" << std::endl;
      return 1;
    }
  }

  double result[2][nparameters] = {};
  for (int timing = 0; timing < 2; ++timing) {
    if (!fit_timing(timing, events[timing], initial_offset[timing],
                    pitch, depth, slope, power, max_calls, print_level,
                    min_channels, grid_steps, result[timing]))
      return 1;
  }

  write_outputs(output, calibration_output, initial_offset,
                pitch, depth, slope, power, result, events);

  std::cout << "ROOT output:        " << output << std::endl;
  std::cout << "model output:       " << calibration_output << std::endl;
  return 0;
}
