#include "../../macros/lib/trigger_reader.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <boost/program_options.hpp>

#include <algorithm>
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

const int eo2do[nchannels] = {22, 20, 18, 16, 24, 26, 28, 30,
                              25, 27, 29, 31, 23, 21, 19, 17,
                              9,  11, 13, 15, 7,  5,  3,  1,
                              6,  4,  2,  0,  8,  10, 12, 14};

struct event_t {
  std::vector<int> channel;
  std::vector<double> time;
};

struct fit_result_t {
  bool ok = false;
  double x = 0.;
  double y = 0.;
  double t0 = 0.;
  double rms = 0.;
  int n = 0;
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

int electronics_channel(const hit_t &hit)
{
  int channel = hit.pixel + 4 * hit.column;
  if (channel < 0 || channel >= nchannels)
    return -1;
  return channel;
}

double channel_x(int doch, double pitch)
{
  return pitch * (doch % 4);
}

double channel_y(int doch, double pitch)
{
  return pitch * (doch / 4);
}

double path_length(double x, double y, int doch, double pitch, double thickness)
{
  double dx = channel_x(doch, pitch) - x;
  double dy = channel_y(doch, pitch) - y;
  return 2. * std::sqrt(dx * dx + dy * dy + thickness * thickness);
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
        std::cerr << "ERROR: missing TIMING" << timing
                  << " offset for DO " << channel << std::endl;
        ok = false;
      }
    }
  }
  return ok;
}

double objective(const event_t &event,
                 const double offset[nchannels],
                 double pitch,
                 double thickness,
                 double slope,
                 double x,
                 double y,
                 double &t0)
{
  double sum = 0.;
  int n = 0;
  for (size_t i = 0; i < event.channel.size(); ++i) {
    int ch = event.channel[i];
    sum += event.time[i] - offset[ch] - slope * path_length(x, y, ch, pitch, thickness);
    ++n;
  }

  if (n == 0) {
    t0 = 0.;
    return std::numeric_limits<double>::max();
  }

  t0 = sum / n;

  double sum2 = 0.;
  for (size_t i = 0; i < event.channel.size(); ++i) {
    int ch = event.channel[i];
    double residual = event.time[i] - offset[ch] -
                      (t0 + slope * path_length(x, y, ch, pitch, thickness));
    sum2 += residual * residual;
  }

  return sum2 / n;
}

fit_result_t fit_event(const event_t &event,
                       const double offset[nchannels],
                       double pitch,
                       double thickness,
                       double slope,
                       int min_channels)
{
  fit_result_t result;
  result.n = event.channel.size();
  if ((int)event.channel.size() < min_channels)
    return result;

  double xmin = -0.5 * pitch;
  double xmax = 3.5 * pitch;
  double ymin = -0.5 * pitch;
  double ymax = 7.5 * pitch;

  double best_x = 0.;
  double best_y = 0.;
  double best_t0 = 0.;
  double best = std::numeric_limits<double>::max();

  auto scan = [&](double cx, double cy, double half_range, int nstep) {
    for (int ix = 0; ix < nstep; ++ix) {
      double x = cx - half_range + 2. * half_range * ix / std::max(1, nstep - 1);
      x = std::clamp(x, xmin, xmax);
      for (int iy = 0; iy < nstep; ++iy) {
        double y = cy - half_range + 2. * half_range * iy / std::max(1, nstep - 1);
        y = std::clamp(y, ymin, ymax);
        double t0 = 0.;
        double value = objective(event, offset, pitch, thickness, slope, x, y, t0);
        if (value < best) {
          best = value;
          best_x = x;
          best_y = y;
          best_t0 = t0;
        }
      }
    }
  };

  scan(0.5 * (xmin + xmax), 0.5 * (ymin + ymax), std::max(xmax - xmin, ymax - ymin), 17);
  scan(best_x, best_y, 1.5 * pitch, 17);
  scan(best_x, best_y, 0.35 * pitch, 17);
  scan(best_x, best_y, 0.08 * pitch, 17);

  result.ok = true;
  result.x = best_x;
  result.y = best_y;
  result.t0 = best_t0;
  result.rms = std::sqrt(best);
  return result;
}

} // namespace

int main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string calibration;
  std::string output;
  int max_events = 0;
  int min_channels = 8;
  double pitch = 0.37;
  double thickness = 2.;
  double slope = 0.05;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input", po::value<std::string>(&input)->required(), "triggered-frame input ROOT file")
    ("calibration", po::value<std::string>(&calibration)->required(), "fixed TIMING [CHANNEL] calibration file")
    ("output", po::value<std::string>(&output)->default_value("timing_position_fit.root"), "output ROOT file")
    ("max-events", po::value<int>(&max_events)->default_value(max_events), "maximum frames to process, 0 means all")
    ("min-channels", po::value<int>(&min_channels)->default_value(min_channels), "minimum channels required per TIMING fit")
    ("pitch", po::value<double>(&pitch)->default_value(pitch), "channel pitch in cm")
    ("thickness", po::value<double>(&thickness)->default_value(thickness), "scintillator thickness in cm")
    ("slope", po::value<double>(&slope)->default_value(slope), "light propagation slope in clock/cm");

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

  double offset[2][nchannels] = {};
  if (!load_offsets(calibration, offset))
    return 1;

  trigger_reader_t reader;
  if (!reader.open(input))
    return 1;

  auto fout = TFile::Open(output.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << output << std::endl;
    return 1;
  }

  TH2D *hYvsX[2] = {};
  TH1D *hX[2] = {};
  TH1D *hY[2] = {};
  TH1D *hT0[2] = {};
  TH1D *hResidual[2] = {};
  TH2D *hResidualVsN[2] = {};

  for (int timing = 0; timing < 2; ++timing) {
    hYvsX[timing] = new TH2D(TString::Format("hYvsX_timing%d", timing),
                             ";x [cm];y [cm]",
                             120, -0.5 * pitch, 3.5 * pitch,
                             200, -0.5 * pitch, 7.5 * pitch);
    hX[timing] = new TH1D(TString::Format("hX_timing%d", timing),
                          ";x [cm];events", 120, -0.5 * pitch, 3.5 * pitch);
    hY[timing] = new TH1D(TString::Format("hY_timing%d", timing),
                          ";y [cm];events", 200, -0.5 * pitch, 7.5 * pitch);
    hT0[timing] = new TH1D(TString::Format("hT0_timing%d", timing),
                           ";fitted t0 [clock];events", 4096, -32., 32.);
    hResidual[timing] = new TH1D(TString::Format("hResidual_timing%d", timing),
                                 ";fit residual RMS [clock];events", 1000, 0., 2.);
    hResidualVsN[timing] = new TH2D(TString::Format("hResidualVsN_timing%d", timing),
                                    ";number of channels;fit residual RMS [clock]",
                                    33, -0.5, 32.5, 1000, 0., 2.);
  }

  auto hDeltaTiming = new TH1D("hDeltaTiming", ";t0 TIMING0 - t0 TIMING1 [clock];events",
                               2048, -32., 32.);
  auto hX0vsX1 = new TH2D("hX0vsX1", ";TIMING0 x [cm];TIMING1 x [cm]",
                          120, -0.5 * pitch, 3.5 * pitch,
                          120, -0.5 * pitch, 3.5 * pitch);
  auto hY0vsY1 = new TH2D("hY0vsY1", ";TIMING0 y [cm];TIMING1 y [cm]",
                          200, -0.5 * pitch, 7.5 * pitch,
                          200, -0.5 * pitch, 7.5 * pitch);
  auto hPositionDistance = new TH1D("hPositionDistance",
                                    ";distance between TIMING0 and TIMING1 fitted positions [cm];events",
                                    1000, 0., 4.);

  Long64_t frames = 0;
  Long64_t accepted[2] = {};
  Long64_t accepted_pair = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++frames;

      event_t events[2];
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
        for (int ch = 0; ch < nchannels; ++ch) {
          if (!have[timing][ch])
            continue;
          events[timing].channel.push_back(ch);
          events[timing].time.push_back(time[timing][ch]);
        }
      }

      fit_result_t fit[2];
      for (int timing = 0; timing < 2; ++timing) {
        fit[timing] = fit_event(events[timing], offset[timing],
                                pitch, thickness, slope, min_channels);
        if (!fit[timing].ok)
          continue;
        ++accepted[timing];
        hYvsX[timing]->Fill(fit[timing].x, fit[timing].y);
        hX[timing]->Fill(fit[timing].x);
        hY[timing]->Fill(fit[timing].y);
        hT0[timing]->Fill(fit[timing].t0);
        hResidual[timing]->Fill(fit[timing].rms);
        hResidualVsN[timing]->Fill(fit[timing].n, fit[timing].rms);
      }

      if (fit[0].ok && fit[1].ok) {
        ++accepted_pair;
        hDeltaTiming->Fill(fit[0].t0 - fit[1].t0);
        hX0vsX1->Fill(fit[0].x, fit[1].x);
        hY0vsY1->Fill(fit[0].y, fit[1].y);
        double dx = fit[0].x - fit[1].x;
        double dy = fit[0].y - fit[1].y;
        hPositionDistance->Fill(std::sqrt(dx * dx + dy * dy));
      }

      if (max_events > 0 && frames >= max_events)
        break;
    }
    if (max_events > 0 && frames >= max_events)
      break;
  }

  fout->cd();
  for (int timing = 0; timing < 2; ++timing) {
    hYvsX[timing]->Write();
    hX[timing]->Write();
    hY[timing]->Write();
    hT0[timing]->Write();
    hResidual[timing]->Write();
    hResidualVsN[timing]->Write();
  }
  hDeltaTiming->Write();
  hX0vsX1->Write();
  hY0vsY1->Write();
  hPositionDistance->Write();
  fout->Close();

  std::cout << "frames processed:        " << frames << std::endl;
  std::cout << "accepted TIMING0:        " << accepted[0] << std::endl;
  std::cout << "accepted TIMING1:        " << accepted[1] << std::endl;
  std::cout << "accepted pairs:          " << accepted_pair << std::endl;
  std::cout << "pitch [cm]:              " << pitch << std::endl;
  std::cout << "thickness [cm]:          " << thickness << std::endl;
  std::cout << "slope [clock/cm]:        " << slope << std::endl;
  std::cout << "ROOT output:             " << output << std::endl;
  return 0;
}
