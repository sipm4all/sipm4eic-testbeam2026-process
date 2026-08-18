#include "../lib/trigger_reader.h"

#include <TCanvas.h>
#include <TH2D.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

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

bool load_offsets(const char *filename, double offset[2][nchannels])
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
        std::cerr << "ERROR: missing offset for TIMING" << timing
                  << " DO " << channel << std::endl;
        ok = false;
      }
    }
  }

  return ok;
}

void setup_hist(TH2D *hist)
{
  hist->SetStats(false);
  hist->SetMinimum(-2.);
  hist->SetMaximum(8.);
  for (int ix = 1; ix <= 4; ++ix)
    hist->GetXaxis()->SetBinLabel(ix, TString::Format("x%d", ix - 1));
  for (int iy = 1; iy <= 8; ++iy)
    hist->GetYaxis()->SetBinLabel(iy, TString::Format("y%d", iy - 1));
}

void draw_labels(const bool have[nchannels], const double value[nchannels])
{
  TLatex text;
  text.SetTextAlign(22);
  text.SetTextSize(0.035);
  text.SetTextColor(kBlack);

  for (int channel = 0; channel < nchannels; ++channel) {
    if (!have[channel])
      continue;

    int x = channel % 4;
    int y = channel / 4;
    text.DrawLatex(x, y, TString::Format("%.2f", value[channel]));
  }
}

} // namespace

void timing_event_display(const char *filename,
                          const char *calibfilename = "timing_offsets_from_cherenkov860.conf",
                          int first_event = 0,
                          int max_events = 0,
                          bool relative_time = true,
                          bool require_trigger = false)
{
  double offset[2][nchannels] = {};
  if (!load_offsets(calibfilename, offset))
    return;

  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  gStyle->SetOptStat(0);
  auto canvas = new TCanvas("cTimingEventDisplay", "TIMING event display", 1100, 850);
  canvas->Divide(2, 2);

  auto hTiming0 = new TH2D("hTiming0Event",
                           "TIMING0;DO x;DO y;time [clock]",
                           4, -0.5, 3.5, 8, -0.5, 7.5);
  auto hTiming1 = new TH2D("hTiming1Event",
                           "TIMING1;DO x;DO y;time [clock]",
                           4, -0.5, 3.5, 8, -0.5, 7.5);
  auto hTime0 = new TH1D("hTiming0EventTime",
                         "TIMING0 time distribution;time [clock];channels",
                         100, -2., 8.);
  auto hTime1 = new TH1D("hTiming1EventTime",
                         "TIMING1 time distribution;time [clock];channels",
                         100, -2., 8.);
  setup_hist(hTiming0);
  setup_hist(hTiming1);
  hTime0->SetStats(false);
  hTime1->SetStats(false);

  Long64_t iframe_global = 0;
  int displayed = 0;

  while (reader.next_spill()) {
    while (reader.next_frame()) {
      if (require_trigger && reader.trigger_hits().empty()) {
        ++iframe_global;
        continue;
      }

      bool have[2][nchannels] = {};
      double time[2][nchannels] = {};
      double reference = std::numeric_limits<double>::max();

      for (const auto &hit : reader.timing_hits()) {
        int timing = timing_detector(hit);
        if (timing < 0)
          continue;

        int eoch = electronics_channel(hit);
        if (eoch < 0)
          continue;

        int doch = eo2do[eoch];
        double calibrated = hit.time - offset[timing][doch];
        if (!have[timing][doch] || calibrated < time[timing][doch]) {
          have[timing][doch] = true;
          time[timing][doch] = calibrated;
        }
      }

      bool any = false;
      for (int timing = 0; timing < 2; ++timing) {
        for (int channel = 0; channel < nchannels; ++channel) {
          if (!have[timing][channel])
            continue;
          any = true;
          if (time[timing][channel] < reference)
            reference = time[timing][channel];
        }
      }

      if (!any) {
        ++iframe_global;
        continue;
      }

      if ((int)iframe_global < first_event) {
        ++iframe_global;
        continue;
      }

      hTiming0->Reset("ICES");
      hTiming1->Reset("ICES");
      hTime0->Reset("ICES");
      hTime1->Reset("ICES");

      double display_time[2][nchannels] = {};
      for (int channel = 0; channel < nchannels; ++channel) {
        int x = channel % 4;
        int y = channel / 4;

        if (have[0][channel]) {
          display_time[0][channel] =
            relative_time ? time[0][channel] - reference : time[0][channel];
          hTiming0->SetBinContent(x + 1, y + 1, display_time[0][channel]);
          hTime0->Fill(display_time[0][channel]);
        }
        if (have[1][channel]) {
          display_time[1][channel] =
            relative_time ? time[1][channel] - reference : time[1][channel];
          hTiming1->SetBinContent(x + 1, y + 1, display_time[1][channel]);
          hTime1->Fill(display_time[1][channel]);
        }
      }

      auto title = TString::Format("spill %d frame %d global %lld  reference %.6f",
                                   reader.spill_id(), reader.frame_index(),
                                   iframe_global, reference);
      hTiming0->SetTitle(TString::Format("TIMING0 %s;DO x;DO y;%s",
                                         title.Data(),
                                         relative_time ? "t - t_{first} [clock]" : "calibrated time [clock]"));
      hTiming1->SetTitle(TString::Format("TIMING1 %s;DO x;DO y;%s",
                                         title.Data(),
                                         relative_time ? "t - t_{first} [clock]" : "calibrated time [clock]"));
      hTime0->SetTitle(TString::Format("TIMING0 time distribution %s;%s;channels",
                                       title.Data(),
                                       relative_time ? "t - t_{first} [clock]" : "calibrated time [clock]"));
      hTime1->SetTitle(TString::Format("TIMING1 time distribution %s;%s;channels",
                                       title.Data(),
                                       relative_time ? "t - t_{first} [clock]" : "calibrated time [clock]"));

      canvas->cd(1);
      hTiming0->Draw("COLZ");
      draw_labels(have[0], display_time[0]);
      canvas->cd(2);
      hTiming1->Draw("COLZ");
      draw_labels(have[1], display_time[1]);
      canvas->cd(3);
      hTime0->Draw("HIST");
      canvas->cd(4);
      hTime1->Draw("HIST");
      canvas->Modified();
      canvas->Update();

      std::cout << "event " << displayed
                << " spill=" << reader.spill_id()
                << " frame=" << reader.frame_index()
                << " global=" << iframe_global
                << " reference=" << std::setprecision(12) << reference
                << " trigger_hits=" << reader.trigger_hits().size()
                << std::endl;

      ++displayed;
      ++iframe_global;

      if (max_events > 0 && displayed >= max_events)
        return;

      std::cout << "press Enter for next event, q + Enter to stop: " << std::flush;
      std::string line;
      std::getline(std::cin, line);
      if (!line.empty() && (line[0] == 'q' || line[0] == 'Q'))
        return;
    }
  }
}
