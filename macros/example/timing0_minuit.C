#include "../lib/trigger_reader.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMinuit.h>
#include <TString.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace timing0_minuit_detail {

constexpr int nchannels = 32;
constexpr int reference_eoch = 0;

const int eo2do[nchannels] = {22, 20, 18, 16, 24, 26, 28, 30,
                              25, 27, 29, 31, 23, 21, 19, 17,
                              9,  11, 13, 15, 7,  5,  3,  1,
                              6,  4,  2,  0,  8,  10, 12, 14};

struct pair_t {
  int a = 0;
  int b = 0;
  double dt = 0.;
};

std::vector<pair_t> pairs;

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

int do_x(int doch)
{
  return doch % 4;
}

int do_y(int doch)
{
  return doch / 4;
}

int do_distance(int eoch_a, int eoch_b)
{
  int doch_a = eo2do[eoch_a];
  int doch_b = eo2do[eoch_b];
  return std::abs(do_x(doch_a) - do_x(doch_b)) +
         std::abs(do_y(doch_a) - do_y(doch_b));
}

int par_to_eoch(int ipar)
{
  return ipar < reference_eoch ? ipar : ipar + 1;
}

double corrected_dt(const pair_t &pair, const double *offset)
{
  return pair.dt - offset[pair.a] + offset[pair.b];
}

double objective(const double *offset)
{
  if (pairs.empty())
    return 1.e30;

  double chi2 = 0.;
  for (const auto &pair : pairs) {
    double dt = corrected_dt(pair, offset);
    chi2 += dt * dt;
  }
  return chi2 / pairs.size();
}

void fcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  double offset[nchannels] = {};
  offset[reference_eoch] = 0.;
  for (int ipar = 0; ipar < nchannels - 1; ++ipar)
    offset[par_to_eoch(ipar)] = par[ipar];
  f = objective(offset);
}

void write_channel_calibration(const char *filename, const double *offset)
{
  if (!filename || std::string(filename).empty())
    return;

  std::ofstream out(filename);
  if (!out) {
    std::cerr << "ERROR: could not create calibration file: " << filename << std::endl;
    return;
  }

  out << "[CHANNEL]\n";
  out << "# device fifo column pixel offset\n";
  out << "# TIMING0 offsets from timing0_minuit.C\n";
  out << "# calibrated_time = raw_time - offset\n";

  out << std::setprecision(12);
  for (int eoch = 0; eoch < nchannels; ++eoch) {
    int column = eoch / 4;
    int pixel = eoch % 4;
    for (int fifo = 0; fifo < 4; ++fifo)
      out << "200 " << fifo << " " << column << " " << pixel
          << " " << offset[eoch] << "\n";
  }
}

} // namespace timing0_minuit_detail

void timing0_minuit(const char *filename,
                    const char *outfilename = "timing0_minuit.root",
                    const char *calibfilename = "timing0_offsets.conf",
                    int max_calls = 5000,
                    double delta_range = 32.)
{
  using namespace timing0_minuit_detail;

  trigger_reader_t reader;
  if (!reader.open(filename))
    return;

  pairs.clear();

  Long64_t nframes = 0;
  Long64_t nframes_with_timing0 = 0;
  Long64_t nframes_used = 0;

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

        if (!have[eoch] || hit.time < time[eoch]) {
          have[eoch] = true;
          time[eoch] = hit.time;
        }
      }

      bool any = false;
      for (int eoch = 0; eoch < nchannels; ++eoch)
        any = any || have[eoch];
      if (any)
        ++nframes_with_timing0;

      bool used = false;
      for (int a = 0; a < nchannels; ++a) {
        if (!have[a])
          continue;
        for (int b = a + 1; b < nchannels; ++b) {
          if (!have[b])
            continue;
          if (do_distance(a, b) > 1)
            continue;

          pair_t pair;
          pair.a = a;
          pair.b = b;
          pair.dt = time[a] - time[b];
          pairs.push_back(pair);
          used = true;
        }
      }

      if (used)
        ++nframes_used;
    }
  }

  if (pairs.empty()) {
    std::cerr << "ERROR: no TIMING0 neighbour pairs found" << std::endl;
    return;
  }

  double offset[nchannels] = {};
  double objective_before = objective(offset);

  TMinuit minuit(nchannels - 1);
  minuit.SetFCN(fcn);

  Double_t arglist[2] = {-1., 0.};
  Int_t ierflg = 0;
  minuit.mnexcm("SET PRINT", arglist, 1, ierflg);

  for (int ipar = 0; ipar < nchannels - 1; ++ipar) {
    int eoch = par_to_eoch(ipar);
    auto name = TString::Format("offset_eoch_%02d", eoch);
    minuit.DefineParameter(ipar, name.Data(),
                           0., 0.01, 0., 0.);
  }

  arglist[0] = max_calls;
  arglist[1] = 0.01;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

  offset[reference_eoch] = 0.;
  for (int ipar = 0; ipar < nchannels - 1; ++ipar) {
    int eoch = par_to_eoch(ipar);
    Double_t value = 0.;
    Double_t error = 0.;
    minuit.GetParameter(ipar, value, error);
    offset[eoch] = value;
  }

  double objective_after = objective(offset);

  auto fout = TFile::Open(outfilename, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return;
  }

  auto hDeltaBefore = new TH1D("hDeltaBefore", "", 2048, -delta_range, delta_range);
  auto hDeltaAfter = new TH1D("hDeltaAfter", "", 2048, -delta_range, delta_range);
  auto hOffsetEo = new TH1D("hOffsetEo", "", nchannels, 0., nchannels);
  auto hOffsetDo = new TH1D("hOffsetDo", "", nchannels, 0., nchannels);
  auto hPairCountDo = new TH2D("hPairCountDo", "", nchannels, 0., nchannels,
                               nchannels, 0., nchannels);

  TH1D *hPairBefore[nchannels][nchannels] = {};
  TH1D *hPairAfter[nchannels][nchannels] = {};
  for (int a = 0; a < nchannels; ++a) {
    for (int b = a + 1; b < nchannels; ++b) {
      if (do_distance(a, b) > 1)
        continue;
      int doa = eo2do[a];
      int dob = eo2do[b];
      int lo = std::min(doa, dob);
      int hi = std::max(doa, dob);
      auto name_before = TString::Format("deltat_%02d_%02d_timing0_before", lo, hi);
      auto name_after = TString::Format("deltat_%02d_%02d_timing0_after", lo, hi);
      hPairBefore[a][b] = new TH1D(name_before.Data(), "", 2048, -delta_range, delta_range);
      hPairAfter[a][b] = new TH1D(name_after.Data(), "", 2048, -delta_range, delta_range);
    }
  }

  for (int eoch = 0; eoch < nchannels; ++eoch) {
    int doch = eo2do[eoch];
    hOffsetEo->SetBinContent(eoch + 1, offset[eoch]);
    hOffsetDo->SetBinContent(doch + 1, offset[eoch]);
  }

  for (const auto &pair : pairs) {
    double before = pair.dt;
    double after = corrected_dt(pair, offset);
    int doa = eo2do[pair.a];
    int dob = eo2do[pair.b];

    hDeltaBefore->Fill(before);
    hDeltaAfter->Fill(after);
    hPairCountDo->Fill(doa, dob);
    hPairCountDo->Fill(dob, doa);

    hPairBefore[pair.a][pair.b]->Fill(before);
    hPairAfter[pair.a][pair.b]->Fill(after);
  }

  hDeltaBefore->Write();
  hDeltaAfter->Write();
  hOffsetEo->Write();
  hOffsetDo->Write();
  hPairCountDo->Write();
  for (int a = 0; a < nchannels; ++a) {
    for (int b = a + 1; b < nchannels; ++b) {
      if (!hPairBefore[a][b])
        continue;
      hPairBefore[a][b]->Write();
      hPairAfter[a][b]->Write();
    }
  }
  fout->Close();

  write_channel_calibration(calibfilename, offset);

  std::cout << "frames processed:          " << nframes << std::endl;
  std::cout << "frames with TIMING0:       " << nframes_with_timing0 << std::endl;
  std::cout << "frames used:               " << nframes_used << std::endl;
  std::cout << "neighbour pair entries:    " << pairs.size() << std::endl;
  std::cout << "reference offset fixed:    EOCH " << reference_eoch << " = 0" << std::endl;
  std::cout << "Minuit status:             " << (ierflg == 0 ? "OK" : "ERROR") << std::endl;
  std::cout << "objective before/after:    " << objective_before << " / "
            << objective_after << std::endl;
  std::cout << "ROOT output:               " << outfilename << std::endl;
  std::cout << "calibration output:        " << calibfilename << std::endl;
}
