#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

void irt_theta_track(const char *filename,
                     double min_time, double max_time,
                     double min_theta, double max_theta,
                     double min_x0, double max_x0,
                     double min_y0, double max_y0,
                     double min_r, double max_r,
                     int min_hits, int max_hits,
                     Long64_t max_frames = -1)
{
  TFile input(filename, "READ");
  if (input.IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << filename << '\n';
    return;
  }
  auto *irt = dynamic_cast<TTree *>(input.Get("irt"));
  auto *ring = dynamic_cast<TTree *>(input.Get("ring"));
  if (!irt || !ring || irt->GetEntries() != ring->GetEntries()) {
    std::cerr << "ERROR: input must contain aligned irt and ring trees\n";
    return;
  }

  UShort_t nhits = 0;
  Float_t theta[65535] = {};
  Float_t time[65535] = {};
  UChar_t nrings = 0;
  Float_t x0[256] = {}, y0[256] = {}, radius[256] = {};
  if (irt->SetBranchAddress("nhits", &nhits) < 0 ||
      irt->SetBranchAddress("theta", theta) < 0 ||
      irt->SetBranchAddress("time", time) < 0 ||
      ring->SetBranchAddress("nrings", &nrings) < 0 ||
      ring->SetBranchAddress("x0", x0) < 0 ||
      ring->SetBranchAddress("y0", y0) < 0 ||
      ring->SetBranchAddress("r", radius) < 0) {
    std::cerr << "ERROR: input is missing a required branch\n";
    return;
  }

  const std::string output_name = std::string(filename) + ".irt_theta_track.root";
  TFile output(output_name.c_str(), "RECREATE");
  TH2D hMeanVsHits("hMeanThetaVsHits",
                   ";selected hits;Gaussian mean #theta (mrad)",
                   std::max(100, max_hits > 0 ? max_hits : 100),
                   0., std::max(100, max_hits > 0 ? max_hits : 100),
                   200, 35., 45.);
  TH1D hMeans("hFittedThetaMean", ";Gaussian mean #theta (mrad);frames", 200, 35., 45.);
  TH1D hWidths("hFittedThetaSigma", ";Gaussian #sigma #theta (mrad);frames", 100, 0., 20.);
  TH1D frame_histogram("hThetaFrame", "", 100,
                       min_theta * 1000., max_theta * 1000.);
  frame_histogram.SetDirectory(nullptr);
  TF1 gaussian("gaussian", "gaus", min_theta * 1000., max_theta * 1000.);
  gaussian.SetParLimits(1, min_theta * 1000., max_theta * 1000.);
  gaussian.SetParLimits(2, 0.7, 2.0);
  Long64_t selected_events = 0;
  Long64_t fitted_events = 0;
  Long64_t mean_limit_overshoots = 0;
  Long64_t sigma_limit_overshoots = 0;

  const Long64_t entries = max_frames < 0 ? irt->GetEntries() :
      std::min<Long64_t>(irt->GetEntries(), max_frames);
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (irt->GetEntry(entry) <= 0 || ring->GetEntry(entry) <= 0) continue;
    if (nrings != 1) continue;
    if (x0[0] < min_x0 || x0[0] > max_x0 ||
        y0[0] < min_y0 || y0[0] > max_y0 ||
        radius[0] < min_r || radius[0] > max_r)
      continue;

    std::vector<double> values;
    values.reserve(nhits);
    for (unsigned int i = 0; i < nhits; ++i) {
      if (time[i] < min_time || time[i] > max_time ||
          theta[i] < min_theta || theta[i] > max_theta)
        continue;
      values.push_back(1000. * theta[i]);
    }
    if (static_cast<int>(values.size()) < min_hits ||
        (max_hits >= 0 && static_cast<int>(values.size()) > max_hits))
      continue;
    ++selected_events;

    frame_histogram.Reset();
    for (const double value : values) frame_histogram.Fill(value);
    gaussian.SetParameters(frame_histogram.GetMaximum(), 38.,
                           1.3);
    gaussian.SetParLimits(1, min_theta * 1000., max_theta * 1000.);
    gaussian.SetParLimits(2, 0.7, 2.0);
    const auto fit = frame_histogram.Fit(&gaussian, "QLS");
    if (!fit.Get() || !fit->IsValid()) continue;
    ++fitted_events;
    const double raw_mean = fit->Parameter(1);
    const double raw_sigma = std::abs(fit->Parameter(2));
    if (raw_mean < min_theta * 1000. || raw_mean > max_theta * 1000.)
      ++mean_limit_overshoots;
    if (raw_sigma < 0.7 || raw_sigma > 2.0)
      ++sigma_limit_overshoots;
    const double mean = std::clamp(raw_mean,
                                   min_theta * 1000., max_theta * 1000.);
    const double sigma = std::clamp(raw_sigma, 0.7, 2.0);
    hMeanVsHits.Fill(static_cast<double>(values.size()), mean);
    hMeans.Fill(mean);
    hWidths.Fill(sigma);
  }

  hMeanVsHits.Write();
  hMeans.Write();
  hWidths.Write();
  TCanvas canvas("cThetaTrack", "", 800, 800);
  canvas.SetMargin(0.15, 0.15, 0.15, 0.15);
  hMeanVsHits.Draw("colz");
  canvas.Write();
  output.Close();
  std::cout << "selected events: " << selected_events
            << "\nfitted events: " << fitted_events
            << "\nraw mean limit overshoots: " << mean_limit_overshoots
            << "\nraw sigma limit overshoots: " << sigma_limit_overshoots
            << "\noutput: " << output_name << '\n';
}
