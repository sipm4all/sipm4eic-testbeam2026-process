#include "../lib/trigger_reader.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TF2.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <TRandom.h>
#include <Math/Functor.h>
#include <Fit/Fitter.h>
#include <map>
#include <array>
#include <cmath>
#include <iostream>
#include <string>

std::map<std::string, std::array<double, 2>>
hitmap_fit(std::string input_filename, std::string output_filename)
{
  trigger_reader_t reader;
  if (!reader.open(input_filename)) return {};

  TH2F hmap("hMap", ";x (mm);y (mm)", 396, -99., 99., 396, -99., 99.);
  TH2F hmap_vis("hMap_vis", ";x (mm);y (mm)", 396, -99., 99., 396, -99., 99.);
  std::map<std::array<float, 2>, double> hitmap;
  std::map<std::array<float, 2>, double> fit_hitmap;
  int selected_frames = 0;
  while (reader.next_spill()) {
    while (reader.next_frame()) {
      ++selected_frames;
      for (const auto &hit : reader.cherenkov_hits()) {
        if (!std::isfinite(hit.x) || !std::isfinite(hit.y) ||
            !std::isfinite(hit.time)) continue;
        const std::array<float, 2> pos = {static_cast<float>(hit.x), static_cast<float>(hit.y)};
        hitmap[pos] += 1.;
        hmap.Fill(hit.x, hit.y);
        hmap_vis.Fill(gRandom->Uniform(hit.x - 1.5, hit.x + 1.5),
                      gRandom->Uniform(hit.y - 1.5, hit.y + 1.5));
        if (hit.x >= -36.)
          fit_hitmap[pos] += 1.;
      }
    }
  }
  if (selected_frames == 0 || fit_hitmap.empty()) return {};
  hmap.Scale(1. / selected_frames);

  TF2 model("fmodel",
    "[0]/(2.*3.1415927*[3]*[4])*TMath::Gaus(TMath::Sqrt("
    "(((x-[1])*TMath::Cos([5])+(y-[2])*TMath::Sin([5]))/[3])^2+"
    "((-(x-[1])*TMath::Sin([5])+(y-[2])*TMath::Cos([5]))/[4])^2),1.,[6],1)"
    "+[7]/(200.*200.)", -99., 99., -99., 99.);
  auto chi2 = [&](const double *p) {
    for (int i = 0; i < 8; ++i) model.SetParameter(i, p[i]);
    double sum = 0.;
    for (const auto &[pos, count] : fit_hitmap) {
      const double prediction = 1. - TMath::PoissonI(0,
        model.Integral(pos[0] - 1.5, pos[0] + 1.5,
                       pos[1] - 1.5, pos[1] + 1.5));
      const double value = count / selected_frames;
      const double error = std::sqrt(count) / selected_frames;
      if (error > 0.) sum += std::pow((value - prediction) / error, 2);
    }
    return sum;
  };
  const double start[8] = {10., 30., -30., 45., 35., 0., .05, 0.};
  ROOT::Fit::Fitter fitter;
  ROOT::Math::Functor fcn(chi2, 8);
  fitter.SetFCN(fcn, start);
  const char *names[] = {"Nsig", "X0", "Y0", "A", "B", "theta", "sigmaRho", "Nbkg"};
  const double limits[][2] = {{0,100},{-100,100},{-100,100},{0,100},{0,100},
                              {-TMath::Pi()/2,TMath::Pi()/2},{.005,.20},{0,100}};
  for (int i = 0; i < 8; ++i) {
    fitter.Config().ParSettings(i).SetName(names[i]);
    fitter.Config().ParSettings(i).SetLimits(limits[i][0], limits[i][1]);
  }
  const bool fit_ok = fitter.FitFCN();
  const auto result = fitter.Result();
  result.Print(std::cout);
  if (!fit_ok) return {};

  TFile output(output_filename.c_str(), "RECREATE");
  auto *results = new TH1F("hResults", "fit results", 8, 0, 8);
  for (int i = 0; i < 8; ++i) {
    results->GetXaxis()->SetBinLabel(i + 1, names[i]);
    results->SetBinContent(i + 1, result.Parameter(i));
    results->SetBinError(i + 1, result.ParError(i));
  }
  TCanvas canvas("cResults", "hitmap ellipse fit", 1600, 800);
  canvas.Divide(2, 1);
  canvas.cd(1)->SetLogz();
  hmap_vis.Draw("col");
  const double angle = result.Parameter(5) * 180. / TMath::Pi();
  for (int sign : {-1, 1}) {
    auto *ellipse = new TEllipse(result.Parameter(1), result.Parameter(2),
      result.Parameter(3) * (1. + sign * result.Parameter(6)),
      result.Parameter(4) * (1. + sign * result.Parameter(6)), 0., 360., angle);
    ellipse->SetFillStyle(0); ellipse->SetLineColor(kBlack); ellipse->Draw("same");
  }
  canvas.cd(2);
  TLatex text; text.SetNDC();
  for (int i = 0; i < 8; ++i)
    text.DrawLatex(.15, .9 - .09 * i,
      Form("%s = %.6g #pm %.3g", names[i], result.Parameter(i), result.ParError(i)));
  results->Write();
  hmap.Write(); hmap_vis.Write(); canvas.Write();
  canvas.SaveAs((input_filename + ".hitmap_fit.png").c_str());
  output.Close();
  return {{"Nsig", {result.Parameter(0), result.ParError(0)}},
          {"X0", {result.Parameter(1), result.ParError(1)}},
          {"Y0", {result.Parameter(2), result.ParError(2)}},
          {"A", {result.Parameter(3), result.ParError(3)}},
          {"B", {result.Parameter(4), result.ParError(4)}},
          {"theta", {result.Parameter(5), result.ParError(5)}},
          {"sigmaRho", {result.Parameter(6), result.ParError(6)}},
          {"Nbkg", {result.Parameter(7), result.ParError(7)}}};
}
