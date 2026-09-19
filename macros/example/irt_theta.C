#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TPaveStats.h>
#include <TTree.h>
#include <TStyle.h>

#include <iostream>
#include <limits>

void irt_theta(const char *filename,
               double min_time,
               double max_time,
               Long64_t spill = -1,
               Long64_t frame = -1)
{
  TFile input(filename, "READ");
  if (input.IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << filename << '\n';
    return;
  }
  auto *irt = dynamic_cast<TTree *>(input.Get("irt"));
  auto *frames = dynamic_cast<TTree *>(input.Get("frames"));
  if (!irt || !frames) {
    std::cerr << "ERROR: input must contain irt and frames trees\n";
    return;
  }
  if (irt->GetEntries() != frames->GetEntries()) {
    std::cerr << "ERROR: irt and frames entry counts differ\n";
    return;
  }

  UShort_t nhits = 0;
  Float_t theta[65535] = {};
  Float_t time[65535] = {};
  UInt_t spill_id = 0;
  if (irt->SetBranchAddress("nhits", &nhits) < 0 ||
      irt->SetBranchAddress("theta", theta) < 0 ||
      irt->SetBranchAddress("time", time) < 0 ||
      frames->SetBranchAddress("spill", &spill_id) < 0) {
    std::cerr << "ERROR: input is missing a required branch\n";
    return;
  }

  auto histogram = new TH1D("hCherenkovThetaSelected",
                            ";#theta (mrad);hits",
                            100, 0., 100.);
  auto histogram_all = new TH1D("hCherenkovThetaAll",
                                ";#theta (mrad);hits",
                                100, 0., 100.);
  histogram->SetDirectory(nullptr);
  histogram_all->SetDirectory(nullptr);
  histogram_all->SetLineColor(kBlack);
  histogram_all->SetLineWidth(1);
  histogram->SetLineColor(kAzure - 3);
  histogram->SetLineWidth(2);
  histogram->SetFillColor(kAzure - 3);
  histogram->SetFillStyle(3002);
  histogram->SetStats(true);
  histogram_all->SetStats(false);
  histogram->GetXaxis()->SetTitleOffset(1.5);
  histogram->GetYaxis()->SetTitleOffset(1.5);
  histogram_all->GetXaxis()->SetTitleOffset(1.5);
  histogram_all->GetYaxis()->SetTitleOffset(1.5);
  Long64_t selected_frames = 0;
  Long64_t selected_hits = 0;
  Long64_t frame_in_spill = 0;
  UInt_t previous_spill = std::numeric_limits<UInt_t>::max();
  for (Long64_t entry = 0; entry < irt->GetEntries(); ++entry) {
    if (frames->GetEntry(entry) <= 0) continue;
    if (spill_id != previous_spill) {
      frame_in_spill = 0;
      previous_spill = spill_id;
    }
    const Long64_t this_frame = frame_in_spill++;
    if (spill >= 0 && spill_id != spill) continue;
    if (frame >= 0 && this_frame != frame) continue;
    if (irt->GetEntry(entry) <= 0) continue;
    ++selected_frames;
    for (unsigned int i = 0; i < nhits; ++i) {
      histogram_all->Fill(1000. * theta[i]);
      if (time[i] < min_time || time[i] > max_time) continue;
      histogram->Fill(1000. * theta[i]);
      ++selected_hits;
    }
  }
  histogram->SetTitle(Form("spill %lld frame %lld;#theta (mrad);hits",
                           spill, frame));
  histogram_all->SetTitle(Form("spill %lld frame %lld;#theta (mrad);hits",
                               spill, frame));

  auto canvas = new TCanvas("cIRTTheta", "", 800, 800);
  canvas->SetMargin(0.15, 0.15, 0.15, 0.15);
  gStyle->SetOptStat(1110);
  canvas->cd();
  histogram->SetMinimum(0.);
  histogram->SetMaximum(1.1 * histogram_all->GetMaximum());
  histogram->Draw("hist");
  histogram_all->Draw("hist same");
  canvas->Modified();
  canvas->Update();
  if (auto *stats = dynamic_cast<TPaveStats *>(histogram->FindObject("stats"))) {
    stats->SetTextColor(kAzure - 3);
    stats->SetLineColor(kAzure - 3);
    stats->SetLineWidth(2);
  }
  canvas->Modified();
  canvas->Update();
  std::cout << "selected frames: " << selected_frames
            << "\nselected hits: " << selected_hits << '\n';
}
