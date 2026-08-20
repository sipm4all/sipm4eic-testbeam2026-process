#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

std::string
default_calibration_filename(const std::string &root_filename)
{
  auto position = root_filename.rfind(".root");
  if (position == std::string::npos)
    return root_filename + ".conf";
  return root_filename.substr(0, position) + ".conf";
}

void
channel_address(int global_channel, int &device, int &fifo, int &column, int &pixel)
{
  device = 192 + global_channel / 256;
  int local_channel = global_channel % 256;
  int chip = local_channel / 32;
  int within_chip = local_channel % 32;
  column = within_chip / 4;
  pixel = within_chip % 4;
  fifo = 4 * chip + column / 2;
}

}

void
fit_calib(std::string infilename,
          std::string outfilename,
          std::string conffilename = "")
{
  auto fin = TFile::Open(infilename.c_str());
  auto hin = (TH2 *)fin->Get("hDeltaT");

  auto c = new TCanvas("c", "c", 800, 800);
  auto hCalib = new TH1F("hCalib", "", 2048, 0., 2048);
  auto hMethod = new TH1I("hMethod", "calibration method;channel;method", 2048, 0., 2048);

  auto fgaus = new TF1("fgaus", "[0] * TMath::Gaus(x, [1], [2])", -2., 2.);
  fgaus->SetParLimits(1, -2., 2.);
  fgaus->SetParLimits(2, 0.1, 0.2);

  double integrale;
  std::array<int, 2> bins;
  bins[0] = hin->GetYaxis()->FindBin(-2. + 0.001);
  bins[1] = hin->GetYaxis()->FindBin(2. - 0.001);
  for (int ich = 0; ich < 2048; ++ich) {
    auto hpy = hin->ProjectionY("hpy", ich + 1, ich + 1);
    hpy->GetXaxis()->SetRange(bins[0], bins[1]);
    auto integral = hpy->IntegralAndError(bins[0], bins[1], integrale);
    if (integral <= 0. || integrale <= 0.) {
      delete hpy;
      continue;
    }

    auto mean = hpy->GetMean();
    auto mean_error = hpy->GetMeanError();
    auto integraler = integrale / integral;
    std::cout << ich << ": integrale / integral = " << integraler << std::endl;

    // These thresholds correspond approximately to effective statistics of
    // 11, 44, 156, and 625 entries through 1 / integraler^2.
    if (integraler > 0.30) {
      hMethod->SetBinContent(ich + 1, 0);
      delete hpy;
      continue;
    }

    if (integraler > 0.15) {
      hCalib->SetBinContent(ich + 1, mean);
      hCalib->SetBinError(ich + 1, mean_error);
      hMethod->SetBinContent(ich + 1, 1);
      delete hpy;
      continue;
    }

    int method = 4;
    if (integraler > 0.08) {
      hpy->Rebin(4);
      method = 2;
    } else if (integraler > 0.04) {
      hpy->Rebin(2);
      method = 3;
    }

    auto bmax = hpy->GetMaximumBin();
    auto cmax = hpy->GetBinCenter(bmax);
    auto vmax = hpy->GetBinContent(bmax);
    fgaus->SetParameter(0, vmax);
    fgaus->SetParameter(1, cmax);
    fgaus->SetParameter(2, 0.15);
    auto fit_status = hpy->Fit(fgaus, "0q", "", cmax - 0.15, cmax + 0.15);
    auto fit_position = fgaus->GetParameter(1);
    auto fit_error = fgaus->GetParError(1);
    if (fit_status != 0 || !TMath::Finite(fit_position) || !TMath::Finite(fit_error)) {
      hCalib->SetBinContent(ich + 1, mean);
      hCalib->SetBinError(ich + 1, mean_error);
      hMethod->SetBinContent(ich + 1, 1);
    } else {
      hCalib->SetBinContent(ich + 1, fit_position);
      hCalib->SetBinError(ich + 1, fit_error);
      hMethod->SetBinContent(ich + 1, method);
    }
    delete hpy;
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  hCalib->Write();
  hMethod->Write();

  if (conffilename.empty())
    conffilename = default_calibration_filename(outfilename);

  std::ofstream calibration(conffilename);
  if (!calibration) {
    std::cerr << "ERROR: could not create calibration file: " << conffilename << std::endl;
  } else {
    calibration << "# Channel offsets derived by fit_calib.C\n"
                << "# Input: " << infilename << "\n"
                << "# calibrated_time = raw_time - offset\n"
                << "# Method 0 channels are omitted because they were skipped.\n"
                << "[CHANNEL]\n"
                << "# device fifo column pixel offset\n"
                << std::fixed << std::setprecision(9);

    for (int ich = 0; ich < 2048; ++ich) {
      if (hMethod->GetBinContent(ich + 1) <= 0.)
        continue;

      int device = 0;
      int fifo = 0;
      int column = 0;
      int pixel = 0;
      channel_address(ich, device, fifo, column, pixel);
      calibration << device << ' ' << fifo << ' ' << column << ' ' << pixel << ' '
                  << hCalib->GetBinContent(ich + 1) << '\n';
    }
  }

  auto hAligned = (TH2 *)hin->Clone("hDeltaT_aligned");
  hAligned->SetTitle("aligned delta-t;channel;delta-t - fitted channel offset");
  hAligned->Reset("ICES");

  for (int ich = 0; ich < 2048; ++ich) {
    bool calibrated = hMethod->GetBinContent(ich + 1) > 0.;
    double offset = calibrated ? hCalib->GetBinContent(ich + 1) : 0.;
    for (int iy = 1; iy <= hin->GetNbinsY(); ++iy) {
      double content = hin->GetBinContent(ich + 1, iy);
      if (content == 0.)
        continue;

      double aligned_delta = hin->GetYaxis()->GetBinCenter(iy) - offset;
      int aligned_bin = hAligned->GetYaxis()->FindBin(aligned_delta);
      if (aligned_bin < 1 || aligned_bin > hAligned->GetNbinsY())
        continue;

      hAligned->AddBinContent(ich + 1, aligned_bin, content);
      double error = hin->GetBinError(ich + 1, iy);
      double previous_error = hAligned->GetBinError(ich + 1, aligned_bin);
      hAligned->SetBinError(ich + 1, aligned_bin,
                            TMath::Sqrt(previous_error * previous_error + error * error));
    }
  }

  hAligned->Write();

  // Search the aligned data for a FIFO-wide secondary population displaced
  // by one native time unit. Eight global channels correspond to one FIFO.
  // The input histograms are normalized, so total_weight is used as a
  // statistics proxy; 0.001 corresponds to roughly 245 events in the usual
  // 244877-frame calibration sample.
  const double fifo_min_weight = 0.001;
  const double fifo_min_fraction = 0.05;
  const double fifo_window = 0.20;
  std::vector<int> flagged_fifos;
  for (int fifo = 0; fifo < 256; ++fifo) {
    int first_channel = 8 * fifo;
    if (first_channel >= hAligned->GetNbinsX())
      break;

    double total_weight = 0.;
    double zero_weight = 0.;
    double plus_one_weight = 0.;
    double minus_one_weight = 0.;
    for (int iy = 1; iy <= hAligned->GetNbinsY(); ++iy) {
      double delta = hAligned->GetYaxis()->GetBinCenter(iy);
      double weight = hAligned->Integral(first_channel + 1,
                                         std::min(first_channel + 8,
                                                  hAligned->GetNbinsX()),
                                         iy, iy);
      total_weight += weight;
      if (std::abs(delta) < fifo_window)
        zero_weight += weight;
      if (std::abs(delta - 1.) < fifo_window)
        plus_one_weight += weight;
      if (std::abs(delta + 1.) < fifo_window)
        minus_one_weight += weight;
    }

    if (total_weight < fifo_min_weight)
      continue;

    double plus_fraction = plus_one_weight / total_weight;
    double minus_fraction = minus_one_weight / total_weight;
    if (plus_fraction >= fifo_min_fraction || minus_fraction >= fifo_min_fraction) {
      flagged_fifos.push_back(fifo);
      std::cout << "FIFO " << fifo
                << ": total=" << total_weight
                << " near0=" << zero_weight
                << " plus1=" << plus_one_weight
                << " minus1=" << minus_one_weight
                << " fractions=(" << plus_fraction
                << "," << minus_fraction << ")" << std::endl;

      auto hFifo = (TH1 *)hAligned->ProjectionY(Form("hDeltaT_fifo%03d", fifo),
                                                first_channel + 1,
                                                std::min(first_channel + 8,
                                                         hAligned->GetNbinsX()));
      hFifo->SetTitle(Form("aligned delta-t, FIFO %d;delta-t;weighted entries", fifo));
      hFifo->Write();
      delete hFifo;
    }
  }
  std::cout << "FIFO one-cycle diagnostics: " << flagged_fifos.size()
            << " selected" << std::endl;

  fout->Close();
  fin->Close();
}
