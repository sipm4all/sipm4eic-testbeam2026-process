#include <Math/Functor.h>
#include <Fit/Fitter.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "data_word.h"

struct channel_key_t {
  int device;
  int fifo;
  int column;
  int pixel;

  bool operator<(const channel_key_t &rhs) const
  {
    if (device != rhs.device) return device < rhs.device;
    if (fifo != rhs.fifo) return fifo < rhs.fifo;
    if (column != rhs.column) return column < rhs.column;
    return pixel < rhs.pixel;
  }
};

struct dcalib_pair_t {
  int tdc;
  int fine;
  int prev_tdc;
  int prev_fine;
  double coarse_delta;
};

struct fit_result_t {
  double par[9] = {0.};
  double err[9] = {0.};
  bool ok = false;
};

const int default_expected_period = 10000;
const int default_min_pairs = 100;
const int safe_margin = 10;

static std::string
channel_name(const channel_key_t &key)
{
  return "device" + std::to_string(key.device) +
         "_fifo" + std::to_string(key.fifo) +
         "_col" + std::to_string(key.column) +
         "_pix" + std::to_string(key.pixel);
}

static std::map<channel_key_t, std::vector<dcalib_pair_t>>
load_pairs(const std::string &infilename)
{
  auto *fin = TFile::Open(infilename.c_str());
  if (!fin || fin->IsZombie())
    throw std::runtime_error("cannot open input file: " + infilename);

  auto *tin = dynamic_cast<TTree *>(fin->Get("alcor"));
  if (!tin)
    throw std::runtime_error("cannot find TTree named 'alcor'");

  const auto nev = tin->GetEntries();
  data_t hit{};
  hit.link_to_tree(tin);

  std::map<channel_key_t, std::vector<dcalib_pair_t>> pairs;
  std::map<channel_key_t, data_t> previous;

  for (Long64_t iev = 0; iev < nev; ++iev) {
    auto bytes = tin->GetEntry(iev);
    if (bytes <= 0)
      throw std::runtime_error("ROOT GetEntry failed at entry " + std::to_string(iev));

    if (hit.is_start_spill() || hit.is_end_spill()) {
      previous.clear();
      continue;
    }

    if (!hit.is_alcor_hit())
      continue;

    if (hit.tdc < 0 || hit.tdc >= 4) {
      std::cerr << "WARNING: ignoring hit with invalid TDC"
                << " device=" << hit.device
                << " fifo=" << hit.fifo
                << " column=" << hit.column
                << " pixel=" << hit.pixel
                << " tdc=" << hit.tdc << std::endl;
      continue;
    }

    channel_key_t key{hit.device, hit.fifo, hit.column, hit.pixel};
    auto prev = previous.find(key);
    if (prev != previous.end()) {
      const int coarse = hit.coarse + data_t::rollover_to_clock * hit.rollover;
      const int prev_coarse = prev->second.coarse + data_t::rollover_to_clock * prev->second.rollover;

      dcalib_pair_t pair{};
      pair.tdc = hit.tdc;
      pair.fine = hit.fine;
      pair.prev_tdc = prev->second.tdc;
      pair.prev_fine = prev->second.fine;
      pair.coarse_delta = (double)(coarse - prev_coarse);
      pairs[key].push_back(pair);
    }

    previous[key] = hit;
  }

  fin->Close();
  return pairs;
}

static fit_result_t
fit_channel(const std::vector<dcalib_pair_t> &pairs, int expected_period)
{
  auto chi2 = [&](const double *par) {
    double f = 0.;
    for (const auto &p : pairs) {
      if (std::fabs(p.coarse_delta - (double)expected_period) >= (double)safe_margin)
        continue;

      const double corr = par[p.tdc] + (double)p.fine * par[p.tdc + 4];
      const double prev_corr = par[p.prev_tdc] + (double)p.prev_fine * par[p.prev_tdc + 4];
      const double delta = p.coarse_delta - (corr - prev_corr) - par[8];
      f += delta * delta;
    }
    return f;
  };

  ROOT::Math::Functor fcn(chi2, 9);
  ROOT::Fit::Fitter fitter;
  double pStart[9] = {0., 0., 0., 0., 0.01575, 0.01575, 0.01575, 0.01575,
                      (double)expected_period};
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("off_0");
  fitter.Config().ParSettings(1).SetName("off_1");
  fitter.Config().ParSettings(2).SetName("off_2");
  fitter.Config().ParSettings(3).SetName("off_3");
  fitter.Config().ParSettings(4).SetName("iif_0");
  fitter.Config().ParSettings(5).SetName("iif_1");
  fitter.Config().ParSettings(6).SetName("iif_2");
  fitter.Config().ParSettings(7).SetName("iif_3");
  fitter.Config().ParSettings(8).SetName("period");

  fitter.Config().ParSettings(0).Fix();
  fitter.Config().ParSettings(1).SetLimits(-0.32, 0.32);
  fitter.Config().ParSettings(2).SetLimits(-0.32, 0.32);
  fitter.Config().ParSettings(3).SetLimits(-0.32, 0.32);
  fitter.Config().ParSettings(4).SetLimits(0.0144, 0.0176);
  fitter.Config().ParSettings(5).SetLimits(0.0144, 0.0176);
  fitter.Config().ParSettings(6).SetLimits(0.0144, 0.0176);
  fitter.Config().ParSettings(7).SetLimits(0.0144, 0.0176);

  fit_result_t out;
  out.ok = fitter.FitFCN();
  auto result = fitter.Result();
  for (int ipar = 0; ipar < 9; ++ipar) {
    out.par[ipar] = result.Parameter(ipar);
    out.err[ipar] = result.ParError(ipar);
  }
  return out;
}

static void
write_histograms(const channel_key_t &key, const std::vector<dcalib_pair_t> &pairs,
                 const fit_result_t &fit, int expected_period, TFile *fout)
{
  auto *dir = fout->mkdir(channel_name(key).c_str());
  dir->cd();

  auto *hParam = new TH1F("hParam", "", 9, 0., 9.);
  for (int ipar = 0; ipar < 9; ++ipar) {
    hParam->SetBinContent(ipar + 1, fit.par[ipar]);
    hParam->SetBinError(ipar + 1, fit.err[ipar]);
  }

  auto *hDelta = new TH1F("hDelta", "", 1250, -10., 10.);
  for (const auto &p : pairs) {
    const double corr = fit.par[p.tdc] + (double)p.fine * fit.par[p.tdc + 4];
    const double prev_corr = fit.par[p.prev_tdc] + (double)p.prev_fine * fit.par[p.prev_tdc + 4];
    const double time_delta = p.coarse_delta - (corr - prev_corr);
    const double delta = time_delta - fit.par[8];
    hDelta->Fill(delta);
  }

  hParam->Write();
  hDelta->Write();
  fout->cd();
}

static bool
dcalib(const std::string &infilename, const std::string &outfilename,
       const std::string &calibfilename, int expected_period, int min_pairs)
{
  auto pairs = load_pairs(infilename);
  if (pairs.empty())
    throw std::runtime_error("no usable ALCOR hit pairs found");

  auto *fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie())
    throw std::runtime_error("cannot open output file: " + outfilename);

  std::ofstream calib(calibfilename);
  if (!calib)
    throw std::runtime_error("cannot open calibration output file: " + calibfilename);

  calib << "[TDC]\n";
  calib << "# device fifo column pixel tdc off iif\n";
  calib << "# generated by dcalib from " << infilename << "\n\n";
  calib << std::setprecision(10);

  int nfound = pairs.size();
  int nfitted = 0;
  int nskipped = 0;
  int nrows = 0;

  for (auto &item : pairs) {
    auto &key = item.first;
    auto &channel_pairs = item.second;

    if ((int)channel_pairs.size() < min_pairs) {
      std::cerr << "WARNING: skipping " << channel_name(key)
                << ": only " << channel_pairs.size()
                << " pairs, min-pairs=" << min_pairs << std::endl;
      ++nskipped;
      continue;
    }

    auto fit = fit_channel(channel_pairs, expected_period);
    if (!fit.ok)
      std::cerr << "WARNING: fit did not converge cleanly for " << channel_name(key) << std::endl;

    write_histograms(key, channel_pairs, fit, expected_period, fout);

    for (int tdc = 0; tdc < 4; ++tdc) {
      calib << key.device << ' '
            << key.fifo << ' '
            << key.column << ' '
            << key.pixel << ' '
            << tdc << ' '
            << fit.par[tdc] << ' '
            << fit.par[tdc + 4] << '\n';
      ++nrows;
    }
    calib << '\n';
    ++nfitted;
  }

  fout->Write();
  fout->Close();
  calib.close();

  std::cout << "channels found:       " << nfound << std::endl;
  std::cout << "channels fitted:      " << nfitted << std::endl;
  std::cout << "channels skipped:     " << nskipped << std::endl;
  std::cout << "calibration rows:     " << nrows << std::endl;
  return true;
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  int expected_period = default_expected_period;
  int min_pairs = default_min_pairs;
  std::string input_filename;
  std::string output_filename;
  std::string calib_filename;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print help message")
    ("input,i", po::value<std::string>(&input_filename)->required(), "Input ROOT data file")
    ("output,o", po::value<std::string>(&output_filename)->required(), "Output ROOT diagnostic file")
    ("calibration-output,c", po::value<std::string>(&calib_filename)->required(), "Output TDC calibration text file")
    ("period,p", po::value<int>(&expected_period)->default_value(default_expected_period), "Expected period in coarse clock cycles")
    ("min-pairs,m", po::value<int>(&min_pairs)->default_value(default_min_pairs), "Minimum adjacent-hit pairs required per channel");

  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cout << "usage: " << argv[0]
                << " --input input.root --output dcalib.root --calibration-output tdc.conf\n\n"
                << desc << std::endl;
      return 0;
    }

    po::notify(vm);

    if (expected_period <= 0)
      throw std::runtime_error("period must be positive");
    if (min_pairs <= 0)
      throw std::runtime_error("min-pairs must be positive");

    return dcalib(input_filename, output_filename, calib_filename,
                  expected_period, min_pairs) ? 0 : 1;
  } catch (const std::exception &ex) {
    std::cerr << "error: " << ex.what() << std::endl;
    std::cerr << desc << std::endl;
    return 1;
  }
}
