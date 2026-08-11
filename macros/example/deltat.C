#include <TFile.h>
#include <TH2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

struct category_reader_t {
  TTreeReaderArray<int> *frame_start;
  TTreeReaderArray<int> *frame_nhits;
  TTreeReaderArray<int> *device;
  TTreeReaderArray<int> *fifo;
  TTreeReaderArray<int> *type;
  TTreeReaderArray<int> *column;
  TTreeReaderArray<int> *pixel;
  TTreeReaderArray<int> *tdc;
  TTreeReaderArray<int> *rollover;
  TTreeReaderArray<int> *coarse;
  TTreeReaderArray<int> *fine;
  TTreeReaderArray<double> *time;

  category_reader_t(TTreeReader &reader, const std::string &prefix)
    : frame_start(new TTreeReaderArray<int>(reader, (prefix + "_frame_start").c_str())),
      frame_nhits(new TTreeReaderArray<int>(reader, (prefix + "_frame_nhits").c_str())),
      device(new TTreeReaderArray<int>(reader, (prefix + "_device").c_str())),
      fifo(new TTreeReaderArray<int>(reader, (prefix + "_fifo").c_str())),
      type(new TTreeReaderArray<int>(reader, (prefix + "_type").c_str())),
      column(new TTreeReaderArray<int>(reader, (prefix + "_column").c_str())),
      pixel(new TTreeReaderArray<int>(reader, (prefix + "_pixel").c_str())),
      tdc(new TTreeReaderArray<int>(reader, (prefix + "_tdc").c_str())),
      rollover(new TTreeReaderArray<int>(reader, (prefix + "_rollover").c_str())),
      coarse(new TTreeReaderArray<int>(reader, (prefix + "_coarse").c_str())),
      fine(new TTreeReaderArray<int>(reader, (prefix + "_fine").c_str())),
      time(new TTreeReaderArray<double>(reader, (prefix + "_time").c_str()))
  {
  }

  double hit_time(int i) const
  {
    return (*time)[i];
  }

  int channel(int i) const
  {
    int chip = (*fifo)[i] / 4;
    int local_channel = (*pixel)[i] + 4 * (*column)[i] + 32 * chip;
    return local_channel + 256 * ((*device)[i] - 192);
  }

  static bool match_field(int value, int requested)
  {
    return requested < 0 || value == requested;
  }

  bool match(int i, int trigger_type, int trigger_device, int trigger_fifo,
             int trigger_column, int trigger_pixel) const
  {
    return match_field((*type)[i], trigger_type) &&
           match_field((*device)[i], trigger_device) &&
           match_field((*fifo)[i], trigger_fifo) &&
           match_field((*column)[i], trigger_column) &&
           match_field((*pixel)[i], trigger_pixel);
  }
};

void
deltat(const std::string filename,
       int trigger_type,
       int trigger_device,
       int trigger_fifo,
       int trigger_column,
       int trigger_pixel,
       const std::string outfilename = "deltat.root")
{
  auto fin = TFile::Open(filename.c_str());
  if (!fin || fin->IsZombie()) {
    std::cerr << " --- could not open input file: " << filename << std::endl;
    return;
  }

  auto tin = (TTree *)fin->Get("frames");
  if (!tin) {
    std::cerr << " --- could not find 'frames' tree in input file" << std::endl;
    fin->Close();
    return;
  }

  auto scan = [&](TH2D *hist, double &dtmin, double &dtmax, bool &have_dt, int &ntriggers, int &nfills) {
    TTreeReader reader(tin);
    TTreeReaderValue<int> nframes(reader, "nframes");
    category_reader_t trigger(reader, "trigger");
    category_reader_t timing(reader, "timing");
    category_reader_t cherenkov(reader, "cherenkov");
    category_reader_t *cats[] = {&trigger, &timing, &cherenkov};

    while (reader.Next()) {
      for (int iframe = 0; iframe < *nframes; ++iframe) {
        category_reader_t *trigger_cat = nullptr;
        int itrigger = -1;
        for (auto cat : cats) {
          int first = (*cat->frame_start)[iframe];
          int nhits = (*cat->frame_nhits)[iframe];
          for (int ihit = 0; ihit < nhits; ++ihit) {
            int i = first + ihit;
            if (!cat->match(i, trigger_type, trigger_device, trigger_fifo, trigger_column, trigger_pixel))
              continue;
            trigger_cat = cat;
            itrigger = i;
            break;
          }
          if (itrigger >= 0)
            break;
        }

        if (itrigger < 0)
          continue;

        auto event_time = trigger_cat->hit_time(itrigger);
        ++ntriggers;

        for (auto cat : cats) {
          int first = (*cat->frame_start)[iframe];
          int nhits = (*cat->frame_nhits)[iframe];
          for (int ihit = 0; ihit < nhits; ++ihit) {
            int i = first + ihit;
            if (cat == trigger_cat && i == itrigger)
              continue;

            auto delta_t = cat->hit_time(i) - event_time;
            if (hist) {
              hist->Fill(cat->channel(i), delta_t);
              ++nfills;
            } else if (!have_dt) {
              dtmin = dtmax = delta_t;
              have_dt = true;
            } else {
              if (delta_t < dtmin) dtmin = delta_t;
              if (delta_t > dtmax) dtmax = delta_t;
            }
          }
        }
      }
    }
  };

  double dtmin = 0.;
  double dtmax = 0.;
  bool have_dt = false;
  int ntriggers = 0;
  int nfills = 0;
  scan(nullptr, dtmin, dtmax, have_dt, ntriggers, nfills);

  if (!have_dt) {
    dtmin = -1.;
    dtmax = 1.;
  } else if (dtmin == dtmax) {
    dtmin -= 1.;
    dtmax += 1.;
  } else {
    auto margin = 0.001 * (dtmax - dtmin);
    dtmin -= margin;
    dtmax += margin;
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    fin->Close();
    return;
  }

  constexpr int min_device = 192;
  constexpr int max_device = 200;
  constexpr int channels_per_device = 256;
  constexpr int nchannels = (max_device - min_device + 1) * channels_per_device;
  int nbins = std::max(1, (int)std::ceil(16. * (dtmax - dtmin)));
  auto hDeltaT = new TH2D("hDeltaT", "", nchannels, 0., nchannels, nbins, dtmin, dtmax);

  ntriggers = 0;
  nfills = 0;
  scan(hDeltaT, dtmin, dtmax, have_dt, ntriggers, nfills);

  if (ntriggers == 0)
    std::cerr << " --- no matching trigger hit found inside the stored frames" << std::endl;
  std::cout << " --- triggers found: " << ntriggers << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;

  hDeltaT->Sumw2();
  if (ntriggers > 0)
    hDeltaT->Scale(1. / ntriggers);

  hDeltaT->Write();
  fout->Close();
  fin->Close();
}


struct hit_ref_t {
  category_reader_t *cat;
  int index;
};

void
deltat(const std::string filename,
       int target_type,
       int target_device,
       int target_fifo,
       int target_column,
       int target_pixel,
       int reference_type,
       int reference_device,
       int reference_fifo,
       int reference_column,
       int reference_pixel,
       const std::string outfilename = "deltat.root")
{
  auto fin = TFile::Open(filename.c_str());
  if (!fin || fin->IsZombie()) {
    std::cerr << " --- could not open input file: " << filename << std::endl;
    return;
  }

  auto tin = (TTree *)fin->Get("frames");
  if (!tin) {
    std::cerr << " --- could not find 'frames' tree in input file" << std::endl;
    fin->Close();
    return;
  }

  auto scan = [&](TH2D *hist, TH2D **hist_tdc, double &dtmin, double &dtmax,
                  bool &have_dt, int &nframes_used, int &nfills) {
    TTreeReader reader(tin);
    TTreeReaderValue<int> nframes(reader, "nframes");
    category_reader_t trigger(reader, "trigger");
    category_reader_t timing(reader, "timing");
    category_reader_t cherenkov(reader, "cherenkov");
    category_reader_t *cats[] = {&trigger, &timing, &cherenkov};

    while (reader.Next()) {
      for (int iframe = 0; iframe < *nframes; ++iframe) {
        std::vector<hit_ref_t> targets;
        std::vector<hit_ref_t> references;

        for (auto cat : cats) {
          int first = (*cat->frame_start)[iframe];
          int nhits = (*cat->frame_nhits)[iframe];
          for (int ihit = 0; ihit < nhits; ++ihit) {
            int i = first + ihit;
            if (cat->match(i, target_type, target_device, target_fifo, target_column, target_pixel))
              targets.push_back({cat, i});
            if (cat->match(i, reference_type, reference_device, reference_fifo, reference_column, reference_pixel))
              references.push_back({cat, i});
          }
        }

        if (targets.empty() || references.empty())
          continue;

        ++nframes_used;
        for (const auto &target : targets) {
          for (const auto &reference : references) {
            if (target.cat == reference.cat && target.index == reference.index)
              continue;

            auto delta_t = target.cat->hit_time(target.index) - reference.cat->hit_time(reference.index);
            if (hist) {
              int channel = target.cat->channel(target.index);
              hist->Fill(channel, delta_t);
              int tdc = (*target.cat->tdc)[target.index];
              if (tdc >= 0 && tdc < 4 && hist_tdc && hist_tdc[tdc])
                hist_tdc[tdc]->Fill(channel, delta_t);
              ++nfills;
            } else if (!have_dt) {
              dtmin = dtmax = delta_t;
              have_dt = true;
            } else {
              if (delta_t < dtmin) dtmin = delta_t;
              if (delta_t > dtmax) dtmax = delta_t;
            }
          }
        }
      }
    }
  };

  double dtmin = 0.;
  double dtmax = 0.;
  bool have_dt = false;
  int nframes_used = 0;
  int nfills = 0;
  scan(nullptr, nullptr, dtmin, dtmax, have_dt, nframes_used, nfills);

  if (!have_dt) {
    dtmin = -1.;
    dtmax = 1.;
  } else if (dtmin == dtmax) {
    dtmin -= 1.;
    dtmax += 1.;
  } else {
    auto margin = 0.001 * (dtmax - dtmin);
    dtmin -= margin;
    dtmax += margin;
  }

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    fin->Close();
    return;
  }

  constexpr int min_device = 192;
  constexpr int max_device = 200;
  constexpr int channels_per_device = 256;
  constexpr int nchannels = (max_device - min_device + 1) * channels_per_device;
  int nbins = std::max(1, (int)std::ceil(16. * (dtmax - dtmin)));
  auto hDeltaT = new TH2D("hDeltaT", "", nchannels, 0., nchannels, nbins, dtmin, dtmax);

  TH2D *hDeltaT_tdc[4] = {nullptr, nullptr, nullptr, nullptr};
  for (int itdc = 0; itdc < 4; ++itdc)
    hDeltaT_tdc[itdc] = new TH2D(Form("hDeltaT_tdc%d", itdc), "", nchannels, 0., nchannels, nbins, dtmin, dtmax);

  nframes_used = 0;
  nfills = 0;
  scan(hDeltaT, hDeltaT_tdc, dtmin, dtmax, have_dt, nframes_used, nfills);

  if (nframes_used == 0)
    std::cerr << " --- no frames with both target and reference hits found" << std::endl;
  std::cout << " --- frames with target/reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;

  hDeltaT->Sumw2();
  if (nframes_used > 0)
    hDeltaT->Scale(1. / nframes_used);
  hDeltaT->Write();

  for (int itdc = 0; itdc < 4; ++itdc) {
    hDeltaT_tdc[itdc]->Sumw2();
    if (nframes_used > 0)
      hDeltaT_tdc[itdc]->Scale(1. / nframes_used);
    hDeltaT_tdc[itdc]->Write();
  }

  fout->Close();
  fin->Close();
}
