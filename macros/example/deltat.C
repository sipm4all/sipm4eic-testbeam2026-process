#include <TFile.h>
#include <TH2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <iostream>
#include <string>

struct category_reader_t {
  TTreeReaderArray<int> *frame_start;
  TTreeReaderArray<int> *frame_nhits;
  TTreeReaderArray<int> *device;
  TTreeReaderArray<int> *fifo;
  TTreeReaderArray<int> *type;
  TTreeReaderArray<int> *column;
  TTreeReaderArray<int> *pixel;
  TTreeReaderArray<int> *rollover;
  TTreeReaderArray<int> *coarse;
  TTreeReaderArray<int> *fine;

  category_reader_t(TTreeReader &reader, const std::string &prefix)
    : frame_start(new TTreeReaderArray<int>(reader, (prefix + "_frame_start").c_str())),
      frame_nhits(new TTreeReaderArray<int>(reader, (prefix + "_frame_nhits").c_str())),
      device(new TTreeReaderArray<int>(reader, (prefix + "_device").c_str())),
      fifo(new TTreeReaderArray<int>(reader, (prefix + "_fifo").c_str())),
      type(new TTreeReaderArray<int>(reader, (prefix + "_type").c_str())),
      column(new TTreeReaderArray<int>(reader, (prefix + "_column").c_str())),
      pixel(new TTreeReaderArray<int>(reader, (prefix + "_pixel").c_str())),
      rollover(new TTreeReaderArray<int>(reader, (prefix + "_rollover").c_str())),
      coarse(new TTreeReaderArray<int>(reader, (prefix + "_coarse").c_str())),
      fine(new TTreeReaderArray<int>(reader, (prefix + "_fine").c_str()))
  {
  }

  double time(int i) const
  {
    auto out = (*coarse)[i] + 32768. * (*rollover)[i];
    if ((*type)[i] == 1)
      out -= 0.0157 * (*fine)[i];
    return out;
  }

  int channel(int i) const
  {
    int chip = (*fifo)[i] / 4;
    int local_channel = (*pixel)[i] + 4 * (*column)[i] + 32 * chip;
    return local_channel + 256 * ((*device)[i] - 192);
  }

  bool match(int i, int trigger_type, int trigger_device, int trigger_fifo,
             int trigger_column, int trigger_pixel) const
  {
    return (*type)[i] == trigger_type &&
           (*device)[i] == trigger_device &&
           (*fifo)[i] == trigger_fifo &&
           (*column)[i] == trigger_column &&
           (*pixel)[i] == trigger_pixel;
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

        auto event_time = trigger_cat->time(itrigger);
        ++ntriggers;

        for (auto cat : cats) {
          int first = (*cat->frame_start)[iframe];
          int nhits = (*cat->frame_nhits)[iframe];
          for (int ihit = 0; ihit < nhits; ++ihit) {
            int i = first + ihit;
            if (cat == trigger_cat && i == itrigger)
              continue;

            auto delta_t = cat->time(i) - event_time;
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
