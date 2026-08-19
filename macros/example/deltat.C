#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

constexpr int deltat_nbins = 2048;
constexpr double deltat_min = -32.;
constexpr double deltat_max = 32.;
constexpr int spill_nbins = 100;
constexpr double spill_min = 0.;
constexpr double spill_max = 100.;

struct field_selector_t {
  int type;
  int device;
  int fifo;
  int column;
  int pixel;

  field_selector_t(int type_, int device_, int fifo_, int column_, int pixel_)
    : type(type_), device(device_), fifo(fifo_), column(column_), pixel(pixel_)
  {
  }
};

struct channel_selector_t {
  int type;
  int channel;

  channel_selector_t(int type_, int channel_)
    : type(type_), channel(channel_)
  {
  }
};

struct timing_reference_t {
  std::string branch;

  timing_reference_t(const std::string &branch_ = "T")
    : branch(branch_)
  {
  }
};

struct ring_selection_t {
  int min_inliers;
  int max_inliers;
  int iterations;
  double tolerance;
  double min_radius;
  double max_radius;

  ring_selection_t(int min_inliers_ = 8,
                   int max_inliers_ = std::numeric_limits<int>::max(),
                   int iterations_ = 256,
                   double tolerance_ = 3.5,
                   double min_radius_ = 1.,
                   double max_radius_ = 200.)
    : min_inliers(min_inliers_), max_inliers(max_inliers_), iterations(iterations_),
      tolerance(tolerance_),
      min_radius(min_radius_), max_radius(max_radius_)
  {
  }
};

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
  TTreeReaderArray<double> *x;
  TTreeReaderArray<double> *y;

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
      time(new TTreeReaderArray<double>(reader, (prefix + "_time").c_str())),
      x(nullptr),
      y(nullptr)
  {
    if (reader.GetTree()->GetBranch((prefix + "_x").c_str()))
      x = new TTreeReaderArray<double>(reader, (prefix + "_x").c_str());
    if (reader.GetTree()->GetBranch((prefix + "_y").c_str()))
      y = new TTreeReaderArray<double>(reader, (prefix + "_y").c_str());
  }

  double hit_time(int i) const
  {
    return (*time)[i];
  }

  bool has_xy() const
  {
    return x && y;
  }

  double hit_x(int i) const
  {
    return (*x)[i];
  }

  double hit_y(int i) const
  {
    return (*y)[i];
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

  bool match(int i, field_selector_t selector) const
  {
    return match_field((*type)[i], selector.type) &&
           match_field((*device)[i], selector.device) &&
           match_field((*fifo)[i], selector.fifo) &&
           match_field((*column)[i], selector.column) &&
           match_field((*pixel)[i], selector.pixel);
  }

  bool match(int i, channel_selector_t selector) const
  {
    int requested_type = selector.type;
    int requested_channel = selector.channel;

    if (requested_type == 9)
      return match(i, field_selector_t{9, requested_channel, 32, -1, -1});

    if (!match_field((*type)[i], requested_type))
      return false;

    if ((*type)[i] != 1)
      return false;

    return match_field(channel(i), requested_channel);
  }
};

double
timing_reference_value(const std::string &branch,
                       int iframe,
                       const TTreeReaderArray<double> &T0,
                       const TTreeReaderArray<double> &T1,
                       const TTreeReaderArray<double> &T)
{
  if (branch == "T0")
    return T0[iframe];
  if (branch == "T1")
    return T1[iframe];
  if (branch == "T")
    return T[iframe];

  std::cerr << " --- unsupported timing reference branch: " << branch
            << " (valid choices are T, T0, T1)" << std::endl;
  return std::numeric_limits<double>::quiet_NaN();
}

void
deltat(const std::string filename,
       field_selector_t reference,
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

  auto scan = [&](TH2D *hist, int &ntriggers, int &nfills) {
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
            if (!cat->match(i, reference))
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
            }
          }
        }
      }
    }
  };

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
  auto hDeltaT = new TH2D("hDeltaT", "", nchannels, 0., nchannels, deltat_nbins, deltat_min, deltat_max);

  int ntriggers = 0;
  int nfills = 0;
  scan(hDeltaT, ntriggers, nfills);

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

struct ring_point_t {
  double x;
  double y;
};

struct ring_t {
  double x = 0.;
  double y = 0.;
  double radius = 0.;
  int inliers = 0;
  double residual = std::numeric_limits<double>::infinity();
  bool valid = false;

  bool contains(double px, double py, double tolerance) const
  {
    if (!valid)
      return false;
    return std::abs(std::hypot(px - x, py - y) - radius) <= tolerance;
  }
};

bool
circle_from_three(const ring_point_t &a,
                  const ring_point_t &b,
                  const ring_point_t &c,
                  ring_t &ring)
{
  double determinant = 2. * (a.x * (b.y - c.y) +
                             b.x * (c.y - a.y) +
                             c.x * (a.y - b.y));
  if (std::abs(determinant) < 1.e-9)
    return false;

  double aa = a.x * a.x + a.y * a.y;
  double bb = b.x * b.x + b.y * b.y;
  double cc = c.x * c.x + c.y * c.y;
  ring.x = (aa * (b.y - c.y) + bb * (c.y - a.y) + cc * (a.y - b.y)) / determinant;
  ring.y = (aa * (c.x - b.x) + bb * (a.x - c.x) + cc * (b.x - a.x)) / determinant;
  ring.radius = std::hypot(a.x - ring.x, a.y - ring.y);
  return std::isfinite(ring.x) && std::isfinite(ring.y) && std::isfinite(ring.radius);
}

ring_t
fit_ring(const category_reader_t &cherenkov,
         int iframe,
         int spill,
         const ring_selection_t &selection)
{
  ring_t best;
  if (!cherenkov.has_xy())
    return best;

  std::vector<ring_point_t> points;
  int first = (*cherenkov.frame_start)[iframe];
  int nhits = (*cherenkov.frame_nhits)[iframe];
  points.reserve(nhits);
  for (int ihit = 0; ihit < nhits; ++ihit) {
    int index = first + ihit;
    double x = cherenkov.hit_x(index);
    double y = cherenkov.hit_y(index);
    if (std::isfinite(x) && std::isfinite(y))
      points.push_back({x, y});
  }

  if (points.size() < 3)
    return best;

  std::mt19937 generator(static_cast<unsigned>(0x9e3779b9u ^
                                                static_cast<unsigned>(spill * 65537 + iframe)));
  std::uniform_int_distribution<int> pick(0, points.size() - 1);
  for (int iteration = 0; iteration < selection.iterations; ++iteration) {
    int ia = pick(generator);
    int ib = pick(generator);
    int ic = pick(generator);
    if (ia == ib || ia == ic || ib == ic)
      continue;

    ring_t candidate;
    if (!circle_from_three(points[ia], points[ib], points[ic], candidate))
      continue;
    if (candidate.radius < selection.min_radius || candidate.radius > selection.max_radius)
      continue;

    int inliers = 0;
    double residual = 0.;
    for (const auto &point : points) {
      double distance = std::abs(std::hypot(point.x - candidate.x,
                                             point.y - candidate.y) - candidate.radius);
      if (distance <= selection.tolerance) {
        ++inliers;
        residual += distance;
      }
    }

    if (inliers > best.inliers ||
        (inliers == best.inliers && residual < best.residual)) {
      candidate.inliers = inliers;
      candidate.residual = residual;
      best = candidate;
    }
  }

  best.valid = best.inliers >= selection.min_inliers &&
               best.inliers <= selection.max_inliers;
  return best;
}

void
deltat(const std::string filename,
       field_selector_t target_selector,
       timing_reference_t timing_reference,
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

  if (!tin->GetBranch("timing_valid") ||
      !tin->GetBranch("T0") ||
      !tin->GetBranch("T1") ||
      !tin->GetBranch("T")) {
    std::cerr << " --- input file does not contain timing estimator branches" << std::endl;
    fin->Close();
    return;
  }

  if (timing_reference.branch != "T" &&
      timing_reference.branch != "T0" &&
      timing_reference.branch != "T1") {
    std::cerr << " --- unsupported timing reference branch: " << timing_reference.branch
              << " (valid choices are T, T0, T1)" << std::endl;
    fin->Close();
    return;
  }

  auto scan = [&](TH2D *hist, TH2D **hist_tdc, TH2D *hist_spill, int &nframes_used, int &nfills) {
    TTreeReader reader(tin);
    TTreeReaderValue<int> spill_id(reader, "id");
    TTreeReaderValue<int> nframes(reader, "nframes");
    TTreeReaderArray<int> timing_valid(reader, "timing_valid");
    TTreeReaderArray<double> T0(reader, "T0");
    TTreeReaderArray<double> T1(reader, "T1");
    TTreeReaderArray<double> T(reader, "T");
    category_reader_t trigger(reader, "trigger");
    category_reader_t timing(reader, "timing");
    category_reader_t cherenkov(reader, "cherenkov");
    category_reader_t *cats[] = {&trigger, &timing, &cherenkov};

    while (reader.Next()) {
      for (int iframe = 0; iframe < *nframes; ++iframe) {
        if (timing_valid[iframe] == 0)
          continue;

        auto reference_time = timing_reference_value(timing_reference.branch, iframe, T0, T1, T);
        if (!std::isfinite(reference_time))
          continue;

        std::vector<hit_ref_t> targets;
        for (auto cat : cats) {
          int first = (*cat->frame_start)[iframe];
          int nhits = (*cat->frame_nhits)[iframe];
          for (int ihit = 0; ihit < nhits; ++ihit) {
            int i = first + ihit;
            if (cat->match(i, target_selector))
              targets.push_back({cat, i});
          }
        }

        if (targets.empty())
          continue;

        ++nframes_used;
        for (const auto &target : targets) {
          auto delta_t = target.cat->hit_time(target.index) - reference_time;
          if (hist) {
            int channel = target.cat->channel(target.index);
            hist->Fill(channel, delta_t);
            if (hist_spill)
              hist_spill->Fill(*spill_id, delta_t);
            int tdc = (*target.cat->tdc)[target.index];
            if (tdc >= 0 && tdc < 4 && hist_tdc && hist_tdc[tdc])
              hist_tdc[tdc]->Fill((*target.cat->fine)[target.index], delta_t);
            ++nfills;
          }
        }
      }
    }
  };

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
  auto hDeltaT = new TH2D("hDeltaT", "", nchannels, 0., nchannels, deltat_nbins, deltat_min, deltat_max);
  auto hDeltaT_spill = new TH2D("hDeltaT_spill", "", spill_nbins, spill_min, spill_max, deltat_nbins, deltat_min, deltat_max);

  TH2D *hDeltaT_tdc[4] = {nullptr, nullptr, nullptr, nullptr};
  for (int itdc = 0; itdc < 4; ++itdc)
    hDeltaT_tdc[itdc] = new TH2D(Form("hDeltaT_tdc%d", itdc), "", 256, 0., 256., deltat_nbins, deltat_min, deltat_max);

  int nframes_used = 0;
  int nfills = 0;
  scan(hDeltaT, hDeltaT_tdc, hDeltaT_spill, nframes_used, nfills);

  if (nframes_used == 0)
    std::cerr << " --- no valid timing-reference frames with target hits found" << std::endl;
  std::cout << " --- timing reference: " << timing_reference.branch << std::endl;
  std::cout << " --- frames with target/timing-reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;

  hDeltaT->Sumw2();
  if (nframes_used > 0)
    hDeltaT->Scale(1. / nframes_used);
  hDeltaT->Write();

  hDeltaT_spill->Sumw2();
  if (nframes_used > 0)
    hDeltaT_spill->Scale(1. / nframes_used);
  hDeltaT_spill->Write();

  for (int itdc = 0; itdc < 4; ++itdc) {
    hDeltaT_tdc[itdc]->Sumw2();
    if (nframes_used > 0)
      hDeltaT_tdc[itdc]->Scale(1. / nframes_used);
    hDeltaT_tdc[itdc]->Write();
  }

  fout->Close();
  fin->Close();
}

void
deltat_timing_reference_impl(const std::string filename,
                             channel_selector_t target_selector,
                             timing_reference_t timing_reference,
                             const ring_selection_t *ring_selection,
                             const std::string outfilename)
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

  if (!tin->GetBranch("timing_valid") ||
      !tin->GetBranch("T0") ||
      !tin->GetBranch("T1") ||
      !tin->GetBranch("T")) {
    std::cerr << " --- input file does not contain timing estimator branches" << std::endl;
    fin->Close();
    return;
  }

  if (ring_selection &&
      (!tin->GetBranch("cherenkov_x") || !tin->GetBranch("cherenkov_y"))) {
    std::cerr << " --- ring selection requires cherenkov_x and cherenkov_y branches" << std::endl;
    fin->Close();
    return;
  }

  if (timing_reference.branch != "T" &&
      timing_reference.branch != "T0" &&
      timing_reference.branch != "T1") {
    std::cerr << " --- unsupported timing reference branch: " << timing_reference.branch
              << " (valid choices are T, T0, T1)" << std::endl;
    fin->Close();
    return;
  }

  auto scan = [&](TH2D *hist, TH2D **hist_tdc, TH2D *hist_spill,
                  TH1D *hist_ring_radius, TH1D *hist_ring_inliers,
                  int &nframes_used, int &nfills, int &nring_frames) {
    TTreeReader reader(tin);
    TTreeReaderValue<int> spill_id(reader, "id");
    TTreeReaderValue<int> nframes(reader, "nframes");
    TTreeReaderArray<int> timing_valid(reader, "timing_valid");
    TTreeReaderArray<double> T0(reader, "T0");
    TTreeReaderArray<double> T1(reader, "T1");
    TTreeReaderArray<double> T(reader, "T");
    category_reader_t trigger(reader, "trigger");
    category_reader_t timing(reader, "timing");
    category_reader_t cherenkov(reader, "cherenkov");
    category_reader_t *cats[] = {&trigger, &timing, &cherenkov};

    while (reader.Next()) {
      for (int iframe = 0; iframe < *nframes; ++iframe) {
        if (timing_valid[iframe] == 0)
          continue;

        auto reference_time = timing_reference_value(timing_reference.branch, iframe, T0, T1, T);
        if (!std::isfinite(reference_time))
          continue;

        ring_t ring;
        if (ring_selection) {
          ring = fit_ring(cherenkov, iframe, *spill_id, *ring_selection);
          if (!ring.valid)
            continue;
          ++nring_frames;
          if (hist_ring_radius)
            hist_ring_radius->Fill(ring.radius);
          if (hist_ring_inliers)
            hist_ring_inliers->Fill(ring.inliers);
        }

        std::vector<hit_ref_t> targets;
        for (auto cat : cats) {
          int first = (*cat->frame_start)[iframe];
          int nhits = (*cat->frame_nhits)[iframe];
          for (int ihit = 0; ihit < nhits; ++ihit) {
            int i = first + ihit;
            if (cat->match(i, target_selector) &&
                (!ring_selection ||
                 (cat == &cherenkov && ring.contains(cat->hit_x(i), cat->hit_y(i),
                                                     ring_selection->tolerance))))
              targets.push_back({cat, i});
          }
        }

        if (targets.empty())
          continue;

        ++nframes_used;
        for (const auto &target : targets) {
          auto delta_t = target.cat->hit_time(target.index) - reference_time;
          if (hist) {
            int channel = target.cat->channel(target.index);
            hist->Fill(channel, delta_t);
            if (hist_spill)
              hist_spill->Fill(*spill_id, delta_t);
            int tdc = (*target.cat->tdc)[target.index];
            if (tdc >= 0 && tdc < 4 && hist_tdc && hist_tdc[tdc])
              hist_tdc[tdc]->Fill((*target.cat->fine)[target.index], delta_t);
            ++nfills;
          }
        }
      }
    }
  };

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
  auto hDeltaT = new TH2D("hDeltaT", "", nchannels, 0., nchannels, deltat_nbins, deltat_min, deltat_max);
  auto hDeltaT_spill = new TH2D("hDeltaT_spill", "", spill_nbins, spill_min, spill_max, deltat_nbins, deltat_min, deltat_max);

  TH1D *hRingRadius = nullptr;
  TH1D *hRingInliers = nullptr;
  if (ring_selection) {
    hRingRadius = new TH1D("hRingRadius", "RANSAC ring radius;radius [mm];frames",
                           200, ring_selection->min_radius, ring_selection->max_radius);
    hRingInliers = new TH1D("hRingInliers", "RANSAC ring inliers;inlier hits;frames",
                            256, -0.5, 255.5);
  }

  TH2D *hDeltaT_tdc[4] = {nullptr, nullptr, nullptr, nullptr};
  for (int itdc = 0; itdc < 4; ++itdc)
    hDeltaT_tdc[itdc] = new TH2D(Form("hDeltaT_tdc%d", itdc), "", 256, 0., 256., deltat_nbins, deltat_min, deltat_max);

  int nframes_used = 0;
  int nfills = 0;
  int nring_frames = 0;
  scan(hDeltaT, hDeltaT_tdc, hDeltaT_spill, hRingRadius, hRingInliers,
       nframes_used, nfills, nring_frames);

  if (nframes_used == 0)
    std::cerr << " --- no valid timing-reference frames with target hits found" << std::endl;
  std::cout << " --- timing reference: " << timing_reference.branch << std::endl;
  std::cout << " --- frames with target/timing-reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;
  if (ring_selection)
    std::cout << " --- frames with accepted Cherenkov rings: " << nring_frames << std::endl;

  hDeltaT->Sumw2();
  if (nframes_used > 0)
    hDeltaT->Scale(1. / nframes_used);
  hDeltaT->Write();

  hDeltaT_spill->Sumw2();
  if (nframes_used > 0)
    hDeltaT_spill->Scale(1. / nframes_used);
  hDeltaT_spill->Write();

  for (int itdc = 0; itdc < 4; ++itdc) {
    hDeltaT_tdc[itdc]->Sumw2();
    if (nframes_used > 0)
      hDeltaT_tdc[itdc]->Scale(1. / nframes_used);
    hDeltaT_tdc[itdc]->Write();
  }

  if (hRingRadius)
    hRingRadius->Write();
  if (hRingInliers)
    hRingInliers->Write();

  fout->Close();
  fin->Close();
}

void
deltat(const std::string filename,
       channel_selector_t target_selector,
       timing_reference_t timing_reference,
       const std::string outfilename = "deltat.root")
{
  deltat_timing_reference_impl(filename, target_selector, timing_reference, nullptr, outfilename);
}

void
deltat(const std::string filename,
       channel_selector_t target_selector,
       timing_reference_t timing_reference,
       ring_selection_t ring_selection,
       const std::string outfilename = "deltat.root")
{
  deltat_timing_reference_impl(filename, target_selector, timing_reference,
                                &ring_selection, outfilename);
}

void
deltat(const std::string filename,
       field_selector_t target_selector,
       field_selector_t reference_selector,
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

  auto scan = [&](TH2D *hist, TH2D **hist_tdc, TH2D *hist_spill, int &nframes_used, int &nfills) {
    TTreeReader reader(tin);
    TTreeReaderValue<int> spill_id(reader, "id");
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
            if (cat->match(i, target_selector))
              targets.push_back({cat, i});
            if (cat->match(i, reference_selector))
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
              if (hist_spill)
                hist_spill->Fill(*spill_id, delta_t);
              int tdc = (*target.cat->tdc)[target.index];
              if (tdc >= 0 && tdc < 4 && hist_tdc && hist_tdc[tdc])
                hist_tdc[tdc]->Fill((*target.cat->fine)[target.index], delta_t);
              ++nfills;
            }
          }
        }
      }
    }
  };

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
  auto hDeltaT = new TH2D("hDeltaT", "", nchannels, 0., nchannels, deltat_nbins, deltat_min, deltat_max);
  auto hDeltaT_spill = new TH2D("hDeltaT_spill", "", spill_nbins, spill_min, spill_max, deltat_nbins, deltat_min, deltat_max);

  TH2D *hDeltaT_tdc[4] = {nullptr, nullptr, nullptr, nullptr};
  for (int itdc = 0; itdc < 4; ++itdc)
    hDeltaT_tdc[itdc] = new TH2D(Form("hDeltaT_tdc%d", itdc), "", 256, 0., 256., deltat_nbins, deltat_min, deltat_max);

  int nframes_used = 0;
  int nfills = 0;
  scan(hDeltaT, hDeltaT_tdc, hDeltaT_spill, nframes_used, nfills);

  if (nframes_used == 0)
    std::cerr << " --- no frames with both target and reference hits found" << std::endl;
  std::cout << " --- frames with target/reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;

  hDeltaT->Sumw2();
  if (nframes_used > 0)
    hDeltaT->Scale(1. / nframes_used);
  hDeltaT->Write();

  hDeltaT_spill->Sumw2();
  if (nframes_used > 0)
    hDeltaT_spill->Scale(1. / nframes_used);
  hDeltaT_spill->Write();

  for (int itdc = 0; itdc < 4; ++itdc) {
    hDeltaT_tdc[itdc]->Sumw2();
    if (nframes_used > 0)
      hDeltaT_tdc[itdc]->Scale(1. / nframes_used);
    hDeltaT_tdc[itdc]->Write();
  }

  fout->Close();
  fin->Close();
}


void
deltat(const std::string filename,
       channel_selector_t target_selector,
       channel_selector_t reference_selector,
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

  auto scan = [&](TH2D *hist, TH2D **hist_tdc, TH2D *hist_spill, int &nframes_used, int &nfills) {
    TTreeReader reader(tin);
    TTreeReaderValue<int> spill_id(reader, "id");
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
            if (cat->match(i, target_selector))
              targets.push_back({cat, i});
            if (cat->match(i, reference_selector))
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
              if (hist_spill)
                hist_spill->Fill(*spill_id, delta_t);
              int tdc = (*target.cat->tdc)[target.index];
              if (tdc >= 0 && tdc < 4 && hist_tdc && hist_tdc[tdc])
                hist_tdc[tdc]->Fill((*target.cat->fine)[target.index], delta_t);
              ++nfills;
            }
          }
        }
      }
    }
  };

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
  auto hDeltaT = new TH2D("hDeltaT", "", nchannels, 0., nchannels, deltat_nbins, deltat_min, deltat_max);
  auto hDeltaT_spill = new TH2D("hDeltaT_spill", "", spill_nbins, spill_min, spill_max, deltat_nbins, deltat_min, deltat_max);

  TH2D *hDeltaT_tdc[4] = {nullptr, nullptr, nullptr, nullptr};
  for (int itdc = 0; itdc < 4; ++itdc)
    hDeltaT_tdc[itdc] = new TH2D(Form("hDeltaT_tdc%d", itdc), "", 256, 0., 256., deltat_nbins, deltat_min, deltat_max);

  int nframes_used = 0;
  int nfills = 0;
  scan(hDeltaT, hDeltaT_tdc, hDeltaT_spill, nframes_used, nfills);

  if (nframes_used == 0)
    std::cerr << " --- no frames with both target and reference hits found" << std::endl;
  std::cout << " --- frames with target/reference hits: " << nframes_used << std::endl;
  std::cout << " --- histogram fills: " << nfills << std::endl;

  hDeltaT->Sumw2();
  if (nframes_used > 0)
    hDeltaT->Scale(1. / nframes_used);
  hDeltaT->Write();

  hDeltaT_spill->Sumw2();
  if (nframes_used > 0)
    hDeltaT_spill->Scale(1. / nframes_used);
  hDeltaT_spill->Write();

  for (int itdc = 0; itdc < 4; ++itdc) {
    hDeltaT_tdc[itdc]->Sumw2();
    if (nframes_used > 0)
      hDeltaT_tdc[itdc]->Scale(1. / nframes_used);
    hDeltaT_tdc[itdc]->Write();
  }

  fout->Close();
  fin->Close();
}

void
deltat(const std::string filename,
       int trigger_type,
       int trigger_device,
       int trigger_fifo,
       int trigger_column,
       int trigger_pixel,
       const std::string outfilename = "deltat.root")
{
  deltat(filename,
         field_selector_t{trigger_type, trigger_device, trigger_fifo, trigger_column, trigger_pixel},
         outfilename);
}

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
  deltat(filename,
         field_selector_t{target_type, target_device, target_fifo, target_column, target_pixel},
         field_selector_t{reference_type, reference_device, reference_fifo, reference_column, reference_pixel},
         outfilename);
}
