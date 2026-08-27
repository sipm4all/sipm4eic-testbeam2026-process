#pragma once

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

struct source_t {
  int device = 0;
  int fifo = 0;
};

struct hit_t {
  int device = 0;
  int fifo = 0;
  int type = 0;
  int counter = 0;
  int column = 0;
  int pixel = 0;
  int tdc = 0;
  int rollover = 0;
  int coarse = 0;
  int fine = 0;
  double time = 0.;
  double x = std::numeric_limits<double>::quiet_NaN();
  double y = std::numeric_limits<double>::quiet_NaN();
};

struct ring_t {
  double x0 = 0.;
  double y0 = 0.;
  double radius = 0.;
  double eccentricity = 0.;
  double phi = 0.;
  double time = 0.;
  int ninliers = 0;
};

class trigger_reader_t {
public:
  trigger_reader_t() = default;

  bool open(const std::string &filename,
            const std::string &ring_name = "ring");

  bool next_spill();
  bool next_frame();

  int spill_id() const;
  int frame_index() const;
  int nframes() const;
  int nsources() const;
  bool has_timing() const;
  bool has_rings() const;
  int timing_valid() const;
  double T0() const;
  double sigma0() const;
  double T1() const;
  double sigma1() const;
  double T() const;
  double sigmaT() const;

  // These vectors describe the current frame and remain valid until the next
  // state-changing call.
  const std::vector<source_t> &sources() const;
  const std::vector<hit_t> &trigger_hits() const;
  const std::vector<hit_t> &timing_hits() const;
  const std::vector<hit_t> &cherenkov_hits() const;
  const std::vector<ring_t> &rings() const;

private:
  static constexpr int maxhits_ = 65535;
  static constexpr int maxrings_ = 8;
  static constexpr int maxsources_ = 4096;

  struct category_t {
    bool trigger = false;
    UShort_t nhits = 0;
    std::vector<UChar_t> device;
    std::vector<UChar_t> fifo;
    std::vector<UChar_t> column;
    std::vector<UChar_t> pixel;
    std::vector<UChar_t> tdc;
    std::vector<UShort_t> counter;
    std::vector<UShort_t> rollover;
    std::vector<UShort_t> coarse;
    std::vector<UShort_t> fine;
    std::vector<Float_t> time;
    std::vector<Float_t> x;
    std::vector<Float_t> y;

    void reset();
    bool bind(TTree *tree, bool is_trigger);
    bool load_hits(std::vector<hit_t> &out) const;
  };

  static bool has_branch(TTree *tree, const char *name);
  static bool bind(TTree *tree, const char *name, void *address);
  bool load_frame(Long64_t entry);
  bool load_spill_sources();
  void reset();

  std::unique_ptr<TFile> file_;
  TTree *frames_tree_ = nullptr;
  TTree *trigger_tree_ = nullptr;
  TTree *timing_tree_ = nullptr;
  TTree *cherenkov_tree_ = nullptr;
  TTree *ring_tree_ = nullptr;
  TTree *meta_tree_ = nullptr;
  std::string ring_name_ = "ring";

  UInt_t frame_spill_ = 0;
  Double_t frame_time_ = 0.;
  bool has_frame_time_ = false;

  category_t trigger_;
  category_t timing_;
  category_t cherenkov_;

  bool has_timing_ = false;
  bool timing_valid_ = false;
  Float_t T0_ = 0.;
  Float_t sigma0_ = 0.;
  Float_t T1_ = 0.;
  Float_t sigma1_ = 0.;
  Float_t T_ = 0.;
  Float_t sigmaT_ = 0.;

  bool has_rings_ = false;
  UChar_t nring_ = 0;
  Float_t ring_x0_[maxrings_];
  Float_t ring_y0_[maxrings_];
  Float_t ring_r_[maxrings_];
  Float_t ring_e_[maxrings_];
  Float_t ring_phi_[maxrings_];
  Float_t ring_time_[maxrings_];
  UShort_t ring_ninliers_[maxrings_];

  bool has_meta_ = false;
  Long64_t meta_entries_ = 0;
  Long64_t meta_entry_ = 0;
  int meta_spill_ = 0;
  int meta_counter_ = 0;
  int meta_nsources_ = 0;
  int meta_source_device_[maxsources_];
  int meta_source_fifo_[maxsources_];

  Long64_t entries_ = 0;
  Long64_t next_entry_ = 0;
  Long64_t spill_begin_ = -1;
  Long64_t spill_end_ = -1;
  Long64_t current_entry_ = -1;
  int spill_id_ = 0;
  int frame_index_ = -1;
  int nframes_ = 0;

  std::vector<source_t> sources_;
  std::vector<hit_t> trigger_hits_;
  std::vector<hit_t> timing_hits_;
  std::vector<hit_t> cherenkov_hits_;
  std::vector<ring_t> rings_;
};

inline bool
trigger_reader_t::has_branch(TTree *tree, const char *name)
{
  return tree && tree->GetBranch(name);
}

inline bool
trigger_reader_t::bind(TTree *tree, const char *name, void *address)
{
  if (!has_branch(tree, name)) {
    std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
    return false;
  }
  if (tree->SetBranchAddress(name, address) < 0) {
    std::cerr << "ERROR: could not bind branch '" << name << "'" << std::endl;
    return false;
  }
  return true;
}

inline void
trigger_reader_t::category_t::reset()
{
  nhits = 0;
  device.clear();
  fifo.clear();
  column.clear();
  pixel.clear();
  tdc.clear();
  counter.clear();
  rollover.clear();
  coarse.clear();
  fine.clear();
  time.clear();
  x.clear();
  y.clear();
}

inline bool
trigger_reader_t::category_t::bind(TTree *tree, bool is_trigger)
{
  trigger = is_trigger;
  const int n = trigger_reader_t::maxhits_;
  device.resize(n);
  rollover.resize(n);
  coarse.resize(n);
  time.resize(n);

  if (!trigger) {
    fifo.resize(n);
    column.resize(n);
    pixel.resize(n);
    tdc.resize(n);
    fine.resize(n);
    x.resize(n);
    y.resize(n);
  }
  if (!trigger)
    counter.clear();
  else
    counter.resize(n);

  if (!trigger_reader_t::bind(tree, "nhits", &nhits) ||
      !trigger_reader_t::bind(tree, "device", device.data()) ||
      !trigger_reader_t::bind(tree, "rollover", rollover.data()) ||
      !trigger_reader_t::bind(tree, "coarse", coarse.data()) ||
      !trigger_reader_t::bind(tree, "time", time.data()))
    return false;

  if (trigger) {
    return trigger_reader_t::bind(tree, "counter", counter.data());
  }

  return trigger_reader_t::bind(tree, "fifo", fifo.data()) &&
         trigger_reader_t::bind(tree, "column", column.data()) &&
         trigger_reader_t::bind(tree, "pixel", pixel.data()) &&
         trigger_reader_t::bind(tree, "tdc", tdc.data()) &&
         trigger_reader_t::bind(tree, "fine", fine.data()) &&
         trigger_reader_t::bind(tree, "x", x.data()) &&
         trigger_reader_t::bind(tree, "y", y.data());
}

inline bool
trigger_reader_t::category_t::load_hits(std::vector<hit_t> &out) const
{
  out.clear();
  if (nhits > trigger_reader_t::maxhits_) {
    std::cerr << "ERROR: category contains too many hits: " << nhits << std::endl;
    return false;
  }

  out.reserve(nhits);
  for (unsigned int i = 0; i < nhits; ++i) {
    hit_t hit;
    hit.device = device[i];
    hit.type = trigger ? 9 : 1;
    hit.counter = trigger ? counter[i] : 0;
    hit.rollover = rollover[i];
    hit.coarse = coarse[i];
    hit.time = time[i];
    if (trigger) {
      hit.fifo = 0;
      hit.column = 0;
      hit.pixel = 0;
      hit.tdc = 0;
    } else {
      hit.fifo = fifo[i];
      hit.column = column[i];
      hit.pixel = pixel[i];
      hit.tdc = tdc[i];
      hit.fine = fine[i];
      hit.x = x[i];
      hit.y = y[i];
    }
    out.push_back(hit);
  }
  return true;
}

inline void
trigger_reader_t::reset()
{
  trigger_hits_.clear();
  timing_hits_.clear();
  cherenkov_hits_.clear();
  sources_.clear();
  rings_.clear();

  file_.reset();
  frames_tree_ = nullptr;
  trigger_tree_ = nullptr;
  timing_tree_ = nullptr;
  cherenkov_tree_ = nullptr;
  ring_tree_ = nullptr;
  meta_tree_ = nullptr;
  trigger_.reset();
  timing_.reset();
  cherenkov_.reset();

  has_frame_time_ = false;
  has_timing_ = false;
  has_rings_ = false;
  has_meta_ = false;
  meta_entries_ = 0;
  meta_entry_ = 0;
  entries_ = 0;
  next_entry_ = 0;
  spill_begin_ = -1;
  spill_end_ = -1;
  current_entry_ = -1;
  spill_id_ = 0;
  frame_index_ = -1;
  nframes_ = 0;
}

inline bool
trigger_reader_t::open(const std::string &filename,
                       const std::string &ring_name)
{
  reset();
  if (ring_name.empty()) {
    std::cerr << "ERROR: ring tree name is empty" << std::endl;
    return false;
  }
  ring_name_ = ring_name;
  file_.reset(TFile::Open(filename.c_str(), "READ"));
  if (!file_ || file_->IsZombie()) {
    std::cerr << "ERROR: could not open input file '" << filename << "'" << std::endl;
    reset();
    return false;
  }

  frames_tree_ = (TTree *)file_->Get("frames");
  trigger_tree_ = (TTree *)file_->Get("trigger");
  timing_tree_ = (TTree *)file_->Get("timing");
  cherenkov_tree_ = (TTree *)file_->Get("cherenkov");
  if (!frames_tree_ || !trigger_tree_ || !timing_tree_ || !cherenkov_tree_) {
    std::cerr << "ERROR: input must contain frames, trigger, timing, and cherenkov trees"
              << std::endl;
    reset();
    return false;
  }

  if (!bind(frames_tree_, "spill", &frame_spill_)) {
    reset();
    return false;
  }
  has_frame_time_ = has_branch(frames_tree_, "time");
  if (has_frame_time_ && !bind(frames_tree_, "time", &frame_time_)) {
    reset();
    return false;
  }

  if (!trigger_.bind(trigger_tree_, true) ||
      !timing_.bind(timing_tree_, false) ||
      !cherenkov_.bind(cherenkov_tree_, false)) {
    reset();
    return false;
  }

  entries_ = frames_tree_->GetEntries();
  if (trigger_tree_->GetEntries() != entries_ ||
      timing_tree_->GetEntries() != entries_ ||
      cherenkov_tree_->GetEntries() != entries_) {
    std::cerr << "ERROR: frame-tree entry-count mismatch" << std::endl;
    reset();
    return false;
  }

  const bool timing_any = timing_tree_->GetBranch("timing_valid") ||
                          timing_tree_->GetBranch("T0") ||
                          timing_tree_->GetBranch("sigma0") ||
                          timing_tree_->GetBranch("T1") ||
                          timing_tree_->GetBranch("sigma1") ||
                          timing_tree_->GetBranch("T") ||
                          timing_tree_->GetBranch("sigmaT");
  const bool timing_all = timing_tree_->GetBranch("timing_valid") &&
                          timing_tree_->GetBranch("T0") &&
                          timing_tree_->GetBranch("sigma0") &&
                          timing_tree_->GetBranch("T1") &&
                          timing_tree_->GetBranch("sigma1") &&
                          timing_tree_->GetBranch("T") &&
                          timing_tree_->GetBranch("sigmaT");
  if (timing_any && !timing_all) {
    std::cerr << "ERROR: timing estimator branches are incomplete" << std::endl;
    reset();
    return false;
  }
  has_timing_ = timing_all;
  if (has_timing_) {
    if (!bind(timing_tree_, "timing_valid", &timing_valid_) ||
        !bind(timing_tree_, "T0", &T0_) ||
        !bind(timing_tree_, "sigma0", &sigma0_) ||
        !bind(timing_tree_, "T1", &T1_) ||
        !bind(timing_tree_, "sigma1", &sigma1_) ||
        !bind(timing_tree_, "T", &T_) ||
        !bind(timing_tree_, "sigmaT", &sigmaT_)) {
      reset();
      return false;
    }
  }

  ring_tree_ = (TTree *)file_->Get(ring_name_.c_str());
  if (ring_tree_) {
    if (ring_tree_->GetEntries() != entries_ ||
        !bind(ring_tree_, "nrings", &nring_) ||
        !bind(ring_tree_, "x0", ring_x0_) ||
        !bind(ring_tree_, "y0", ring_y0_) ||
        !bind(ring_tree_, "r", ring_r_) ||
        !bind(ring_tree_, "e", ring_e_) ||
        !bind(ring_tree_, "phi", ring_phi_) ||
        !bind(ring_tree_, "time", ring_time_) ||
        !bind(ring_tree_, "ninliers", ring_ninliers_)) {
      std::cerr << "ERROR: invalid ring tree" << std::endl;
      reset();
      return false;
    }
    has_rings_ = true;
  }

  meta_tree_ = (TTree *)file_->Get("spill_participation");
  if (meta_tree_) {
    if (!bind(meta_tree_, "spill", &meta_spill_) ||
        !bind(meta_tree_, "counter", &meta_counter_) ||
        !bind(meta_tree_, "nsources", &meta_nsources_) ||
        !bind(meta_tree_, "source_device", meta_source_device_) ||
        !bind(meta_tree_, "source_fifo", meta_source_fifo_)) {
      reset();
      return false;
    }
    meta_entries_ = meta_tree_->GetEntries();
    has_meta_ = true;
  }

  return true;
}

inline bool
trigger_reader_t::load_spill_sources()
{
  sources_.clear();
  if (!has_meta_)
    return true;

  // The metadata tree contains one entry for every input spill, whereas the
  // frames tree contains entries only for spills with accepted frames. Find
  // the metadata entry by raw spill counter so empty trigger spills do not
  // shift the metadata/frame association.
  bool found = false;
  while (meta_entry_ < meta_entries_) {
    const Long64_t entry = meta_entry_;
    if (meta_tree_->GetEntry(entry) <= 0) {
      std::cerr << "ERROR: failed to read spill_participation entry "
                << entry << std::endl;
      return false;
    }
    ++meta_entry_;

    if (meta_counter_ < spill_id_)
      continue;
    if (meta_counter_ > spill_id_) {
      std::cerr << "ERROR: missing spill_participation entry for spill "
                << spill_id_ << "; next metadata counter=" << meta_counter_
                << std::endl;
      return false;
    }
    found = true;
    break;
  }

  if (!found) {
    std::cerr << "ERROR: missing spill_participation entry for spill "
              << spill_id_ << std::endl;
    return false;
  }
  if (meta_nsources_ < 0 || meta_nsources_ > maxsources_) {
    std::cerr << "ERROR: invalid spill_participation nsources=" << meta_nsources_ << std::endl;
    return false;
  }
  sources_.reserve(meta_nsources_);
  for (int i = 0; i < meta_nsources_; ++i)
    sources_.push_back({meta_source_device_[i], meta_source_fifo_[i]});
  return true;
}

inline bool
trigger_reader_t::load_frame(Long64_t entry)
{
  if (frames_tree_->GetEntry(entry) <= 0 ||
      trigger_tree_->GetEntry(entry) <= 0 ||
      timing_tree_->GetEntry(entry) <= 0 ||
      cherenkov_tree_->GetEntry(entry) <= 0) {
    std::cerr << "ERROR: failed to read frame entry " << entry << std::endl;
    return false;
  }
  if (has_rings_ && ring_tree_->GetEntry(entry) <= 0) {
    std::cerr << "ERROR: failed to read ring entry " << entry << std::endl;
    return false;
  }
  if (!trigger_.load_hits(trigger_hits_) ||
      !timing_.load_hits(timing_hits_) ||
      !cherenkov_.load_hits(cherenkov_hits_))
    return false;

  rings_.clear();
  if (has_rings_) {
    if (nring_ > maxrings_) {
      std::cerr << "ERROR: invalid nrings=" << static_cast<int>(nring_) << std::endl;
      return false;
    }
    rings_.reserve(nring_);
    for (unsigned int i = 0; i < nring_; ++i)
      rings_.push_back({ring_x0_[i], ring_y0_[i], ring_r_[i], ring_e_[i],
                        ring_phi_[i], ring_time_[i], ring_ninliers_[i]});
  }
  current_entry_ = entry;
  return true;
}

inline bool
trigger_reader_t::next_spill()
{
  trigger_hits_.clear();
  timing_hits_.clear();
  cherenkov_hits_.clear();
  rings_.clear();
  frame_index_ = -1;

  if (!frames_tree_ || next_entry_ >= entries_)
    return false;

  if (frames_tree_->GetEntry(next_entry_) <= 0) {
    std::cerr << "ERROR: failed to read frame entry " << next_entry_ << std::endl;
    return false;
  }
  spill_begin_ = next_entry_;
  spill_id_ = static_cast<int>(frame_spill_);
  spill_end_ = spill_begin_ + 1;
  while (spill_end_ < entries_) {
    if (frames_tree_->GetEntry(spill_end_) <= 0) {
      std::cerr << "ERROR: failed to read frame entry " << spill_end_ << std::endl;
      return false;
    }
    if (frame_spill_ != static_cast<UInt_t>(spill_id_))
      break;
    ++spill_end_;
  }
  next_entry_ = spill_end_;
  current_entry_ = spill_begin_ - 1;
  nframes_ = static_cast<int>(spill_end_ - spill_begin_);
  return load_spill_sources();
}

inline bool
trigger_reader_t::next_frame()
{
  if (!frames_tree_ || spill_begin_ < 0 || current_entry_ + 1 >= spill_end_)
    return false;
  ++current_entry_;
  ++frame_index_;
  return load_frame(current_entry_);
}

inline int trigger_reader_t::spill_id() const { return spill_id_; }
inline int trigger_reader_t::frame_index() const { return frame_index_; }
inline int trigger_reader_t::nframes() const { return nframes_; }
inline int trigger_reader_t::nsources() const { return static_cast<int>(sources_.size()); }
inline bool trigger_reader_t::has_timing() const { return has_timing_; }
inline bool trigger_reader_t::has_rings() const { return has_rings_; }

inline int
trigger_reader_t::timing_valid() const
{
  return has_timing_ && frame_index_ >= 0 && timing_valid_ ? 1 : 0;
}

inline double trigger_reader_t::T0() const
{
  return has_timing_ && frame_index_ >= 0 ? T0_ : std::numeric_limits<double>::quiet_NaN();
}
inline double trigger_reader_t::sigma0() const
{
  return has_timing_ && frame_index_ >= 0 ? sigma0_ : std::numeric_limits<double>::quiet_NaN();
}
inline double trigger_reader_t::T1() const
{
  return has_timing_ && frame_index_ >= 0 ? T1_ : std::numeric_limits<double>::quiet_NaN();
}
inline double trigger_reader_t::sigma1() const
{
  return has_timing_ && frame_index_ >= 0 ? sigma1_ : std::numeric_limits<double>::quiet_NaN();
}
inline double trigger_reader_t::T() const
{
  return has_timing_ && frame_index_ >= 0 ? T_ : std::numeric_limits<double>::quiet_NaN();
}
inline double trigger_reader_t::sigmaT() const
{
  return has_timing_ && frame_index_ >= 0 ? sigmaT_ : std::numeric_limits<double>::quiet_NaN();
}

inline const std::vector<source_t> &trigger_reader_t::sources() const { return sources_; }
inline const std::vector<hit_t> &trigger_reader_t::trigger_hits() const { return trigger_hits_; }
inline const std::vector<hit_t> &trigger_reader_t::timing_hits() const { return timing_hits_; }
inline const std::vector<hit_t> &trigger_reader_t::cherenkov_hits() const { return cherenkov_hits_; }
inline const std::vector<ring_t> &trigger_reader_t::rings() const { return rings_; }
