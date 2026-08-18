#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

struct source_t {
  int device;
  int fifo;
};

struct hit_t {
  int device;
  int fifo;
  int type;
  int counter;
  int column;
  int pixel;
  int tdc;
  int rollover;
  int coarse;
  int fine;
  double time;
  double x;
  double y;
};

class trigger_reader_t {
public:
  trigger_reader_t() = default;

  bool open(const std::string &filename);

  bool next_spill();
  bool next_frame();

  int spill_id() const;
  int frame_index() const;
  int nframes() const;
  int nsources() const;

  // Source and hit vectors describe the current spill/frame and remain valid
  // until open(), next_spill(), or next_frame() is called again.
  const std::vector<source_t> &sources() const;
  const std::vector<hit_t> &trigger_hits() const;
  const std::vector<hit_t> &timing_hits() const;
  const std::vector<hit_t> &cherenkov_hits() const;

private:
  struct category_t {
    std::unique_ptr<TTreeReaderValue<int>> nhits;
    std::unique_ptr<TTreeReaderArray<int>> frame_start;
    std::unique_ptr<TTreeReaderArray<int>> frame_nhits;
    std::unique_ptr<TTreeReaderArray<int>> device;
    std::unique_ptr<TTreeReaderArray<int>> fifo;
    std::unique_ptr<TTreeReaderArray<int>> type;
    std::unique_ptr<TTreeReaderArray<int>> counter;
    std::unique_ptr<TTreeReaderArray<int>> column;
    std::unique_ptr<TTreeReaderArray<int>> pixel;
    std::unique_ptr<TTreeReaderArray<int>> tdc;
    std::unique_ptr<TTreeReaderArray<int>> rollover;
    std::unique_ptr<TTreeReaderArray<int>> coarse;
    std::unique_ptr<TTreeReaderArray<int>> fine;
    std::unique_ptr<TTreeReaderArray<double>> time;
    std::unique_ptr<TTreeReaderArray<double>> x;
    std::unique_ptr<TTreeReaderArray<double>> y;

    void reset();
    bool bind(TTree *tree, TTreeReader &reader, const std::string &prefix);
    bool fill(std::vector<hit_t> &out, int iframe, int nframes,
              const std::string &prefix) const;
  };

  void reset();
  static bool has_branch(TTree *tree, const std::string &name);

  std::unique_ptr<TFile> file_;
  TTree *tree_ = nullptr;
  TTree *meta_tree_ = nullptr;
  std::unique_ptr<TTreeReader> reader_;
  std::unique_ptr<TTreeReader> meta_reader_;
  std::unique_ptr<TTreeReaderValue<int>> id_branch_;
  std::unique_ptr<TTreeReaderValue<int>> nframes_branch_;
  std::unique_ptr<TTreeReaderValue<int>> meta_counter_branch_;
  std::unique_ptr<TTreeReaderValue<int>> meta_nsources_branch_;
  std::unique_ptr<TTreeReaderArray<int>> meta_source_device_branch_;
  std::unique_ptr<TTreeReaderArray<int>> meta_source_fifo_branch_;

  category_t trigger_;
  category_t timing_;
  category_t cherenkov_;

  Long64_t nspills_ = 0;
  Long64_t spill_entry_ = -1;
  int spill_id_ = 0;
  int nframes_ = 0;
  int frame_index_ = -1;

  std::vector<source_t> sources_;
  std::vector<hit_t> trigger_hits_;
  std::vector<hit_t> timing_hits_;
  std::vector<hit_t> cherenkov_hits_;
};

inline bool
trigger_reader_t::has_branch(TTree *tree, const std::string &name)
{
  if (tree && tree->GetBranch(name.c_str()))
    return true;

  std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
  return false;
}

inline void
trigger_reader_t::category_t::reset()
{
  nhits.reset();
  frame_start.reset();
  frame_nhits.reset();
  device.reset();
  fifo.reset();
  type.reset();
  counter.reset();
  column.reset();
  pixel.reset();
  tdc.reset();
  rollover.reset();
  coarse.reset();
  fine.reset();
  time.reset();
  x.reset();
  y.reset();
}

inline bool
trigger_reader_t::category_t::bind(TTree *tree, TTreeReader &reader,
                                   const std::string &prefix)
{
  auto nname = std::string("n") + prefix + "hits";

  std::vector<std::string> names = {
    nname,
    prefix + "_frame_start",
    prefix + "_frame_nhits",
    prefix + "_device",
    prefix + "_fifo",
    prefix + "_type",
    prefix + "_counter",
    prefix + "_column",
    prefix + "_pixel",
    prefix + "_tdc",
    prefix + "_rollover",
    prefix + "_coarse",
    prefix + "_fine",
    prefix + "_time",
    prefix + "_x",
    prefix + "_y"
  };

  for (auto &name : names) {
    if (!trigger_reader_t::has_branch(tree, name))
      return false;
  }

  nhits.reset(new TTreeReaderValue<int>(reader, nname.c_str()));
  frame_start.reset(new TTreeReaderArray<int>(reader, (prefix + "_frame_start").c_str()));
  frame_nhits.reset(new TTreeReaderArray<int>(reader, (prefix + "_frame_nhits").c_str()));
  device.reset(new TTreeReaderArray<int>(reader, (prefix + "_device").c_str()));
  fifo.reset(new TTreeReaderArray<int>(reader, (prefix + "_fifo").c_str()));
  type.reset(new TTreeReaderArray<int>(reader, (prefix + "_type").c_str()));
  counter.reset(new TTreeReaderArray<int>(reader, (prefix + "_counter").c_str()));
  column.reset(new TTreeReaderArray<int>(reader, (prefix + "_column").c_str()));
  pixel.reset(new TTreeReaderArray<int>(reader, (prefix + "_pixel").c_str()));
  tdc.reset(new TTreeReaderArray<int>(reader, (prefix + "_tdc").c_str()));
  rollover.reset(new TTreeReaderArray<int>(reader, (prefix + "_rollover").c_str()));
  coarse.reset(new TTreeReaderArray<int>(reader, (prefix + "_coarse").c_str()));
  fine.reset(new TTreeReaderArray<int>(reader, (prefix + "_fine").c_str()));
  time.reset(new TTreeReaderArray<double>(reader, (prefix + "_time").c_str()));
  x.reset(new TTreeReaderArray<double>(reader, (prefix + "_x").c_str()));
  y.reset(new TTreeReaderArray<double>(reader, (prefix + "_y").c_str()));

  return true;
}

inline bool
trigger_reader_t::category_t::fill(std::vector<hit_t> &out, int iframe,
                                   int nframes, const std::string &prefix) const
{
  out.clear();

  if (iframe < 0 || iframe >= nframes) {
    std::cerr << "ERROR: invalid frame index " << iframe
              << " for category '" << prefix << "'" << std::endl;
    return false;
  }

  int total = **nhits;
  int first = (*frame_start)[iframe];
  int n = (*frame_nhits)[iframe];

  if (total < 0 || first < 0 || n < 0 || first + n > total) {
    std::cerr << "ERROR: invalid '" << prefix << "' frame indexing"
              << " frame=" << iframe
              << " first=" << first
              << " nhits=" << n
              << " total=" << total
              << std::endl;
    return false;
  }

  out.reserve(n);
  for (int i = first; i < first + n; ++i) {
    hit_t hit;
    hit.device = (*device)[i];
    hit.fifo = (*fifo)[i];
    hit.type = (*type)[i];
    hit.counter = (*counter)[i];
    hit.column = (*column)[i];
    hit.pixel = (*pixel)[i];
    hit.tdc = (*tdc)[i];
    hit.rollover = (*rollover)[i];
    hit.coarse = (*coarse)[i];
    hit.fine = (*fine)[i];
    hit.time = (*time)[i];
    hit.x = (*x)[i];
    hit.y = (*y)[i];
    out.push_back(hit);
  }

  return true;
}

inline void
trigger_reader_t::reset()
{
  sources_.clear();
  trigger_hits_.clear();
  timing_hits_.clear();
  cherenkov_hits_.clear();

  trigger_.reset();
  timing_.reset();
  cherenkov_.reset();

  id_branch_.reset();
  nframes_branch_.reset();
  meta_counter_branch_.reset();
  meta_nsources_branch_.reset();
  meta_source_device_branch_.reset();
  meta_source_fifo_branch_.reset();
  meta_reader_.reset();
  reader_.reset();
  meta_tree_ = nullptr;
  tree_ = nullptr;
  file_.reset();

  nspills_ = 0;
  spill_entry_ = -1;
  spill_id_ = 0;
  nframes_ = 0;
  frame_index_ = -1;
}

inline bool
trigger_reader_t::open(const std::string &filename)
{
  reset();

  file_.reset(TFile::Open(filename.c_str(), "READ"));
  if (!file_ || file_->IsZombie()) {
    std::cerr << "ERROR: could not open input file '" << filename << "'" << std::endl;
    reset();
    return false;
  }

  tree_ = (TTree *)file_->Get("frames");
  if (!tree_) {
    std::cerr << "ERROR: could not find 'frames' tree in '" << filename << "'" << std::endl;
    reset();
    return false;
  }

  if (!has_branch(tree_, "id") || !has_branch(tree_, "nframes")) {
    reset();
    return false;
  }

  reader_.reset(new TTreeReader(tree_));
  id_branch_.reset(new TTreeReaderValue<int>(*reader_, "id"));
  nframes_branch_.reset(new TTreeReaderValue<int>(*reader_, "nframes"));

  if (!trigger_.bind(tree_, *reader_, "trigger") ||
      !timing_.bind(tree_, *reader_, "timing") ||
      !cherenkov_.bind(tree_, *reader_, "cherenkov")) {
    reset();
    return false;
  }

  meta_tree_ = (TTree *)file_->Get("spill_participation");
  if (meta_tree_) {
    for (auto name : {"counter", "nsources", "source_device", "source_fifo"}) {
      if (!has_branch(meta_tree_, name)) {
        reset();
        return false;
      }
    }
    if (meta_tree_->GetEntries() != tree_->GetEntries()) {
      std::cerr << "ERROR: frames/spill_participation entry-count mismatch" << std::endl;
      reset();
      return false;
    }
    meta_reader_.reset(new TTreeReader(meta_tree_));
    meta_counter_branch_.reset(new TTreeReaderValue<int>(*meta_reader_, "counter"));
    meta_nsources_branch_.reset(new TTreeReaderValue<int>(*meta_reader_, "nsources"));
    meta_source_device_branch_.reset(new TTreeReaderArray<int>(*meta_reader_, "source_device"));
    meta_source_fifo_branch_.reset(new TTreeReaderArray<int>(*meta_reader_, "source_fifo"));
  }

  nspills_ = tree_->GetEntries();
  spill_entry_ = -1;
  frame_index_ = -1;
  return true;
}

inline bool
trigger_reader_t::next_spill()
{
  sources_.clear();
  trigger_hits_.clear();
  timing_hits_.clear();
  cherenkov_hits_.clear();
  frame_index_ = -1;

  if (!reader_)
    return false;

  if (spill_entry_ + 1 >= nspills_)
    return false;

  if (!reader_->Next()) {
    std::cerr << "ERROR: failed to read spill entry " << (spill_entry_ + 1) << std::endl;
    return false;
  }

  ++spill_entry_;
  spill_id_ = **id_branch_;
  nframes_ = **nframes_branch_;

  if (meta_reader_) {
    if (!meta_reader_->Next()) {
      std::cerr << "ERROR: failed to read spill_participation entry " << spill_entry_ << std::endl;
      return false;
    }
    if (**meta_counter_branch_ != spill_id_) {
      std::cerr << "ERROR: spill_participation counter mismatch"
                << " frame_id=" << spill_id_
                << " meta_counter=" << **meta_counter_branch_ << std::endl;
      return false;
    }
    int nsrc = **meta_nsources_branch_;
    if (nsrc < 0 || meta_source_device_branch_->GetSize() < nsrc ||
        meta_source_fifo_branch_->GetSize() < nsrc) {
      std::cerr << "ERROR: invalid spill_participation source arrays" << std::endl;
      return false;
    }
    sources_.reserve(nsrc);
    for (int i = 0; i < nsrc; ++i) {
      source_t source;
      source.device = (*meta_source_device_branch_)[i];
      source.fifo = (*meta_source_fifo_branch_)[i];
      sources_.push_back(source);
    }
  }

  if (nframes_ < 0) {
    std::cerr << "ERROR: negative nframes in spill entry " << spill_entry_ << std::endl;
    return false;
  }

  return true;
}

inline bool
trigger_reader_t::next_frame()
{
  if (!reader_ || spill_entry_ < 0)
    return false;

  if (frame_index_ + 1 >= nframes_)
    return false;

  ++frame_index_;

  if (!trigger_.fill(trigger_hits_, frame_index_, nframes_, "trigger"))
    return false;
  if (!timing_.fill(timing_hits_, frame_index_, nframes_, "timing"))
    return false;
  if (!cherenkov_.fill(cherenkov_hits_, frame_index_, nframes_, "cherenkov"))
    return false;

  return true;
}

inline int
trigger_reader_t::spill_id() const
{
  return spill_id_;
}

inline int
trigger_reader_t::frame_index() const
{
  return frame_index_;
}

inline int
trigger_reader_t::nframes() const
{
  return nframes_;
}

inline int
trigger_reader_t::nsources() const
{
  return (int)sources_.size();
}

inline const std::vector<source_t> &
trigger_reader_t::sources() const
{
  return sources_;
}

inline const std::vector<hit_t> &
trigger_reader_t::trigger_hits() const
{
  return trigger_hits_;
}

inline const std::vector<hit_t> &
trigger_reader_t::timing_hits() const
{
  return timing_hits_;
}

inline const std::vector<hit_t> &
trigger_reader_t::cherenkov_hits() const
{
  return cherenkov_hits_;
}
