#pragma once

#include <TTree.h>

struct data_t {
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

  double time = 0.;
  double x = 0.;
  double y = 0.;
  bool has_time = false;
  bool has_x = false;
  bool has_y = false;

  static const int rollover_to_clock = 32768;

  bool is_alcor_hit()   const { return type == 1; }
  bool is_trigger_tag() const { return type == 9; }
  bool is_start_spill() const { return type == 7; }
  bool is_end_spill()   const { return type == 15; }
  bool is_data_word()   const { return is_alcor_hit() || is_trigger_tag(); }

  int channel() const { return pixel + 4 * column; }

  void link_to_tree(TTree *t)
  {
    if (!t) return;
    has_time = false;
    has_x = false;
    has_y = false;
    t->SetBranchAddress("device", &device);
    t->SetBranchAddress("fifo", &fifo);
    t->SetBranchAddress("type", &type);
    t->SetBranchAddress("counter", &counter);
    t->SetBranchAddress("column", &column);
    t->SetBranchAddress("pixel", &pixel);
    t->SetBranchAddress("tdc", &tdc);
    t->SetBranchAddress("rollover", &rollover);
    t->SetBranchAddress("coarse", &coarse);
    t->SetBranchAddress("fine", &fine);
    if (t->GetBranch("time")) {
      t->SetBranchAddress("time", &time);
      has_time = true;
    }
    if (t->GetBranch("x")) {
      t->SetBranchAddress("x", &x);
      has_x = true;
    }
    if (t->GetBranch("y")) {
      t->SetBranchAddress("y", &y);
      has_y = true;
    }
  }

  void set_nominal_time()
  {
    time = coarse + rollover_to_clock * rollover;
    if (is_alcor_hit())
      time -= 0.0157 * fine;
  }

  void set_time()
  {
    if (!has_time)
      set_nominal_time();
  }
};
