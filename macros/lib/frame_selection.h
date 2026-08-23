#pragma once

#include "trigger_reader.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

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

class target_t {
public:
  virtual ~target_t() = default;

  // Called once after next_frame(). A target may reject the complete frame.
  virtual bool process(const trigger_reader_t &reader) = 0;
  virtual bool matches(const hit_t &hit) const = 0;
};

class reference_t {
public:
  virtual ~reference_t() = default;

  virtual bool process(const trigger_reader_t &reader) = 0;
  virtual double reference_time() const = 0;
};

class selection_t {
public:
  virtual ~selection_t() = default;

  virtual bool is_selected(const trigger_reader_t &reader) const = 0;
};

namespace frame_selection_detail {

inline bool
match_field(int value, int requested)
{
  return requested < 0 || value == requested;
}

inline int
channel_number(const hit_t &hit)
{
  const int fifo = hit.type == 9 ? 32 : hit.fifo;
  const int chip = fifo / 4;
  const int local = hit.pixel + 4 * hit.column + 32 * chip;
  return local + 256 * (hit.device - 192);
}

inline bool
matches(const hit_t &hit, const field_selector_t &selector)
{
  return match_field(hit.type, selector.type) &&
         match_field(hit.device, selector.device) &&
         match_field(hit.fifo, selector.fifo) &&
         match_field(hit.column, selector.column) &&
         match_field(hit.pixel, selector.pixel);
}

inline bool
matches(const hit_t &hit, const channel_selector_t &selector)
{
  if (selector.type == 9)
    return hit.type == 9 && match_field(hit.device, selector.channel);

  return hit.type == selector.type &&
         match_field(channel_number(hit), selector.channel);
}

template <typename selector_t>
inline bool
find_reference_time(const trigger_reader_t &reader,
                    const selector_t &selector,
                    double &value)
{
  auto scan = [&](const std::vector<hit_t> &hits) {
    for (const auto &hit : hits) {
      if (std::isfinite(hit.time) && matches(hit, selector)) {
        value = hit.time;
        return true;
      }
    }
    return false;
  };

  return scan(reader.trigger_hits()) ||
         scan(reader.timing_hits()) ||
         scan(reader.cherenkov_hits());
}

inline double
ring_residual(const ring_t &ring, const hit_t &hit)
{
  const double dx = hit.x - ring.x0;
  const double dy = hit.y - ring.y0;
  const double c = std::cos(ring.phi);
  const double s = std::sin(ring.phi);
  const double xp = c * dx + s * dy;
  const double yp = -s * dx + c * dy;
  const double minor = ring.radius * std::sqrt(std::max(
    1.e-12, 1. - ring.eccentricity * ring.eccentricity));
  if (ring.radius <= 0. || minor <= 0.)
    return std::numeric_limits<double>::infinity();
  return std::abs(std::hypot(xp / ring.radius, yp / minor) - 1.) * ring.radius;
}

inline bool
ring_contains(const ring_t &ring, const hit_t &hit, double tolerance)
{
  return std::isfinite(hit.x) && std::isfinite(hit.y) &&
         ring_residual(ring, hit) <= tolerance;
}

inline bool
passes_trigger(const trigger_reader_t &reader, bool enabled, int device)
{
  if (!enabled)
    return true;
  for (const auto &hit : reader.trigger_hits()) {
    if (device < 0 || hit.device == device)
      return true;
  }
  return false;
}

} // namespace frame_selection_detail

struct field_target_t : target_t {
  field_selector_t selector;

  field_target_t(int device, int fifo = -1, int column = -1, int pixel = -1)
    : selector{1, device, fifo, column, pixel}
  {
  }

  bool process(const trigger_reader_t &) override { return true; }
  bool matches(const hit_t &hit) const override
  {
    return frame_selection_detail::matches(hit, selector);
  }
};

struct channel_target_t : target_t {
  channel_selector_t selector;

  explicit channel_target_t(int channel)
    : selector{1, channel}
  {
  }

  bool process(const trigger_reader_t &) override { return true; }
  bool matches(const hit_t &hit) const override
  {
    return frame_selection_detail::matches(hit, selector);
  }
};

struct fifo_target_t : target_t {
  field_selector_t selector;

  fifo_target_t(int device, int fifo)
    : selector{1, device, fifo, -1, -1}
  {
  }

  bool process(const trigger_reader_t &) override { return true; }
  bool matches(const hit_t &hit) const override
  {
    return frame_selection_detail::matches(hit, selector);
  }
};

struct trigger_target_t : target_t {
  field_selector_t selector;

  explicit trigger_target_t(int device = -1)
    : selector{9, device, -1, -1, -1}
  {
  }

  bool process(const trigger_reader_t &) override { return true; }
  bool matches(const hit_t &hit) const override
  {
    return frame_selection_detail::matches(hit, selector);
  }
};

struct field_reference_t : reference_t {
  field_selector_t selector;
  double value = std::numeric_limits<double>::quiet_NaN();

  explicit field_reference_t(field_selector_t selector_)
    : selector(selector_)
  {
  }

  bool process(const trigger_reader_t &reader) override
  {
    value = std::numeric_limits<double>::quiet_NaN();
    return frame_selection_detail::find_reference_time(reader, selector, value);
  }

  double reference_time() const override { return value; }
};

struct channel_reference_t : reference_t {
  channel_selector_t selector;
  double value = std::numeric_limits<double>::quiet_NaN();

  explicit channel_reference_t(channel_selector_t selector_)
    : selector(selector_)
  {
  }

  bool process(const trigger_reader_t &reader) override
  {
    value = std::numeric_limits<double>::quiet_NaN();
    return frame_selection_detail::find_reference_time(reader, selector, value);
  }

  double reference_time() const override { return value; }
};

struct timing_reference_t : reference_t {
  std::string branch;
  double value = std::numeric_limits<double>::quiet_NaN();

  explicit timing_reference_t(const std::string &branch_ = "T")
    : branch(branch_)
  {
  }

  bool process(const trigger_reader_t &reader) override
  {
    value = std::numeric_limits<double>::quiet_NaN();
    if (!reader.has_timing()) {
      std::cerr << "ERROR: input file does not contain timing estimator branches"
                << std::endl;
      return false;
    }
    if (reader.timing_valid() == 0)
      return false;

    if (branch == "T0")
      value = reader.T0();
    else if (branch == "T1")
      value = reader.T1();
    else if (branch == "T")
      value = reader.T();
    else {
      std::cerr << "ERROR: unsupported timing reference branch: " << branch
                << " (valid choices are T, T0, T1)" << std::endl;
      return false;
    }
    return std::isfinite(value);
  }

  double reference_time() const override { return value; }
};

struct trigger_reference_t : reference_t {
  int device = -1;
  double value = std::numeric_limits<double>::quiet_NaN();

  explicit trigger_reference_t(int device_ = -1)
    : device(device_)
  {
  }

  bool process(const trigger_reader_t &reader) override
  {
    value = std::numeric_limits<double>::quiet_NaN();
    return frame_selection_detail::find_reference_time(
        reader, field_selector_t{9, device, -1, -1, -1}, value);
  }

  double reference_time() const override { return value; }
};

struct trigger_selection_t : selection_t {
  bool enabled = false;
  int device = -1;

  trigger_selection_t() = default;

  explicit trigger_selection_t(int device_)
    : enabled(true), device(device_)
  {
  }

  bool is_selected(const trigger_reader_t &reader) const override
  {
    return frame_selection_detail::passes_trigger(reader, enabled, device);
  }
};

struct ring_selection_t : selection_t {
  double min_x0 = -std::numeric_limits<double>::infinity();
  double max_x0 =  std::numeric_limits<double>::infinity();
  double min_y0 = -std::numeric_limits<double>::infinity();
  double max_y0 =  std::numeric_limits<double>::infinity();
  double min_r = 0.;
  double max_r = std::numeric_limits<double>::infinity();
  int min_n_rings = 1;
  int max_n_rings = std::numeric_limits<int>::max();

  ring_selection_t(double min_x0_ = -std::numeric_limits<double>::infinity(),
                   double max_x0_ =  std::numeric_limits<double>::infinity(),
                   double min_y0_ = -std::numeric_limits<double>::infinity(),
                   double max_y0_ =  std::numeric_limits<double>::infinity(),
                   double min_r_ = 0.,
                   double max_r_ = std::numeric_limits<double>::infinity(),
                   int min_n_rings_ = 1,
                   int max_n_rings_ = std::numeric_limits<int>::max())
    : min_x0(min_x0_), max_x0(max_x0_), min_y0(min_y0_), max_y0(max_y0_),
      min_r(min_r_), max_r(max_r_), min_n_rings(min_n_rings_),
      max_n_rings(max_n_rings_)
  {
  }

  bool is_selected(const trigger_reader_t &reader) const override
  {
    const int n = static_cast<int>(reader.rings().size());
    if (n < min_n_rings || n > max_n_rings)
      return false;
    if (n == 0)
      return true;

    for (const auto &ring : reader.rings()) {
      if (ring.x0 >= min_x0 && ring.x0 <= max_x0 &&
          ring.y0 >= min_y0 && ring.y0 <= max_y0 &&
          ring.radius >= min_r && ring.radius <= max_r)
        return true;
    }
    return false;
  }
};

struct ring_target_t : target_t {
  int min_inliers = 0;
  int max_inliers = std::numeric_limits<int>::max();
  double tolerance = 3.5;
  double min_radius = 1.;
  double max_radius = 200.;

  std::vector<ring_t> accepted_rings;

  ring_target_t(int min_inliers_ = 0,
                int max_inliers_ = std::numeric_limits<int>::max(),
                double tolerance_ = 3.5,
                double min_radius_ = 1.,
                double max_radius_ = 200.)
    : min_inliers(min_inliers_), max_inliers(max_inliers_),
      tolerance(tolerance_), min_radius(min_radius_), max_radius(max_radius_)
  {
  }

  bool process(const trigger_reader_t &reader) override
  {
    accepted_rings.clear();
    for (const auto &ring : reader.rings()) {
      if (ring.ninliers < min_inliers || ring.ninliers > max_inliers)
        continue;
      if (ring.radius < min_radius || ring.radius > max_radius)
        continue;
      accepted_rings.push_back(ring);
    }
    return !accepted_rings.empty();
  }

  bool matches(const hit_t &hit) const override
  {
    for (const auto &ring : accepted_rings) {
      if (frame_selection_detail::ring_contains(ring, hit, tolerance))
        return true;
    }
    return false;
  }
};

using target_ptr_t = std::shared_ptr<target_t>;
using reference_ptr_t = std::shared_ptr<reference_t>;
using selection_ptr_t = std::shared_ptr<selection_t>;
