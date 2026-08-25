#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "ring-finder-hough-cuda.h"

#include <boost/program_options.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int maxrings = 8;
constexpr double pi = 3.14159265358979323846;
constexpr double gaussian_cut = 4.;

struct point_t {
  double x;
  double y;
  double time;
};

struct circle_t {
  double x = 0.;
  double y = 0.;
  double radius = 0.;
  double ring_time = 0.;
  int inliers = 0;
  double residual_sum = std::numeric_limits<double>::infinity();
  double score = 0.;
  bool valid = false;
};

double
radial_residual(const circle_t &circle, const point_t &point)
{
  return std::hypot(point.x - circle.x, point.y - circle.y) - circle.radius;
}

bool
circle_from_three(const point_t &a, const point_t &b, const point_t &c,
                  circle_t &circle)
{
  const double determinant = 2. * (a.x * (b.y - c.y) +
                                   b.x * (c.y - a.y) +
                                   c.x * (a.y - b.y));
  if (std::abs(determinant) < 1.e-9)
    return false;

  const double aa = a.x * a.x + a.y * a.y;
  const double bb = b.x * b.x + b.y * b.y;
  const double cc = c.x * c.x + c.y * c.y;
  circle.x = (aa * (b.y - c.y) + bb * (c.y - a.y) +
              cc * (a.y - b.y)) / determinant;
  circle.y = (aa * (c.x - b.x) + bb * (a.x - c.x) +
              cc * (b.x - a.x)) / determinant;
  circle.radius = std::hypot(a.x - circle.x, a.y - circle.y);
  return std::isfinite(circle.x) && std::isfinite(circle.y) &&
         std::isfinite(circle.radius);
}

double
mean_time(const std::vector<point_t> &points,
          const std::vector<int> &indices)
{
  if (indices.empty())
    return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.;
  for (int index : indices)
    sum += points[index].time;
  return sum / indices.size();
}

bool
fit_circle(const std::vector<point_t> &points,
            const std::vector<int> &indices,
            circle_t &circle)
{
  if (indices.size() < 3)
    return false;

  double normal[3][3] = {};
  double rhs[3] = {};
  for (int index : indices) {
    const double x = points[index].x;
    const double y = points[index].y;
    const double row[3] = {x, y, 1.};
    const double value = x * x + y * y;
    for (int i = 0; i < 3; ++i) {
      rhs[i] += row[i] * value;
      for (int j = 0; j < 3; ++j)
        normal[i][j] += row[i] * row[j];
    }
  }

  for (int column = 0; column < 3; ++column) {
    int pivot = column;
    for (int row = column + 1; row < 3; ++row) {
      if (std::abs(normal[row][column]) >
          std::abs(normal[pivot][column]))
        pivot = row;
    }
    if (std::abs(normal[pivot][column]) < 1.e-12)
      return false;
    if (pivot != column) {
      for (int j = column; j < 3; ++j)
        std::swap(normal[column][j], normal[pivot][j]);
      std::swap(rhs[column], rhs[pivot]);
    }
    for (int row = column + 1; row < 3; ++row) {
      const double factor = normal[row][column] / normal[column][column];
      for (int j = column; j < 3; ++j)
        normal[row][j] -= factor * normal[column][j];
      rhs[row] -= factor * rhs[column];
    }
  }

  double solution[3] = {};
  for (int row = 2; row >= 0; --row) {
    double value = rhs[row];
    for (int column = row + 1; column < 3; ++column)
      value -= normal[row][column] * solution[column];
    solution[row] = value / normal[row][row];
  }

  circle.x = 0.5 * solution[0];
  circle.y = 0.5 * solution[1];
  const double radius_squared = solution[2] + circle.x * circle.x +
                                circle.y * circle.y;
  if (radius_squared <= 0.)
    return false;
  circle.radius = std::sqrt(radius_squared);
  return std::isfinite(circle.x) && std::isfinite(circle.y) &&
         std::isfinite(circle.radius);
}

std::vector<circle_t>
ransac_seeds(const std::vector<point_t> &points,
            const std::vector<int> &indices,
            int max_seeds, int iterations, int min_inliers,
            double min_x0, double max_x0,
            double min_y0, double max_y0,
            double min_radius, double max_radius,
            double ransac_tolerance, double time_window,
            std::mt19937 &generator)
{
  std::vector<circle_t> hypotheses;
  if (indices.size() < 3 || max_seeds <= 0 || iterations <= 0)
    return hypotheses;

  const double spatial_cut = ransac_tolerance;
  const double time_cut = time_window;
  for (int iteration = 0; iteration < iterations; ++iteration) {
    std::uniform_int_distribution<std::size_t> pick(0, indices.size() - 1);
    const std::size_t ia = pick(generator);
    const double seed_time = points[indices[ia]].time;
    std::vector<int> temporal_indices;
    for (int index : indices) {
      if (std::abs(points[index].time - seed_time) <= time_cut)
        temporal_indices.push_back(index);
    }
    if (temporal_indices.size() < 3)
      continue;

    std::uniform_int_distribution<std::size_t> temporal_pick(
        0, temporal_indices.size() - 1);
    const std::size_t ib = temporal_pick(generator);
    const std::size_t ic = temporal_pick(generator);
    if (temporal_indices[ib] == indices[ia] ||
        temporal_indices[ic] == indices[ia] || ib == ic)
      continue;

    circle_t hypothesis;
    if (!circle_from_three(points[indices[ia]],
                           points[temporal_indices[ib]],
                           points[temporal_indices[ic]], hypothesis) ||
        hypothesis.x < min_x0 || hypothesis.x > max_x0 ||
        hypothesis.y < min_y0 || hypothesis.y > max_y0 ||
        hypothesis.radius < min_radius || hypothesis.radius > max_radius)
      continue;

    std::vector<int> accepted;
    auto collect_inliers = [&](const circle_t &model, double ring_time) {
      accepted.clear();
      for (int index : indices) {
        const double spatial =
            std::abs(radial_residual(model, points[index]));
        const double temporal = std::abs(points[index].time - ring_time);
        if (spatial <= spatial_cut && temporal <= time_cut)
          accepted.push_back(index);
      }
    };

    hypothesis.ring_time = points[indices[ia]].time;
    collect_inliers(hypothesis, hypothesis.ring_time);
    hypothesis.ring_time = mean_time(points, accepted);
    collect_inliers(hypothesis, hypothesis.ring_time);
    if (!fit_circle(points, accepted, hypothesis))
      continue;
    for (int refit = 0; refit < 2; ++refit) {
      hypothesis.ring_time = mean_time(points, accepted);
      collect_inliers(hypothesis, hypothesis.ring_time);
      if (!fit_circle(points, accepted, hypothesis)) {
        accepted.clear();
        break;
      }
    }
    if (accepted.size() < 3 || !std::isfinite(hypothesis.ring_time) ||
        hypothesis.x < min_x0 || hypothesis.x > max_x0 ||
        hypothesis.y < min_y0 || hypothesis.y > max_y0 ||
        hypothesis.radius < min_radius || hypothesis.radius > max_radius)
      continue;

    hypothesis.ring_time = mean_time(points, accepted);
    collect_inliers(hypothesis, hypothesis.ring_time);
    hypothesis.inliers = static_cast<int>(accepted.size());
    hypothesis.residual_sum = 0.;
    for (int index : accepted)
      hypothesis.residual_sum +=
          std::abs(radial_residual(hypothesis, points[index]));
    hypotheses.push_back(hypothesis);
  }

  std::vector<circle_t> seeds;
  seeds.reserve(max_seeds);
  std::vector<int> remaining = indices;
  for (int iseed = 0;
       iseed < max_seeds && static_cast<int>(remaining.size()) >= min_inliers;
       ++iseed) {
    int best_hypothesis = -1;
    int best_inliers = 0;
    double best_residual_sum = std::numeric_limits<double>::infinity();
    std::vector<int> best_inlier_indices;
    for (int ihypothesis = 0;
         ihypothesis < static_cast<int>(hypotheses.size()); ++ihypothesis) {
      const circle_t &hypothesis = hypotheses[ihypothesis];
      std::vector<int> inlier_indices;
      double residual_sum = 0.;
      for (int index : remaining) {
        const double spatial =
            std::abs(radial_residual(hypothesis, points[index]));
        const double temporal =
            std::abs(points[index].time - hypothesis.ring_time);
        if (spatial <= spatial_cut && temporal <= time_cut) {
          inlier_indices.push_back(index);
          residual_sum += spatial;
        }
      }
      const int ninliers = static_cast<int>(inlier_indices.size());
      if (ninliers < min_inliers ||
          ninliers < best_inliers ||
          (ninliers == best_inliers && residual_sum >= best_residual_sum))
        continue;
      best_hypothesis = ihypothesis;
      best_inliers = ninliers;
      best_residual_sum = residual_sum;
      best_inlier_indices.swap(inlier_indices);
    }
    if (best_hypothesis < 0)
      break;

    circle_t seed = hypotheses[best_hypothesis];
    seed.inliers = best_inliers;
    seed.residual_sum = best_residual_sum;
    seeds.push_back(seed);

    std::vector<int> next_remaining;
    next_remaining.reserve(remaining.size() - best_inlier_indices.size());
    for (int index : remaining) {
      if (std::find(best_inlier_indices.begin(),
                    best_inlier_indices.end(), index) ==
          best_inlier_indices.end())
        next_remaining.push_back(index);
      }
    remaining.swap(next_remaining);
  }
  return seeds;
}

bool
has_branch(TTree *tree, const char *name)
{
  if (tree->GetBranch(name))
    return true;
  std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
  return false;
}

bool
grid_bins(double minimum, double maximum, double step, int &bins)
{
  if (!std::isfinite(minimum) || !std::isfinite(maximum) ||
      !std::isfinite(step) || step <= 0. || minimum > maximum)
    return false;
  const double count = std::ceil((maximum - minimum) / step) + 1.;
  if (count < 1. || count > 100000.)
    return false;
  bins = static_cast<int>(count);
  return true;
}

void
align_grid_bounds(double global_min, double global_max, double step,
                  double &local_min, double &local_max)
{
  const double lower = global_min +
                       std::floor((local_min - global_min) / step) * step;
  const double upper = global_min +
                       std::ceil((local_max - global_min) / step) * step;
  local_min = std::max(global_min, lower);
  local_max = std::min(global_max, upper);
}

std::vector<circle_t>
hough_candidates(const std::vector<point_t> &points,
               const std::vector<int> &indices,
               int max_candidates,
               double min_x0, double max_x0, double x0_step,
               double min_y0, double max_y0, double y0_step,
               double min_radius, double max_radius, double radius_step,
               double min_t, double max_t, double t_step,
               double spatial_resolution, double time_resolution)
{
  std::vector<circle_t> empty;
  int nx = 0, ny = 0, nr = 0, nt = 0;
  if (!grid_bins(min_x0, max_x0, x0_step, nx) ||
      !grid_bins(min_y0, max_y0, y0_step, ny) ||
      !grid_bins(min_radius, max_radius, radius_step, nr) ||
      !grid_bins(min_t, max_t, t_step, nt))
    return empty;
  if (static_cast<long long>(nx) * ny > 2000000LL ||
      static_cast<long long>(nr) * nt > 5000000LL) {
    std::cerr << "ERROR: Hough grid is too large; increase the grid steps"
              << std::endl;
    return empty;
  }

  std::vector<double> accumulator(static_cast<std::size_t>(nr) * nt);
  const int pool_limit = std::max(32, max_candidates * 32);
  std::vector<circle_t> peaks;
  peaks.reserve(pool_limit);

  for (int ix = 0; ix < nx; ++ix) {
    const double x0 = min_x0 + ix * x0_step;
    for (int iy = 0; iy < ny; ++iy) {
      const double y0 = min_y0 + iy * y0_step;
      std::fill(accumulator.begin(), accumulator.end(), 0.);

      for (int index : indices) {
        const double radius = std::hypot(points[index].x - x0,
                                         points[index].y - y0);
        const int rfirst = std::max(
            0, static_cast<int>(std::ceil(
                   (radius - gaussian_cut * spatial_resolution - min_radius) /
                   radius_step)));
        const int rlast = std::min(
            nr - 1, static_cast<int>(std::floor(
                        (radius + gaussian_cut * spatial_resolution -
                         min_radius) /
                        radius_step)));
        const int tfirst = std::max(
            0, static_cast<int>(std::ceil(
                   (points[index].time - gaussian_cut * time_resolution -
                    min_t) /
                   t_step)));
        const int tlast = std::min(
            nt - 1, static_cast<int>(std::floor(
                        (points[index].time + gaussian_cut * time_resolution -
                         min_t) /
                        t_step)));
        for (int ir = rfirst; ir <= rlast; ++ir) {
          const double radius_bin = min_radius + ir * radius_step;
          const double spatial_weight = std::exp(
              -0.5 * std::pow((radius - radius_bin) / spatial_resolution, 2)) /
              (std::sqrt(2. * pi) * spatial_resolution);
          for (int it = tfirst; it <= tlast; ++it) {
            const double time_bin = min_t + it * t_step;
            const double time_weight = std::exp(
                -0.5 * std::pow((points[index].time - time_bin) /
                                 time_resolution, 2)) /
                (std::sqrt(2. * pi) * time_resolution);
            accumulator[static_cast<std::size_t>(ir) * nt + it] +=
                spatial_weight * time_weight;
          }
        }
      }

      for (int ir = 0; ir < nr; ++ir) {
        for (int it = 0; it < nt; ++it) {
          const double score = accumulator[static_cast<std::size_t>(ir) * nt + it];
          if (!(score > 0.))
            continue;
          bool local_maximum = true;
          for (int dr = -1; dr <= 1 && local_maximum; ++dr) {
            const int neighbor_r = ir + dr;
            if (neighbor_r < 0 || neighbor_r >= nr)
              continue;
            for (int dt = -1; dt <= 1; ++dt) {
              const int neighbor_t = it + dt;
              if (neighbor_t < 0 || neighbor_t >= nt ||
                  (dr == 0 && dt == 0))
                continue;
              if (accumulator[static_cast<std::size_t>(neighbor_r) * nt +
                              neighbor_t] > score) {
                local_maximum = false;
                break;
              }
            }
          }
          if (!local_maximum)
            continue;

          circle_t peak;
          peak.x = min_x0 + ix * x0_step;
          peak.y = min_y0 + iy * y0_step;
          peak.radius = min_radius + ir * radius_step;
          peak.ring_time = min_t + it * t_step;
          peak.score = score;
          peaks.push_back(peak);
        }
      }
      if (static_cast<int>(peaks.size()) > 2 * pool_limit) {
        std::nth_element(
            peaks.begin(), peaks.begin() + pool_limit, peaks.end(),
            [](const circle_t &a, const circle_t &b) {
              return a.score > b.score;
            });
        peaks.resize(pool_limit);
      }
    }
  }

  std::sort(peaks.begin(), peaks.end(),
            [](const circle_t &a, const circle_t &b) {
              return a.score > b.score;
            });
  std::vector<circle_t> candidates;
  candidates.reserve(max_candidates);
  for (const circle_t &peak : peaks) {
    bool duplicate = false;
    for (const circle_t &candidate : candidates) {
      if (std::abs(peak.x - candidate.x) <= x0_step &&
          std::abs(peak.y - candidate.y) <= y0_step &&
          std::abs(peak.radius - candidate.radius) <= radius_step &&
          std::abs(peak.ring_time - candidate.ring_time) <= t_step) {
        duplicate = true;
        break;
      }
    }
    if (duplicate)
      continue;
    candidates.push_back(peak);
    if (static_cast<int>(candidates.size()) == max_candidates)
      break;
  }
  return candidates;
}

ring_hough_cuda::grid_t
make_cuda_grid(double min_x0, double max_x0, double x0_step,
               double min_y0, double max_y0, double y0_step,
               double min_radius, double max_radius, double radius_step,
               double min_t, double max_t, double t_step,
               double spatial_resolution, double time_resolution)
{
  ring_hough_cuda::grid_t grid{};
  grid.min_x0 = min_x0;
  grid.min_y0 = min_y0;
  grid.min_radius = min_radius;
  grid.min_t = min_t;
  grid.x0_step = x0_step;
  grid.y0_step = y0_step;
  grid.radius_step = radius_step;
  grid.t_step = t_step;
  grid.spatial_resolution = spatial_resolution;
  grid.time_resolution = time_resolution;
  grid_bins(min_x0, max_x0, x0_step, grid.nx);
  grid_bins(min_y0, max_y0, y0_step, grid.ny);
  grid_bins(min_radius, max_radius, radius_step, grid.nr);
  grid_bins(min_t, max_t, t_step, grid.nt);
  return grid;
}

bool
hough_candidates_cuda(const std::vector<point_t> &points,
                       const std::vector<int> &indices,
                       int max_candidates,
                       ring_hough_cuda::engine_t &engine,
                       const ring_hough_cuda::grid_t &grid,
                       std::vector<circle_t> &candidates,
                       std::string &error)
{
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> time;
  x.reserve(indices.size());
  y.reserve(indices.size());
  time.reserve(indices.size());
  for (const int index : indices) {
    x.push_back(static_cast<float>(points[index].x));
    y.push_back(static_cast<float>(points[index].y));
    time.push_back(static_cast<float>(points[index].time));
  }

  std::vector<ring_hough_cuda::candidate_t> results(max_candidates);
  int result_count = 0;
  if (!engine.find(x.data(), y.data(), time.data(),
                  static_cast<int>(indices.size()), max_candidates,
                  results.data(), result_count, error)) {
    return false;
  }
  candidates.clear();
  for (int iresult = 0; iresult < result_count; ++iresult) {
    const auto &result = results[iresult];
    circle_t candidate;
    candidate.x = grid.min_x0 + result.x0_bin * grid.x0_step;
    candidate.y = grid.min_y0 + result.y0_bin * grid.y0_step;
    candidate.radius = grid.min_radius + result.radius_bin * grid.radius_step;
    candidate.ring_time = grid.min_t + result.t_bin * grid.t_step;
    candidate.score = result.score;
    candidates.push_back(candidate);
  }
  return true;
}

double
hough_score(const std::vector<point_t> &points,
            const std::vector<int> &indices,
            const circle_t &circle,
            double spatial_resolution, double time_resolution)
{
  double score = 0.;
  for (const int index : indices) {
    const double spatial = radial_residual(circle, points[index]);
    const double temporal = points[index].time - circle.ring_time;
    if (std::abs(spatial) > gaussian_cut * spatial_resolution ||
        std::abs(temporal) > gaussian_cut * time_resolution)
      continue;
    const double spatial_pull = spatial / spatial_resolution;
    const double time_pull = temporal / time_resolution;
    score += std::exp(-0.5 * spatial_pull * spatial_pull) /
             (std::sqrt(2. * pi) * spatial_resolution) *
             std::exp(-0.5 * time_pull * time_pull) /
             (std::sqrt(2. * pi) * time_resolution);
  }
  return score;
}

bool
interpolate_peak(const std::vector<point_t> &points,
                 const std::vector<int> &indices,
                 double min_x0, double max_x0, double x0_step,
                 double min_y0, double max_y0, double y0_step,
                 double min_radius, double max_radius, double radius_step,
                 double min_t, double max_t, double t_step,
                 double spatial_resolution, double time_resolution,
                 circle_t &candidate)
{
  const double eps = 1.e-9;
  if (candidate.x - x0_step < min_x0 - eps ||
      candidate.x + x0_step > max_x0 + eps ||
      candidate.y - y0_step < min_y0 - eps ||
      candidate.y + y0_step > max_y0 + eps ||
      candidate.radius - radius_step < min_radius - eps ||
      candidate.radius + radius_step > max_radius + eps ||
      candidate.ring_time - t_step < min_t - eps ||
      candidate.ring_time + t_step > max_t + eps)
    return false;

  constexpr int nterms = 15;
  double normal[nterms][nterms + 1] = {};
  for (int ux = -1; ux <= 1; ++ux) {
    for (int uy = -1; uy <= 1; ++uy) {
      for (int ur = -1; ur <= 1; ++ur) {
        for (int ut = -1; ut <= 1; ++ut) {
          const double u[4] = {static_cast<double>(ux),
                               static_cast<double>(uy),
                               static_cast<double>(ur),
                               static_cast<double>(ut)};
          circle_t sample = candidate;
          sample.x += ux * x0_step;
          sample.y += uy * y0_step;
          sample.radius += ur * radius_step;
          sample.ring_time += ut * t_step;
          const double value = hough_score(
              points, indices, sample, spatial_resolution, time_resolution);
          const double feature[nterms] = {
              1., u[0], u[1], u[2], u[3],
              u[0] * u[0], u[1] * u[1], u[2] * u[2], u[3] * u[3],
              u[0] * u[1], u[0] * u[2], u[0] * u[3],
              u[1] * u[2], u[1] * u[3], u[2] * u[3]};
          for (int row = 0; row < nterms; ++row) {
            for (int column = 0; column < nterms; ++column)
              normal[row][column] += feature[row] * feature[column];
            normal[row][nterms] += feature[row] * value;
          }
        }
      }
    }
  }

  double coefficients[nterms] = {};
  for (int column = 0; column < nterms; ++column) {
    int pivot = column;
    for (int row = column + 1; row < nterms; ++row) {
      if (std::abs(normal[row][column]) >
          std::abs(normal[pivot][column]))
        pivot = row;
    }
    if (std::abs(normal[pivot][column]) < 1.e-12)
      return false;
    if (pivot != column) {
      for (int j = column; j <= nterms; ++j)
        std::swap(normal[column][j], normal[pivot][j]);
    }
    for (int row = column + 1; row < nterms; ++row) {
      const double factor = normal[row][column] / normal[column][column];
      for (int j = column; j <= nterms; ++j)
        normal[row][j] -= factor * normal[column][j];
    }
  }
  for (int row = nterms - 1; row >= 0; --row) {
    double value = normal[row][nterms];
    for (int column = row + 1; column < nterms; ++column)
      value -= normal[row][column] * coefficients[column];
    coefficients[row] = value / normal[row][row];
  }

  double gradient[4] = {coefficients[1], coefficients[2],
                        coefficients[3], coefficients[4]};
  double hessian[4][4] = {};
  for (int i = 0; i < 4; ++i)
    hessian[i][i] = 2. * coefficients[5 + i];
  const int cross_terms[4][4] = {
      {-1, 9, 10, 11},
      {9, -1, 12, 13},
      {10, 12, -1, 14},
      {11, 13, 14, -1}};
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j)
      hessian[i][j] = hessian[j][i] = coefficients[cross_terms[i][j]];
  }

  // A concave local quadratic has -H positive definite. Cholesky also
  // provides a stable solve for the sub-bin displacement.
  double lower[4][4] = {};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j <= i; ++j) {
      double value = -hessian[i][j];
      for (int k = 0; k < j; ++k)
        value -= lower[i][k] * lower[j][k];
      if (i == j) {
        if (!(value > 1.e-9) || !std::isfinite(value))
          return false;
        lower[i][j] = std::sqrt(value);
      } else {
        lower[i][j] = value / lower[j][j];
      }
    }
  }

  double forward[4] = {};
  for (int i = 0; i < 4; ++i) {
    double value = gradient[i];
    for (int j = 0; j < i; ++j)
      value -= lower[i][j] * forward[j];
    forward[i] = value / lower[i][i];
  }
  double displacement[4] = {};
  for (int i = 3; i >= 0; --i) {
    double value = forward[i];
    for (int j = i + 1; j < 4; ++j)
      value -= lower[j][i] * displacement[j];
    displacement[i] = value / lower[i][i];
  }
  for (const double value : displacement) {
    if (!std::isfinite(value) || std::abs(value) > 1.)
      return false;
  }

  circle_t refined = candidate;
  refined.x += displacement[0] * x0_step;
  refined.y += displacement[1] * y0_step;
  refined.radius += displacement[2] * radius_step;
  refined.ring_time += displacement[3] * t_step;
  if (refined.radius <= 0. || refined.x < min_x0 || refined.x > max_x0 ||
      refined.y < min_y0 || refined.y > max_y0 ||
      refined.radius < min_radius || refined.radius > max_radius ||
      refined.ring_time < min_t || refined.ring_time > max_t)
    return false;
  refined.score = hough_score(points, indices, refined,
                              spatial_resolution, time_resolution);
  if (!std::isfinite(refined.score))
    return false;
  candidate = refined;
  return true;
}

void
evaluate_candidate(const std::vector<point_t> &points,
                    const std::vector<int> &indices,
                    int min_inliers,
                    double spatial_resolution, double time_resolution,
                    circle_t &circle,
                    std::vector<int> &inlier_indices)
{
  inlier_indices.clear();
  int inliers = 0;
  double residual_sum = 0.;
  for (int index : indices) {
    const double spatial = std::abs(radial_residual(circle, points[index]));
    const double temporal = std::abs(points[index].time - circle.ring_time);
    if (spatial <= gaussian_cut * spatial_resolution &&
        temporal <= gaussian_cut * time_resolution) {
      ++inliers;
      residual_sum += spatial;
      inlier_indices.push_back(index);
    }
  }
  circle.inliers = inliers;
  circle.residual_sum = residual_sum;
  circle.valid = inliers >= min_inliers;
}

int
count_shared_hits(const std::vector<int> &first,
                  const std::vector<int> &second)
{
  int shared = 0;
  std::size_t i = 0;
  std::size_t j = 0;
  while (i < first.size() && j < second.size()) {
    if (first[i] == second[j]) {
      ++shared;
      ++i;
      ++j;
    } else if (first[i] < second[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  return shared;
}

bool
ring_finder_hough(const std::string &filename, const std::string &outfilename,
                  int min_inliers, int max_rings,
                  int max_shared_hits,
                  double min_x0, double max_x0, double x0_step,
                  double min_y0, double max_y0, double y0_step,
                  double min_radius, double max_radius, double radius_step,
                  double min_t, double max_t, double t_step,
                  double spatial_resolution, double time_resolution,
                  int ransac_iterations,
                  double ransac_center_window,
                  double ransac_radius_window,
                  double ransac_time_window,
                  double ransac_tolerance,
                  bool use_gpu)
{
  auto fin = TFile::Open(filename.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << filename << std::endl;
    return false;
  }
  auto frames_in = (TTree *)fin->Get("frames");
  auto cherenkov_in = (TTree *)fin->Get("cherenkov");
  if (!frames_in || !cherenkov_in) {
    std::cerr << "ERROR: input must contain 'frames' and 'cherenkov' trees"
              << std::endl;
    fin->Close();
    return false;
  }
  for (const char *name : {"nhits", "x", "y", "time"}) {
    if (!has_branch(cherenkov_in, name)) {
      fin->Close();
      return false;
    }
  }
  if (fin->Get("ring")) {
    std::cerr << "ERROR: input already contains a 'ring' tree" << std::endl;
    fin->Close();
    return false;
  }

  std::unique_ptr<ring_hough_cuda::engine_t> gpu;
  if (use_gpu) {
#ifdef RING_HOUGH_HAS_CUDA
    gpu = std::make_unique<ring_hough_cuda::engine_t>();
#else
    std::cerr << "ERROR: this build has no CUDA Hough backend; omit --gpu or "
                 "rebuild with CUDA available"
              << std::endl;
    fin->Close();
    return false;
#endif
  }

  const Long64_t entries = frames_in->GetEntries();
  if (cherenkov_in->GetEntries() != entries) {
    std::cerr << "ERROR: frames/cherenkov entry-count mismatch: frames="
              << entries << " cherenkov=" << cherenkov_in->GetEntries()
              << std::endl;
    fin->Close();
    return false;
  }

  TTreeReader reader(cherenkov_in);
  TTreeReaderValue<UShort_t> ncherenkovhits(reader, "nhits");
  TTreeReaderArray<Float_t> x(reader, "x");
  TTreeReaderArray<Float_t> y(reader, "y");
  TTreeReaderArray<Float_t> time(reader, "time");

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename
              << std::endl;
    return false;
  }
  fout->cd();
  auto frames_out = frames_in->CloneTree(0);
  auto trigger_in = (TTree *)fin->Get("trigger");
  auto timing_in = (TTree *)fin->Get("timing");
  auto trigger_out = trigger_in ? trigger_in->CloneTree(0) : nullptr;
  auto timing_out = timing_in ? timing_in->CloneTree(0) : nullptr;
  auto cherenkov_out = cherenkov_in->CloneTree(0);
  if (!frames_out || !cherenkov_out || (trigger_in && !trigger_out) ||
      (timing_in && !timing_out)) {
    std::cerr << "ERROR: failed to create output frame trees" << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  UChar_t nring = 0;
  float ring_x0[maxrings];
  float ring_y0[maxrings];
  float ring_r[maxrings];
  float ring_e[maxrings];
  float ring_phi[maxrings];
  float ring_time[maxrings];
  UShort_t ring_ninliers[maxrings];
  auto ring_out = new TTree("ring", "ring candidates, one entry per frame");
  ring_out->Branch("nring", &nring, "nring/b");
  ring_out->Branch("ring_x0", ring_x0, "ring_x0[nring]/F");
  ring_out->Branch("ring_y0", ring_y0, "ring_y0[nring]/F");
  ring_out->Branch("ring_r", ring_r, "ring_r[nring]/F");
  ring_out->Branch("ring_e", ring_e, "ring_e[nring]/F");
  ring_out->Branch("ring_phi", ring_phi, "ring_phi[nring]/F");
  ring_out->Branch("ring_time", ring_time, "ring_time[nring]/F");
  ring_out->Branch("ring_ninliers", ring_ninliers, "ring_ninliers[nring]/s");

  Long64_t frames = 0;
  Long64_t accepted = 0;
  // Keep several maxima available for validation without making the CUDA
  // reduction scan the full accumulator an excessive number of times.
  const int candidate_limit = std::max(8, max_rings * 4);
  std::mt19937 generator(0x9e3779b9u);
  for (Long64_t iframe = 0; iframe < entries; ++iframe) {
    if (frames_in->GetEntry(iframe) <= 0 || !reader.Next()) {
      std::cerr << "ERROR: failed to read synchronized frame entry " << iframe
                << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    ++frames;
    nring = 0;
    const int nhits = static_cast<int>(*ncherenkovhits);
    if (x.GetSize() < nhits || y.GetSize() < nhits || time.GetSize() < nhits) {
      std::cerr << "ERROR: invalid Cherenkov hit array at frame entry "
                << iframe << " nhits=" << nhits << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    std::vector<point_t> points;
    std::vector<int> indices;
    for (int index = 0; index < nhits; ++index) {
      if (std::isfinite(x[index]) && std::isfinite(y[index]) &&
          std::isfinite(time[index])) {
        points.push_back({x[index], y[index], time[index]});
        indices.push_back(static_cast<int>(points.size() - 1));
      }
    }

    const std::vector<circle_t> seeds = ransac_seeds(
        points, indices, max_rings, ransac_iterations, min_inliers,
        min_x0, max_x0, min_y0, max_y0, min_radius, max_radius,
        ransac_tolerance, ransac_time_window, generator);

    double scan_min_x0 = min_x0;
    double scan_max_x0 = max_x0;
    double scan_min_y0 = min_y0;
    double scan_max_y0 = max_y0;
    double scan_min_radius = min_radius;
    double scan_max_radius = max_radius;
    double scan_min_t = min_t;
    double scan_max_t = max_t;
    if (!seeds.empty()) {
      scan_min_x0 = max_x0;
      scan_max_x0 = min_x0;
      scan_min_y0 = max_y0;
      scan_max_y0 = min_y0;
      scan_min_radius = max_radius;
      scan_max_radius = min_radius;
      scan_min_t = max_t;
      scan_max_t = min_t;
      for (const circle_t &seed : seeds) {
        scan_min_x0 = std::min(scan_min_x0, seed.x - ransac_center_window);
        scan_max_x0 = std::max(scan_max_x0, seed.x + ransac_center_window);
        scan_min_y0 = std::min(scan_min_y0, seed.y - ransac_center_window);
        scan_max_y0 = std::max(scan_max_y0, seed.y + ransac_center_window);
        scan_min_radius = std::min(scan_min_radius,
                                   seed.radius - ransac_radius_window);
        scan_max_radius = std::max(scan_max_radius,
                                   seed.radius + ransac_radius_window);
        scan_min_t = std::min(scan_min_t, seed.ring_time - ransac_time_window);
        scan_max_t = std::max(scan_max_t, seed.ring_time + ransac_time_window);
      }
      scan_min_x0 = std::max(min_x0, scan_min_x0);
      scan_max_x0 = std::min(max_x0, scan_max_x0);
      scan_min_y0 = std::max(min_y0, scan_min_y0);
      scan_max_y0 = std::min(max_y0, scan_max_y0);
      scan_min_radius = std::max(min_radius, scan_min_radius);
      scan_max_radius = std::min(max_radius, scan_max_radius);
      scan_min_t = std::max(min_t, scan_min_t);
      scan_max_t = std::min(max_t, scan_max_t);
      align_grid_bounds(min_x0, max_x0, x0_step, scan_min_x0, scan_max_x0);
      align_grid_bounds(min_y0, max_y0, y0_step, scan_min_y0, scan_max_y0);
      align_grid_bounds(min_radius, max_radius, radius_step,
                        scan_min_radius, scan_max_radius);
      align_grid_bounds(min_t, max_t, t_step, scan_min_t, scan_max_t);
    }

    std::vector<std::vector<int>> accepted_inliers;
    std::vector<circle_t> candidates;
    // The fine scan is intentionally seed-local in the default mode. With
    // quarter-unit grid steps, falling back to the full global volume for a
    // seedless frame would require an impractical accumulator. A global scan
    // remains available explicitly with --ransac-iterations 0.
    const bool scan_allowed = !seeds.empty() || ransac_iterations == 0;
    if (indices.size() >= 3 && scan_allowed && gpu) {
      const auto local_grid = make_cuda_grid(
          scan_min_x0, scan_max_x0, x0_step,
          scan_min_y0, scan_max_y0, y0_step,
          scan_min_radius, scan_max_radius, radius_step,
          scan_min_t, scan_max_t, t_step,
          spatial_resolution, time_resolution);
      std::string error;
      if (!gpu->initialize(local_grid, error)) {
        std::cerr << "ERROR: could not initialize CUDA Hough backend: "
                  << error << std::endl;
        fout->Close();
        fin->Close();
        return false;
      }
      if (!hough_candidates_cuda(points, indices, candidate_limit, *gpu,
                                 local_grid, candidates, error)) {
        std::cerr << "ERROR: CUDA Hough scan failed: " << error << std::endl;
        fout->Close();
        fin->Close();
        return false;
      }
    } else if (indices.size() >= 3 && scan_allowed) {
      candidates = hough_candidates(
          points, indices, candidate_limit,
          scan_min_x0, scan_max_x0, x0_step,
          scan_min_y0, scan_max_y0, y0_step,
          scan_min_radius, scan_max_radius, radius_step,
          scan_min_t, scan_max_t, t_step,
          spatial_resolution, time_resolution);
    }

    for (circle_t &candidate : candidates)
      interpolate_peak(points, indices,
                       scan_min_x0, scan_max_x0, x0_step,
                       scan_min_y0, scan_max_y0, y0_step,
                       scan_min_radius, scan_max_radius, radius_step,
                       scan_min_t, scan_max_t, t_step,
                       spatial_resolution, time_resolution, candidate);

    for (circle_t ring : candidates) {
      std::vector<int> inlier_indices;
      evaluate_candidate(points, indices, min_inliers,
                         spatial_resolution, time_resolution, ring,
                         inlier_indices);
      if (!ring.valid)
        continue;

      bool too_many_shared_hits = false;
      for (const auto &previous_inliers : accepted_inliers) {
        if (count_shared_hits(inlier_indices, previous_inliers) >
            max_shared_hits) {
          too_many_shared_hits = true;
          break;
        }
      }
      if (too_many_shared_hits)
        continue;

      ring_x0[nring] = static_cast<float>(ring.x);
      ring_y0[nring] = static_cast<float>(ring.y);
      ring_r[nring] = static_cast<float>(ring.radius);
      ring_e[nring] = 0.;
      ring_phi[nring] = 0.;
      ring_time[nring] = static_cast<float>(ring.ring_time);
      ring_ninliers[nring] = static_cast<UShort_t>(ring.inliers);
      ++nring;
      ++accepted;
      accepted_inliers.push_back(std::move(inlier_indices));
      if (nring == max_rings)
        break;
    }

    if (trigger_in)
      trigger_in->GetEntry(iframe);
    if (timing_in)
      timing_in->GetEntry(iframe);
    cherenkov_in->GetEntry(iframe);
    frames_out->Fill();
    if (trigger_out)
      trigger_out->Fill();
    if (timing_out)
      timing_out->Fill();
    cherenkov_out->Fill();
    ring_out->Fill();
  }

  for (auto tree : {frames_out, trigger_out, timing_out, cherenkov_out,
                    ring_out}) {
    if (tree && tree->GetEntries() != entries) {
      std::cerr << "ERROR: output tree entry-count mismatch: "
                << tree->GetName() << "=" << tree->GetEntries()
                << " expected=" << entries << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }
  }
  auto participation = (TTree *)fin->Get("spill_participation");
  if (participation) {
    fout->cd();
    if (!participation->CloneTree(-1, "fast")) {
      std::cerr << "ERROR: failed to copy spill_participation tree" << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }
  }
  fout->Write();
  fout->Close();
  fin->Close();
  std::cout << "frames processed: " << frames << std::endl
            << "rings found:      " << accepted << std::endl
            << "output:            " << outfilename << std::endl;
  return true;
}

}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;
  std::string input, output;
  int min_inliers = 8, max_rings = maxrings;
  int max_shared_hits = 2;
  double min_x0 = -100., max_x0 = 100., x0_step = 0.25;
  double min_y0 = -100., max_y0 = 100., y0_step = 0.25;
  double min_radius = 1., max_radius = 200., radius_step = 0.25;
  double min_t = -32., max_t = 32., t_step = 0.25;
  double spatial_resolution = 1.5;
  double time_resolution = 1.;
  int ransac_iterations = 128;
  double ransac_center_window = 10.;
  double ransac_radius_window = 10.;
  double ransac_time_window = 5.;
  double ransac_tolerance = 5.;
  bool use_gpu = false;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input triggered ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("min-inliers", po::value<int>(&min_inliers)->default_value(min_inliers), "minimum spatial/time inliers")
    ("max-rings", po::value<int>(&max_rings)->default_value(max_rings), "maximum rings per frame")
    ("max-shared-hits", po::value<int>(&max_shared_hits)
                             ->default_value(max_shared_hits),
     "maximum hits shared by any two accepted rings")
    ("min-x0", po::value<double>(&min_x0)->default_value(min_x0), "minimum ring center x")
    ("max-x0", po::value<double>(&max_x0)->default_value(max_x0), "maximum ring center x")
    ("x0-step", po::value<double>(&x0_step)->default_value(x0_step), "Hough center x step")
    ("min-y0", po::value<double>(&min_y0)->default_value(min_y0), "minimum ring center y")
    ("max-y0", po::value<double>(&max_y0)->default_value(max_y0), "maximum ring center y")
    ("y0-step", po::value<double>(&y0_step)->default_value(y0_step), "Hough center y step")
    ("min-radius", po::value<double>(&min_radius)->default_value(min_radius), "minimum ring radius")
    ("max-radius", po::value<double>(&max_radius)->default_value(max_radius), "maximum ring radius")
    ("radius-step", po::value<double>(&radius_step)->default_value(radius_step), "Hough radius step")
    ("min-t", po::value<double>(&min_t)->default_value(min_t), "minimum ring time")
    ("max-t", po::value<double>(&max_t)->default_value(max_t), "maximum ring time")
    ("t-step", po::value<double>(&t_step)->default_value(t_step), "Hough ring-time step")
    ("spatial-resolution", po::value<double>(&spatial_resolution)
                               ->default_value(spatial_resolution),
     "Gaussian spatial resolution for Hough weights")
    ("time-resolution", po::value<double>(&time_resolution)
                            ->default_value(time_resolution),
     "Gaussian time resolution for Hough weights")
    ("ransac-iterations", po::value<int>(&ransac_iterations)
                              ->default_value(ransac_iterations),
     "RANSAC iterations used to localize the fine Hough scan")
    ("ransac-center-window", po::value<double>(&ransac_center_window)
                                   ->default_value(ransac_center_window),
     "half-width of the local Hough x0/y0 window around RANSAC seeds")
    ("ransac-radius-window", po::value<double>(&ransac_radius_window)
                                   ->default_value(ransac_radius_window),
     "half-width of the local Hough radius window around RANSAC seeds")
    ("ransac-time-window", po::value<double>(&ransac_time_window)
                                 ->default_value(ransac_time_window),
     "RANSAC time window and local Hough ring-time half-width")
    ("ransac-tolerance", po::value<double>(&ransac_tolerance)
                                ->default_value(ransac_tolerance),
     "RANSAC spatial inlier tolerance in mm")
    ("gpu", po::bool_switch(&use_gpu)->default_value(false),
     "use the CUDA Hough backend (if compiled)")
    ;
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) {
      std::cout << options << std::endl;
      return 0;
    }
    po::notify(vm);
    int unused = 0;
    if (min_inliers < 3 || max_rings < 1 || max_rings > maxrings ||
        max_shared_hits < 0 ||
        ransac_iterations < 0 || ransac_center_window <= 0. ||
        ransac_radius_window <= 0. || ransac_time_window <= 0. ||
        ransac_tolerance <= 0. ||
        min_x0 > max_x0 || min_y0 > max_y0 || min_radius <= 0. ||
        min_radius > max_radius ||
        min_t > max_t || spatial_resolution <= 0. || time_resolution <= 0. ||
        !grid_bins(min_x0, max_x0, x0_step, unused) ||
        !grid_bins(min_y0, max_y0, y0_step, unused) ||
        !grid_bins(min_radius, max_radius, radius_step, unused) ||
        !grid_bins(min_t, max_t, t_step, unused))
      throw std::runtime_error("invalid ring-finder-hough parameter");
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << std::endl
              << options << std::endl;
    return 1;
  }

  return ring_finder_hough(
             input, output, min_inliers, max_rings, max_shared_hits,
             min_x0, max_x0, x0_step, min_y0, max_y0, y0_step,
             min_radius, max_radius, radius_step, min_t, max_t, t_step,
             spatial_resolution, time_resolution,
             ransac_iterations, ransac_center_window,
             ransac_radius_window, ransac_time_window, ransac_tolerance,
             use_gpu)
             ? 0 : 1;
}
