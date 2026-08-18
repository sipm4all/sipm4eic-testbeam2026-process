#pragma once

#include <array>

struct TimingResult {
  double T0;
  double sigma0;
  double T1;
  double sigma1;
  double T;
  double sigmaT;

  double T0_ns() const;
  double T1_ns() const;
  double T_ns() const;
  double sigma0_ps() const;
  double sigma1_ps() const;
  double sigmaT_ps() const;
};

class TimingEstimator {
public:
  static constexpr double native_unit_ns = 3.125;
  static constexpr double native_unit_ps = 3125.;

  TimingEstimator() = default;

  TimingResult estimate(const std::array<double, 64> &times) const;

  struct tree_model_t {
    double baseline;
    int ntrees;
    const int *offset;
    const int *feature;
    const int *left;
    const int *right;
    const unsigned char *is_leaf;
    const unsigned char *missing_left;
    const double *threshold;
    const double *value;
  };

private:
  static void anchor_and_features(const double *times,
                                  std::array<float, 62> &features_a,
                                  std::array<float, 52> &features_b,
                                  double &anchor);
  static void scale_features(const std::array<float, 52> &in,
                             const double *mean,
                             const double *scale,
                             std::array<float, 52> &out);
  static double predict_tree_model(const tree_model_t &model, const float *features);
  static float silu(float x);
  static float softplus(float x);
  static double predict_sigma_net(int detector, const std::array<float, 52> &features);
  static void predict_detector(const double *times,
                               int detector,
                               double &time,
                               double &sigma);
};

inline double
TimingResult::T0_ns() const
{
  return T0 * TimingEstimator::native_unit_ns;
}

inline double
TimingResult::T1_ns() const
{
  return T1 * TimingEstimator::native_unit_ns;
}

inline double
TimingResult::T_ns() const
{
  return T * TimingEstimator::native_unit_ns;
}

inline double
TimingResult::sigma0_ps() const
{
  return sigma0 * TimingEstimator::native_unit_ps;
}

inline double
TimingResult::sigma1_ps() const
{
  return sigma1 * TimingEstimator::native_unit_ps;
}

inline double
TimingResult::sigmaT_ps() const
{
  return sigmaT * TimingEstimator::native_unit_ps;
}
