#include "../include/TimingEstimator.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "TimingEstimatorData.inc"

namespace {

TimingEstimator::tree_model_t
make_model(double baseline,
           int ntrees,
           const int *offset,
           const int *feature,
           const int *left,
           const int *right,
           const unsigned char *is_leaf,
           const unsigned char *missing_left,
           const double *threshold,
           const double *value)
{
  return TimingEstimator::tree_model_t{
    baseline, ntrees, offset, feature, left, right,
    is_leaf, missing_left, threshold, value
  };
}

const TimingEstimator::tree_model_t timeA0 =
  make_model(timing_model_data::timeA0_baseline,
             timing_model_data::timeA0_ntrees,
             timing_model_data::timeA0_offset,
             timing_model_data::timeA0_feature,
             timing_model_data::timeA0_left,
             timing_model_data::timeA0_right,
             timing_model_data::timeA0_is_leaf,
             timing_model_data::timeA0_missing_left,
             timing_model_data::timeA0_threshold,
             timing_model_data::timeA0_value);

const TimingEstimator::tree_model_t timeB0 =
  make_model(timing_model_data::timeB0_baseline,
             timing_model_data::timeB0_ntrees,
             timing_model_data::timeB0_offset,
             timing_model_data::timeB0_feature,
             timing_model_data::timeB0_left,
             timing_model_data::timeB0_right,
             timing_model_data::timeB0_is_leaf,
             timing_model_data::timeB0_missing_left,
             timing_model_data::timeB0_threshold,
             timing_model_data::timeB0_value);

const TimingEstimator::tree_model_t timeA1 =
  make_model(timing_model_data::timeA1_baseline,
             timing_model_data::timeA1_ntrees,
             timing_model_data::timeA1_offset,
             timing_model_data::timeA1_feature,
             timing_model_data::timeA1_left,
             timing_model_data::timeA1_right,
             timing_model_data::timeA1_is_leaf,
             timing_model_data::timeA1_missing_left,
             timing_model_data::timeA1_threshold,
             timing_model_data::timeA1_value);

const TimingEstimator::tree_model_t timeB1 =
  make_model(timing_model_data::timeB1_baseline,
             timing_model_data::timeB1_ntrees,
             timing_model_data::timeB1_offset,
             timing_model_data::timeB1_feature,
             timing_model_data::timeB1_left,
             timing_model_data::timeB1_right,
             timing_model_data::timeB1_is_leaf,
             timing_model_data::timeB1_missing_left,
             timing_model_data::timeB1_threshold,
             timing_model_data::timeB1_value);

} // namespace

void
TimingEstimator::anchor_and_features(const double *times,
                                     std::array<float, 62> &features_a,
                                     std::array<float, 52> &features_b,
                                     double &anchor)
{
  std::array<double, 32> sorted;
  for (int i = 0; i < 32; ++i)
    sorted[i] = times[i];
  std::sort(sorted.begin(), sorted.end());

  anchor = 0.;
  for (int i = 0; i < 10; ++i)
    anchor += sorted[i];
  anchor /= 10.;

  for (int i = 0; i < 32; ++i) {
    float rel = static_cast<float>(times[i] - anchor);
    features_a[i] = rel;
    features_b[i] = rel;
  }

  for (int i = 0; i < 24; ++i)
    features_a[32 + i] = static_cast<float>(sorted[i] - anchor);
  for (int i = 0; i < 20; ++i)
    features_b[32 + i] = static_cast<float>(sorted[i] - anchor);

  features_a[56] = static_cast<float>(sorted[0] - anchor);
  features_a[57] = static_cast<float>(sorted[1] - sorted[0]);
  features_a[58] = static_cast<float>(sorted[3] - sorted[0]);
  features_a[59] = static_cast<float>(sorted[7] - sorted[0]);
  features_a[60] = static_cast<float>(sorted[15] - sorted[0]);
  features_a[61] = static_cast<float>(sorted[23] - sorted[0]);
}

void
TimingEstimator::scale_features(const std::array<float, 52> &in,
                                const double *mean,
                                const double *scale,
                                std::array<float, 52> &out)
{
  for (int i = 0; i < 52; ++i)
    out[i] = static_cast<float>((static_cast<double>(in[i]) - mean[i]) / scale[i]);
}

double
TimingEstimator::predict_tree_model(const tree_model_t &model, const float *features)
{
  double out = model.baseline;

  for (int itree = 0; itree < model.ntrees; ++itree) {
    int base = model.offset[itree];
    int idx = base;
    for (;;) {
      if (model.is_leaf[idx]) {
        out += model.value[idx];
        break;
      }

      int feature = model.feature[idx];
      float value = features[feature];
      bool go_left = false;
      if (std::isnan(value))
        go_left = model.missing_left[idx] != 0;
      else
        go_left = static_cast<double>(value) <= model.threshold[idx];

      idx = base + (go_left ? model.left[idx] : model.right[idx]);
    }
  }

  return out;
}

float
TimingEstimator::silu(float x)
{
  return x / (1.f + std::exp(-x));
}

float
TimingEstimator::softplus(float x)
{
  if (x > 20.f)
    return x;
  if (x < -20.f)
    return std::exp(x);
  return std::log1p(std::exp(x));
}

double
TimingEstimator::predict_sigma_net(int detector, const std::array<float, 52> &features)
{
  const float *w0 = detector == 0 ? timing_model_data::net0_w0 : timing_model_data::net1_w0;
  const float *b0 = detector == 0 ? timing_model_data::net0_b0 : timing_model_data::net1_b0;
  const float *w1 = detector == 0 ? timing_model_data::net0_w1 : timing_model_data::net1_w1;
  const float *b1 = detector == 0 ? timing_model_data::net0_b1 : timing_model_data::net1_b1;
  const float *w2 = detector == 0 ? timing_model_data::net0_w2 : timing_model_data::net1_w2;
  const float *b2 = detector == 0 ? timing_model_data::net0_b2 : timing_model_data::net1_b2;

  std::array<float, 64> h0;
  std::array<float, 32> h1;

  for (int i = 0; i < 64; ++i) {
    float z = b0[i];
    for (int j = 0; j < 52; ++j)
      z += w0[i * 52 + j] * features[j];
    h0[i] = silu(z);
  }

  for (int i = 0; i < 32; ++i) {
    float z = b1[i];
    for (int j = 0; j < 64; ++j)
      z += w1[i * 64 + j] * h0[j];
    h1[i] = silu(z);
  }

  float z = b2[0];
  for (int j = 0; j < 32; ++j)
    z += w2[j] * h1[j];

  return std::sqrt(static_cast<double>(softplus(z)) + 1.e-6) *
         timing_model_data::sigma_scale;
}

void
TimingEstimator::predict_detector(const double *times,
                                  int detector,
                                  double &time,
                                  double &sigma)
{
  std::array<float, 62> features_a;
  std::array<float, 52> features_b;
  std::array<float, 52> scaled_b;
  double anchor = 0.;

  anchor_and_features(times, features_a, features_b, anchor);

  if (detector == 0) {
    scale_features(features_b, timing_model_data::scaler0_mean,
                   timing_model_data::scaler0_scale, scaled_b);
    double pA = predict_tree_model(timeA0, features_a.data());
    double pB = predict_tree_model(timeB0, scaled_b.data());
    time = anchor - (timing_model_data::blend_a * pA +
                     (1. - timing_model_data::blend_a) * pB);
  } else {
    scale_features(features_b, timing_model_data::scaler1_mean,
                   timing_model_data::scaler1_scale, scaled_b);
    double pA = predict_tree_model(timeA1, features_a.data());
    double pB = predict_tree_model(timeB1, scaled_b.data());
    time = anchor + (timing_model_data::blend_a * pA +
                     (1. - timing_model_data::blend_a) * pB);
  }

  sigma = predict_sigma_net(detector, scaled_b);
}

TimingResult
TimingEstimator::estimate(const std::array<double, 64> &times) const
{
  TimingResult result;
  predict_detector(times.data(), 0, result.T0, result.sigma0);
  predict_detector(times.data() + 32, 1, result.T1, result.sigma1);
  result.T = 0.5 * (result.T0 + result.T1);
  result.sigmaT = 0.5 * std::sqrt(result.sigma0 * result.sigma0 +
                                  result.sigma1 * result.sigma1);
  return result;
}
