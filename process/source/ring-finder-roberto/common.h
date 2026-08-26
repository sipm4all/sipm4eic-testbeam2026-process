#pragma once

struct points_t {
  int n;
  float *x;
  float *y;
  float *t;
};

struct map_t {
  float *x;
  float *y;
  float *r;
  float *t;
};

struct hough_t {
  float *h;
  float *rh;
  int *rhi;
  int *nh;
};

struct range_t {
  float x;
  float y;
  float r;
  float t;
};

struct sigma_t {
  float x;
  float t;
};

struct bins_t {
  int x;
  int y;
  int r;
  int t;
};

struct data_t {
  points_t points;
  map_t map;
  hough_t hough;
  range_t min;
  range_t max;
  bins_t bins;
  sigma_t sigma;
};

void hough_init(data_t data);
void hough_transform(data_t data);
void hough_free();
