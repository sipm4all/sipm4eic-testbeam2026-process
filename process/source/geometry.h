#pragma once

#include "data_word.h"

#include <array>
#include <map>
#include <vector>

struct geometry_t {
  static constexpr double timing_pitch = 3.5;

  static bool coordinate(const data_t &data, double &x, double &y)
  {
    x = 0.;
    y = 0.;

    if (!data.is_alcor_hit())
      return false;

    if (data.device == 200)
      return timing_coordinate(data, x, y);

    return cherenkov_coordinate(data, x, y);
  }

private:
  static int timing_do_channel(const data_t &data)
  {
    static const int eo2do[32] = {
      22, 20, 18, 16, 24, 26, 28, 30,
      25, 27, 29, 31, 23, 21, 19, 17,
      9,  11, 13, 15, 7,  5,  3,  1,
      6,  4,  2,  0,  8,  10, 12, 14
    };

    int eoch = data.pixel + 4 * data.column;
    if (eoch < 0 || eoch >= 32)
      return -1;
    return eo2do[eoch];
  }

  static bool timing_coordinate(const data_t &data, double &x, double &y)
  {
    int doch = timing_do_channel(data);
    if (doch < 0)
      return false;

    x = timing_pitch * (doch % 4);
    y = timing_pitch * (doch / 4);
    return true;
  }

  static const std::map<std::array<int, 2>, std::array<int, 2>> &pdu_matrix_map()
  {
    static const std::map<std::array<int, 2>, std::array<int, 2>> map = {
      {{192, 0}, {1, 1}}, {{192, 2}, {1, 2}}, {{192, 4}, {1, 3}}, {{192, 6}, {1, 4}},
      {{193, 0}, {2, 1}}, {{193, 2}, {2, 2}}, {{193, 4}, {2, 3}}, {{193, 6}, {2, 4}},
      {{194, 0}, {3, 1}}, {{194, 2}, {3, 2}}, {{194, 4}, {3, 3}}, {{194, 6}, {3, 4}},
      {{195, 0}, {4, 1}}, {{195, 2}, {4, 2}}, {{195, 4}, {4, 3}}, {{195, 6}, {4, 4}},
      {{196, 0}, {5, 1}}, {{196, 2}, {5, 2}}, {{196, 4}, {5, 3}}, {{196, 6}, {5, 4}},
      {{197, 0}, {6, 1}}, {{197, 2}, {6, 2}}, {{197, 4}, {6, 3}}, {{197, 6}, {6, 4}},
      {{198, 0}, {7, 1}}, {{198, 2}, {7, 2}}, {{198, 4}, {7, 3}}, {{198, 6}, {7, 4}},
      {{199, 0}, {8, 1}}, {{199, 2}, {8, 2}}, {{199, 4}, {8, 3}}, {{199, 6}, {8, 4}},
      {{192, 1}, {1, 1}}, {{192, 3}, {1, 2}}, {{192, 5}, {1, 3}}, {{192, 7}, {1, 4}},
      {{193, 1}, {2, 1}}, {{193, 3}, {2, 2}}, {{193, 5}, {2, 3}}, {{193, 7}, {2, 4}},
      {{194, 1}, {3, 1}}, {{194, 3}, {3, 2}}, {{194, 5}, {3, 3}}, {{194, 7}, {3, 4}},
      {{195, 1}, {4, 1}}, {{195, 3}, {4, 2}}, {{195, 5}, {4, 3}}, {{195, 7}, {4, 4}},
      {{196, 1}, {5, 1}}, {{196, 3}, {5, 2}}, {{196, 5}, {5, 3}}, {{196, 7}, {5, 4}},
      {{197, 1}, {6, 1}}, {{197, 3}, {6, 2}}, {{197, 5}, {6, 3}}, {{197, 7}, {6, 4}},
      {{198, 1}, {7, 1}}, {{198, 3}, {7, 2}}, {{198, 5}, {7, 3}}, {{198, 7}, {7, 4}},
      {{199, 1}, {8, 1}}, {{199, 3}, {8, 2}}, {{199, 5}, {8, 3}}, {{199, 7}, {8, 4}}
    };
    return map;
  }

  static const std::map<int, std::vector<int>> &matrix_mapping()
  {
    static const std::map<int, std::vector<int>> map = {
      {1, {3, 2, 1, 0, 8, 9, 10, 11, 17, 16, 12, 4, 18, 19, 5, 13,
           25, 24, 21, 20, 26, 27, 28, 29, 30, 22, 14, 6, 7, 15, 23, 31,
           35, 34, 33, 32, 36, 37, 38, 39, 43, 42, 41, 40, 44, 45, 46, 47,
           51, 50, 49, 48, 52, 53, 54, 55, 59, 58, 57, 56, 60, 61, 62, 63}},
      {2, {59, 58, 57, 56, 60, 61, 62, 63, 51, 50, 49, 48, 52, 53, 54, 55,
           43, 42, 41, 40, 44, 45, 46, 47, 35, 34, 33, 32, 36, 37, 38, 39,
           0, 8, 16, 24, 25, 17, 26, 27, 29, 28, 1, 9, 30, 31, 18, 19,
           21, 20, 2, 10, 22, 23, 11, 12, 3, 15, 14, 13, 4, 5, 6, 7}},
      {3, {4, 5, 6, 7, 3, 2, 1, 0, 12, 13, 14, 15, 11, 10, 9, 8,
           20, 21, 22, 23, 19, 18, 17, 16, 28, 29, 30, 31, 27, 26, 25, 24,
           46, 63, 55, 47, 54, 62, 39, 38, 34, 35, 36, 37, 33, 32, 45, 44,
           42, 43, 61, 53, 41, 40, 52, 60, 48, 49, 50, 51, 59, 58, 57, 56}},
      {4, {60, 61, 62, 63, 55, 54, 53, 52, 46, 47, 51, 59, 45, 44, 58, 50,
           38, 39, 42, 43, 37, 36, 35, 34, 33, 41, 49, 57, 56, 48, 40, 32,
           28, 29, 30, 31, 27, 26, 25, 24, 20, 21, 22, 23, 19, 18, 17, 16,
           12, 13, 14, 15, 11, 10, 9, 8, 4, 5, 6, 7, 3, 2, 1, 0}}
    };
    return map;
  }

  static const std::map<int, std::array<double, 2>> &placement_xy()
  {
    static const std::map<int, std::array<double, 2>> map = {
      {1, {-82.,  30.}}, {2, {-26.,  35.}}, {3, {30.,  30.}},
      {8, {-82., -26.}},                      {4, {30., -26.}},
      {7, {-82., -82.}}, {6, {-26., -87.}}, {5, {30., -82.}}
    };
    return map;
  }

  static int cherenkov_eo_channel(const data_t &data)
  {
    int chip = data.fifo / 4;
    int chip_side = chip % 2;

    // Decoded 2026 data keep the ALCOR channel convention used elsewhere:
    // eoch = pixel + 4 * column.  The matrix table addresses the two chips
    // on one ALCOR-dual board as eoch=0..63.
    int local_eoch = data.pixel + 4 * data.column;
    int eoch = local_eoch + 32 * chip_side;
    if (local_eoch < 0 || local_eoch >= 32 || eoch < 0 || eoch >= 64)
      return -1;
    return eoch;
  }

  static bool cherenkov_coordinate(const data_t &data, double &x, double &y)
  {
    int chip = data.fifo / 4;
    int eoch = cherenkov_eo_channel(data);
    if (chip < 0 || chip > 7 || eoch < 0)
      return false;

    auto pm = pdu_matrix_map().find({data.device, chip});
    if (pm == pdu_matrix_map().end())
      return false;

    int pdu = pm->second[0];
    int matrix = pm->second[1];
    auto mm = matrix_mapping().find(matrix);
    if (mm == matrix_mapping().end() || eoch >= (int)mm->second.size())
      return false;

    int doch = mm->second[eoch];
    int col = doch / 8;
    int row = doch % 8;

    if (matrix == 2 || matrix == 4)
      row += 8;
    if (matrix == 3 || matrix == 4)
      col += 8;

    col = 15 - col;
    row = 15 - row;

    double xx = 0.05 + 0.1 + 0.2 + 1.5 + 3.2 * col;
    double yy = 0.05 + 0.1 + 0.2 + 1.5 + 3.2 * row;
    if (col > 7)
      xx += 0.3;
    if (row > 7)
      yy += 0.3;

    auto place = placement_xy().find(pdu);
    if (place == placement_xy().end())
      return false;

    x = xx + place->second[0];
    y = yy + place->second[1];
    return true;
  }
};
