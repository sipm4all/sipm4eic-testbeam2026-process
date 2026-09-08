#pragma once
#include <vector>
namespace irt_analysis_cuda {
struct peak_t { float theta, time, score; };
bool available();
bool peak(const std::vector<float> &theta, const std::vector<float> &time,
          peak_t &result);
}
