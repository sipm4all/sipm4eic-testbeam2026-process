#pragma once

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class clock_calibration_t {
public:
  bool load(const std::string &filename, const std::string &run);
  int correction(const std::string &run, int device, int fifo, int spill) const;
  int rules() const { return rules_; }

private:
  struct entry_t {
    std::string run;
    int device;
    int fifo;
    int correction;
    std::set<int> spills;
  };

  std::vector<entry_t> entries_;
  int rules_ = 0;
};

inline bool
clock_calibration_t::load(const std::string &filename, const std::string &run)
{
  entries_.clear();
  rules_ = 0;
  std::ifstream input(filename);
  if (!input) {
    std::cerr << "ERROR: could not open clock correction file: " << filename << std::endl;
    return false;
  }

  std::string line;
  int lineno = 0;
  bool in_clock_section = false;
  std::set<std::tuple<std::string, int, int>> seen;
  try {
    while (std::getline(input, line)) {
      ++lineno;
      auto comment = line.find('#');
      if (comment != std::string::npos)
        line.resize(comment);
      std::istringstream stream(line);
      std::vector<std::string> token;
      std::string word;
      while (stream >> word)
        token.push_back(word);
      if (token.empty())
        continue;
      if (token[0] == "[CLOCK]") {
        if (token.size() != 1)
          throw std::runtime_error("line " + std::to_string(lineno) + ": malformed [CLOCK] header");
        in_clock_section = true;
        continue;
      }
      if (!in_clock_section)
        throw std::runtime_error("line " + std::to_string(lineno) + ": clock row before [CLOCK] section");
      if (token.size() < 5)
        throw std::runtime_error("line " + std::to_string(lineno) + ": [CLOCK] expects run device fifo correction spill...");
      std::size_t pos = 0;
      int device = std::stoi(token[1], &pos);
      if (pos != token[1].size())
        throw std::runtime_error("line " + std::to_string(lineno) + ": malformed device");
      int fifo = std::stoi(token[2], &pos);
      if (pos != token[2].size() || fifo < 0)
        throw std::runtime_error("line " + std::to_string(lineno) + ": malformed fifo");
      int sign = std::stoi(token[3], &pos);
      if (pos != token[3].size() || (sign != -1 && sign != 1))
        throw std::runtime_error("line " + std::to_string(lineno) + ": correction must be +1 or -1");

      entry_t entry{token[0], device, fifo, sign, {}};
      for (std::size_t i = 4; i < token.size(); ++i) {
        int spill = std::stoi(token[i], &pos);
        if (pos != token[i].size() || spill < 0)
          throw std::runtime_error("line " + std::to_string(lineno) + ": malformed spill number");
        if (!entry.spills.insert(spill).second)
          throw std::runtime_error("line " + std::to_string(lineno) + ": duplicate spill number");
      }
      if (entry.spills.empty())
        throw std::runtime_error("line " + std::to_string(lineno) + ": empty spill list");

      // A combined clock file may contain rows for many runs. Only retain
      // rows for the run selected by calibrator's --run argument.
      if (token[0] != run)
        continue;

      auto key = std::make_tuple(run, device, fifo);
      if (!seen.insert(key).second)
        throw std::runtime_error("line " + std::to_string(lineno) + ": duplicate correction for run/device/fifo");
      entry.run = run;
      entries_.push_back(std::move(entry));
      ++rules_;
    }
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    entries_.clear();
    rules_ = 0;
    return false;
  }
  return true;
}

inline int
clock_calibration_t::correction(const std::string &run, int device, int fifo, int spill) const
{
  for (const auto &entry : entries_)
    if (entry.run == run && entry.device == device && entry.fifo == fifo &&
        entry.spills.count(spill) != 0)
      return entry.correction;
  return 0;
}
