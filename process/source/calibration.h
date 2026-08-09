#pragma once

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct tdc_calib_t {
  double off = 0.;
  double iif = 0.;
};

struct channel_calib_t {
  double offset = 0.;
};

struct trigger_calib_t {
  double offset = 0.;
};

class calibration_t {
public:
  bool load(const std::string &filename);

  bool tdc(int device, int fifo, int column, int pixel, int tdc,
           tdc_calib_t &out, std::string &error) const;
  bool channel(int device, int fifo, int column, int pixel,
               channel_calib_t &out, std::string &error) const;
  bool trigger(int device, int fifo, trigger_calib_t &out,
               std::string &error) const;

  int tdc_rules() const { return tdc_rules_.size(); }
  int channel_rules() const { return channel_rules_.size(); }
  int trigger_rules() const { return trigger_rules_.size(); }
  int cached_tdc() const { return tdc_cache_.size(); }
  int cached_channels() const { return channel_cache_.size(); }
  int cached_triggers() const { return trigger_cache_.size(); }

private:
  struct field_t {
    bool wildcard = true;
    int value = 0;
    bool match(int x) const { return wildcard || value == x; }
    int specificity() const { return wildcard ? 0 : 1; }
  };

  struct tdc_rule_t {
    std::array<field_t, 5> key;
    tdc_calib_t calib;
    int line = 0;
  };

  struct channel_rule_t {
    std::array<field_t, 4> key;
    channel_calib_t calib;
    int line = 0;
  };

  struct trigger_rule_t {
    std::array<field_t, 2> key;
    trigger_calib_t calib;
    int line = 0;
  };

  struct tdc_key_t {
    int device;
    int fifo;
    int column;
    int pixel;
    int tdc;
    bool operator==(const tdc_key_t &rhs) const
    {
      return device == rhs.device && fifo == rhs.fifo && column == rhs.column &&
             pixel == rhs.pixel && tdc == rhs.tdc;
    }
  };

  struct channel_key_t {
    int device;
    int fifo;
    int column;
    int pixel;
    bool operator==(const channel_key_t &rhs) const
    {
      return device == rhs.device && fifo == rhs.fifo && column == rhs.column && pixel == rhs.pixel;
    }
  };

  struct trigger_key_t {
    int device;
    int fifo;
    bool operator==(const trigger_key_t &rhs) const
    {
      return device == rhs.device && fifo == rhs.fifo;
    }
  };

  struct tdc_hash_t {
    std::size_t operator()(const tdc_key_t &k) const
    {
      std::size_t h = 1469598103934665603ull;
      auto mix = [&](int v) { h = (h ^ std::hash<int>{}(v)) * 1099511628211ull; };
      mix(k.device); mix(k.fifo); mix(k.column); mix(k.pixel); mix(k.tdc);
      return h;
    }
  };

  struct channel_hash_t {
    std::size_t operator()(const channel_key_t &k) const
    {
      std::size_t h = 1469598103934665603ull;
      auto mix = [&](int v) { h = (h ^ std::hash<int>{}(v)) * 1099511628211ull; };
      mix(k.device); mix(k.fifo); mix(k.column); mix(k.pixel);
      return h;
    }
  };

  struct trigger_hash_t {
    std::size_t operator()(const trigger_key_t &k) const
    {
      std::size_t h = 1469598103934665603ull;
      h = (h ^ std::hash<int>{}(k.device)) * 1099511628211ull;
      h = (h ^ std::hash<int>{}(k.fifo)) * 1099511628211ull;
      return h;
    }
  };

  static field_t parse_field(const std::string &text, int line, const std::string &name);
  static int parse_int(const std::string &text, int line, const std::string &name);
  static double parse_double(const std::string &text, int line, const std::string &name);
  static std::string strip_comment(const std::string &line);
  static bool same_rule(const std::array<field_t, 5> &a, const std::array<field_t, 5> &b);
  static bool same_rule(const std::array<field_t, 4> &a, const std::array<field_t, 4> &b);
  static bool same_rule(const std::array<field_t, 2> &a, const std::array<field_t, 2> &b);
  static std::string tdc_label(const tdc_key_t &key);
  static std::string channel_label(const channel_key_t &key);
  static std::string trigger_label(const trigger_key_t &key);

  std::vector<tdc_rule_t> tdc_rules_;
  std::vector<channel_rule_t> channel_rules_;
  std::vector<trigger_rule_t> trigger_rules_;

  mutable std::unordered_map<tdc_key_t, tdc_calib_t, tdc_hash_t> tdc_cache_;
  mutable std::unordered_map<channel_key_t, channel_calib_t, channel_hash_t> channel_cache_;
  mutable std::unordered_map<trigger_key_t, trigger_calib_t, trigger_hash_t> trigger_cache_;
};

inline std::string
calibration_t::strip_comment(const std::string &line)
{
  auto pos = line.find('#');
  return pos == std::string::npos ? line : line.substr(0, pos);
}

inline int
calibration_t::parse_int(const std::string &text, int line, const std::string &name)
{
  std::size_t pos = 0;
  try {
    int value = std::stoi(text, &pos);
    if (pos != text.size())
      throw std::invalid_argument("trailing characters");
    return value;
  } catch (...) {
    throw std::runtime_error("config line " + std::to_string(line) + ": malformed integer for " + name + ": '" + text + "'");
  }
}

inline double
calibration_t::parse_double(const std::string &text, int line, const std::string &name)
{
  if (text == "*")
    throw std::runtime_error("config line " + std::to_string(line) + ": wildcard is not valid for calibration value " + name);

  std::size_t pos = 0;
  try {
    double value = std::stod(text, &pos);
    if (pos != text.size())
      throw std::invalid_argument("trailing characters");
    return value;
  } catch (...) {
    throw std::runtime_error("config line " + std::to_string(line) + ": malformed floating-point value for " + name + ": '" + text + "'");
  }
}

inline calibration_t::field_t
calibration_t::parse_field(const std::string &text, int line, const std::string &name)
{
  field_t out;
  if (text == "*")
    return out;
  out.wildcard = false;
  out.value = parse_int(text, line, name);
  return out;
}

inline bool
calibration_t::same_rule(const std::array<field_t, 5> &a, const std::array<field_t, 5> &b)
{
  for (int i = 0; i < 5; ++i)
    if (a[i].wildcard != b[i].wildcard || (!a[i].wildcard && a[i].value != b[i].value))
      return false;
  return true;
}

inline bool
calibration_t::same_rule(const std::array<field_t, 4> &a, const std::array<field_t, 4> &b)
{
  for (int i = 0; i < 4; ++i)
    if (a[i].wildcard != b[i].wildcard || (!a[i].wildcard && a[i].value != b[i].value))
      return false;
  return true;
}

inline bool
calibration_t::same_rule(const std::array<field_t, 2> &a, const std::array<field_t, 2> &b)
{
  for (int i = 0; i < 2; ++i)
    if (a[i].wildcard != b[i].wildcard || (!a[i].wildcard && a[i].value != b[i].value))
      return false;
  return true;
}

inline std::string
calibration_t::tdc_label(const tdc_key_t &key)
{
  std::ostringstream out;
  out << "device=" << key.device << " fifo=" << key.fifo
      << " column=" << key.column << " pixel=" << key.pixel
      << " tdc=" << key.tdc;
  return out.str();
}

inline std::string
calibration_t::channel_label(const channel_key_t &key)
{
  std::ostringstream out;
  out << "device=" << key.device << " fifo=" << key.fifo
      << " column=" << key.column << " pixel=" << key.pixel;
  return out.str();
}

inline std::string
calibration_t::trigger_label(const trigger_key_t &key)
{
  std::ostringstream out;
  out << "device=" << key.device << " fifo=" << key.fifo;
  return out.str();
}

inline bool
calibration_t::load(const std::string &filename)
{
  tdc_rules_.clear();
  channel_rules_.clear();
  trigger_rules_.clear();
  tdc_cache_.clear();
  channel_cache_.clear();
  trigger_cache_.clear();

  std::ifstream in(filename);
  if (!in) {
    std::cerr << "ERROR: could not open calibration file: " << filename << std::endl;
    return false;
  }

  enum class section_t { none, tdc, channel, trigger } section = section_t::none;
  std::string line;
  int lineno = 0;

  try {
    while (std::getline(in, line)) {
      ++lineno;
      auto clean = strip_comment(line);
      std::istringstream ss(clean);
      std::vector<std::string> tok;
      std::string word;
      while (ss >> word)
        tok.push_back(word);
      if (tok.empty())
        continue;

      if (tok[0].size() >= 2 && tok[0].front() == '[' && tok[0].back() == ']') {
        if (tok.size() != 1)
          throw std::runtime_error("config line " + std::to_string(lineno) + ": section header must be alone on the line");
        if (tok[0] == "[TDC]") section = section_t::tdc;
        else if (tok[0] == "[CHANNEL]") section = section_t::channel;
        else if (tok[0] == "[TRIGGER]") section = section_t::trigger;
        else throw std::runtime_error("config line " + std::to_string(lineno) + ": unknown section '" + tok[0] + "'");
        continue;
      }

      if (section == section_t::none)
        throw std::runtime_error("config line " + std::to_string(lineno) + ": calibration row before any section");

      if (section == section_t::tdc) {
        if (tok.size() != 7)
          throw std::runtime_error("config line " + std::to_string(lineno) + ": [TDC] expects 7 columns");
        tdc_rule_t rule;
        for (int i = 0; i < 5; ++i)
          rule.key[i] = parse_field(tok[i], lineno, "TDC address");
        if (!rule.key[4].wildcard && (rule.key[4].value < 0 || rule.key[4].value >= 4))
          throw std::runtime_error("config line " + std::to_string(lineno) + ": invalid TDC index " + std::to_string(rule.key[4].value));
        rule.calib.off = parse_double(tok[5], lineno, "off");
        rule.calib.iif = parse_double(tok[6], lineno, "iif");
        rule.line = lineno;
        for (auto &old : tdc_rules_) {
          if (same_rule(old.key, rule.key))
            throw std::runtime_error("config line " + std::to_string(lineno) + ": duplicate [TDC] rule also defined on line " + std::to_string(old.line));
        }
        tdc_rules_.push_back(rule);
      } else if (section == section_t::channel) {
        if (tok.size() != 5)
          throw std::runtime_error("config line " + std::to_string(lineno) + ": [CHANNEL] expects 5 columns");
        channel_rule_t rule;
        for (int i = 0; i < 4; ++i)
          rule.key[i] = parse_field(tok[i], lineno, "CHANNEL address");
        rule.calib.offset = parse_double(tok[4], lineno, "offset");
        rule.line = lineno;
        for (auto &old : channel_rules_) {
          if (same_rule(old.key, rule.key))
            throw std::runtime_error("config line " + std::to_string(lineno) + ": duplicate [CHANNEL] rule also defined on line " + std::to_string(old.line));
        }
        channel_rules_.push_back(rule);
      } else if (section == section_t::trigger) {
        if (tok.size() != 3)
          throw std::runtime_error("config line " + std::to_string(lineno) + ": [TRIGGER] expects 3 columns");
        trigger_rule_t rule;
        for (int i = 0; i < 2; ++i)
          rule.key[i] = parse_field(tok[i], lineno, "TRIGGER address");
        rule.calib.offset = parse_double(tok[2], lineno, "offset");
        rule.line = lineno;
        for (auto &old : trigger_rules_) {
          if (same_rule(old.key, rule.key))
            throw std::runtime_error("config line " + std::to_string(lineno) + ": duplicate [TRIGGER] rule also defined on line " + std::to_string(old.line));
        }
        trigger_rules_.push_back(rule);
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return false;
  }

  return true;
}

inline bool
calibration_t::tdc(int device, int fifo, int column, int pixel, int itdc,
                   tdc_calib_t &out, std::string &error) const
{
  tdc_key_t key{device, fifo, column, pixel, itdc};
  auto cached = tdc_cache_.find(key);
  if (cached != tdc_cache_.end()) {
    out = cached->second;
    return true;
  }

  int best_spec = -1;
  const tdc_rule_t *best = nullptr;
  for (auto &rule : tdc_rules_) {
    if (!rule.key[0].match(device) || !rule.key[1].match(fifo) || !rule.key[2].match(column) ||
        !rule.key[3].match(pixel) || !rule.key[4].match(itdc))
      continue;
    int spec = rule.key[0].specificity() + rule.key[1].specificity() + rule.key[2].specificity() +
               rule.key[3].specificity() + rule.key[4].specificity();
    if (spec > best_spec) {
      best_spec = spec;
      best = &rule;
    } else if (spec == best_spec) {
      error = "ambiguous [TDC] calibration for " + tdc_label(key) +
              ": lines " + std::to_string(best->line) + " and " + std::to_string(rule.line);
      return false;
    }
  }

  if (!best) {
    error = "missing [TDC] calibration for " + tdc_label(key);
    return false;
  }

  out = best->calib;
  tdc_cache_[key] = out;
  return true;
}

inline bool
calibration_t::channel(int device, int fifo, int column, int pixel,
                       channel_calib_t &out, std::string &error) const
{
  channel_key_t key{device, fifo, column, pixel};
  auto cached = channel_cache_.find(key);
  if (cached != channel_cache_.end()) {
    out = cached->second;
    return true;
  }

  int best_spec = -1;
  const channel_rule_t *best = nullptr;
  for (auto &rule : channel_rules_) {
    if (!rule.key[0].match(device) || !rule.key[1].match(fifo) || !rule.key[2].match(column) || !rule.key[3].match(pixel))
      continue;
    int spec = rule.key[0].specificity() + rule.key[1].specificity() + rule.key[2].specificity() + rule.key[3].specificity();
    if (spec > best_spec) {
      best_spec = spec;
      best = &rule;
    } else if (spec == best_spec) {
      error = "ambiguous [CHANNEL] calibration for " + channel_label(key) +
              ": lines " + std::to_string(best->line) + " and " + std::to_string(rule.line);
      return false;
    }
  }

  if (!best) {
    error = "missing [CHANNEL] calibration for " + channel_label(key);
    return false;
  }

  out = best->calib;
  channel_cache_[key] = out;
  return true;
}

inline bool
calibration_t::trigger(int device, int fifo, trigger_calib_t &out,
                       std::string &error) const
{
  trigger_key_t key{device, fifo};
  auto cached = trigger_cache_.find(key);
  if (cached != trigger_cache_.end()) {
    out = cached->second;
    return true;
  }

  int best_spec = -1;
  const trigger_rule_t *best = nullptr;
  for (auto &rule : trigger_rules_) {
    if (!rule.key[0].match(device) || !rule.key[1].match(fifo))
      continue;
    int spec = rule.key[0].specificity() + rule.key[1].specificity();
    if (spec > best_spec) {
      best_spec = spec;
      best = &rule;
    } else if (spec == best_spec) {
      error = "ambiguous [TRIGGER] calibration for " + trigger_label(key) +
              ": lines " + std::to_string(best->line) + " and " + std::to_string(rule.line);
      return false;
    }
  }

  if (!best) {
    error = "missing [TRIGGER] calibration for " + trigger_label(key);
    return false;
  }

  out = best->calib;
  trigger_cache_[key] = out;
  return true;
}
