#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <string>
#include <tuple>
#include <vector>


constexpr int maxframes = 65536;
constexpr int maxhits = 1024;
constexpr int maxspillhits = maxframes * maxhits;
constexpr int maxsources = 4096;

struct spill_participation_t {
  int spill = 0;
  int counter = 0;
  int nsources = 0;
  int source_device[maxsources];
  int source_fifo[maxsources];

  void clear(int _spill, int _counter)
  {
    spill = _spill;
    counter = _counter;
    nsources = 0;
  }

  bool add(int device, int fifo)
  {
    if (nsources >= maxsources)
      return false;
    source_device[nsources] = device;
    source_fifo[nsources] = fifo;
    ++nsources;
    return true;
  }
};

struct field_match_t {
  enum class mode_t { wildcard, exact, range, set };

  mode_t mode = mode_t::wildcard;
  int exact_value = 0;
  int range_min = 0;
  int range_max = 0;
  std::vector<int> values;

  bool match(int value) const
  {
    if (mode == mode_t::wildcard) return true;
    if (mode == mode_t::exact) return value == exact_value;
    if (mode == mode_t::range) return value >= range_min && value <= range_max;
    return std::find(values.begin(), values.end(), value) != values.end();
  }

  std::string str() const
  {
    if (mode == mode_t::wildcard) return "*";
    if (mode == mode_t::exact) return std::to_string(exact_value);
    if (mode == mode_t::range) return "[" + std::to_string(range_min) + "," + std::to_string(range_max) + "]";

    std::string out = "{";
    for (int i = 0; i < (int)values.size(); ++i) {
      if (i > 0) out += ",";
      out += std::to_string(values[i]);
    }
    out += "}";
    return out;
  }
};

struct selector_t {
  field_match_t type;
  field_match_t device;
  field_match_t fifo;
  field_match_t column;
  field_match_t pixel;
  field_match_t counter;
  field_match_t tdc;
  double time_offset = 0.;

  bool match(const data_t &data) const
  {
    return type.match(data.type) &&
           device.match(data.device) &&
           fifo.match(data.fifo) &&
           column.match(data.column) &&
           pixel.match(data.pixel) &&
           counter.match(data.counter) &&
           tdc.match(data.tdc);
  }
};

struct condition_t {
  selector_t selector;
  double dtmin = 0.;
  double dtmax = 0.;
  int min_count = 1;
  bool veto = false;

  double raw_dtmin(double reference_time_offset) const
  {
    return dtmin + selector.time_offset - reference_time_offset;
  }

  double raw_dtmax(double reference_time_offset) const
  {
    return dtmax + selector.time_offset - reference_time_offset;
  }

  double latest_raw_time(double event_time) const
  {
    return event_time + selector.time_offset + dtmax;
  }
};

enum class trigger_mode_t { none, seeded, coincidence };

struct coincidence_t {
  selector_t selector;
  std::string selector_name;
  double dtmin = 0.;
  double dtmax = 0.;
  int min_count = 1;
  bool unique_channel = false;
};

struct trigger_definition_t {
  trigger_mode_t mode = trigger_mode_t::none;
  selector_t seed;
  std::string seed_name;
  bool has_seed = false;
  coincidence_t coincidence;
  bool has_coincidence = false;
  selector_t frame_center_selector;
  std::string frame_center_name;
  bool has_frame_center = false;
  std::vector<condition_t> conditions;

  double frame_event_offset() const
  {
    if (has_frame_center)
      return frame_center_selector.time_offset;
    return 0.;
  }

  double frame_raw_dt() const
  {
    if (mode == trigger_mode_t::seeded) {
      if (has_frame_center)
        return frame_center_selector.time_offset - seed.time_offset;
      return -seed.time_offset;
    }
    return frame_event_offset();
  }

  double seeded_decision_raw_dtmax() const
  {
    double out = 0.;
    for (auto &condition : conditions) {
      auto raw_dtmax = condition.raw_dtmax(seed.time_offset);
      if (raw_dtmax > out)
        out = raw_dtmax;
    }
    return out;
  }

  double decision_latest_raw_time(double event_time) const
  {
    double out = event_time;
    for (auto &condition : conditions) {
      auto latest = condition.latest_raw_time(event_time);
      if (latest > out)
        out = latest;
    }
    return out;
  }

  double history_window(double frame_window) const
  {
    double out = 0.;

    if (mode == trigger_mode_t::seeded) {
      auto frame_dt = frame_raw_dt();
      if (frame_window - frame_dt > out)
        out = frame_window - frame_dt;
      for (auto &condition : conditions) {
        auto raw_dtmin = condition.raw_dtmin(seed.time_offset);
        if (-raw_dtmin > out)
          out = -raw_dtmin;
      }
      return out;
    }

    if (mode == trigger_mode_t::coincidence) {
      auto c = coincidence.selector.time_offset;
      auto frame_offset = frame_event_offset();
      auto frame_history = frame_window + c + coincidence.dtmax - frame_offset;
      if (frame_history > out)
        out = frame_history;
      for (auto &condition : conditions) {
        auto condition_history = c + coincidence.dtmax - condition.selector.time_offset - condition.dtmin;
        if (condition_history > out)
          out = condition_history;
      }
    }

    return out;
  }
};

struct config_parser_t {
  struct token_t { std::string text; int line; };
  std::vector<token_t> tokens;
  int pos = 0;
  std::map<std::string, selector_t> selectors;
  trigger_definition_t trigger;
  bool have_trigger = false;

  [[noreturn]] void error(const std::string &message) const
  {
    int line = pos < (int)tokens.size() ? tokens[pos].line : (tokens.empty() ? 0 : tokens.back().line);
    throw std::runtime_error("config line " + std::to_string(line) + ": " + message);
  }

  bool end() const { return pos >= (int)tokens.size(); }
  std::string peek() const { return end() ? "" : tokens[pos].text; }
  std::string get()
  {
    if (end()) error("unexpected end of file");
    return tokens[pos++].text;
  }
  void expect(const std::string &value)
  {
    auto token = get();
    if (token != value)
      error("expected '" + value + "', found '" + token + "'");
  }

  static bool is_punct(char c)
  {
    return c == '{' || c == '}' || c == '=' || c == '[' || c == ']' || c == ',';
  }

  void tokenize(const std::string &filename)
  {
    std::ifstream in(filename);
    if (!in)
      throw std::runtime_error("could not open config file: " + filename);

    std::string line;
    int iline = 0;
    while (std::getline(in, line)) {
      ++iline;
      auto comment = line.find('#');
      if (comment != std::string::npos)
        line.resize(comment);

      std::string spaced;
      for (char c : line) {
        if (is_punct(c)) {
          spaced += ' ';
          spaced += c;
          spaced += ' ';
        } else {
          spaced += c;
        }
      }

      std::istringstream ss(spaced);
      std::string token;
      while (ss >> token)
        tokens.push_back({token, iline});
    }
  }

  int parse_int_token()
  {
    auto token = get();
    char *endp = nullptr;
    long value = std::strtol(token.c_str(), &endp, 10);
    if (!endp || *endp != 0)
      error("malformed integer value '" + token + "'");
    return (int)value;
  }

  double parse_double_token()
  {
    auto token = get();
    char *endp = nullptr;
    double value = std::strtod(token.c_str(), &endp);
    if (!endp || *endp != 0)
      error("malformed numeric value '" + token + "'");
    return value;
  }

  field_match_t parse_field_match()
  {
    field_match_t match;

    if (peek() == "*") {
      get();
      match.mode = field_match_t::mode_t::wildcard;
      return match;
    }

    if (peek() == "[") {
      get();
      match.mode = field_match_t::mode_t::range;
      match.range_min = parse_int_token();
      expect(",");
      match.range_max = parse_int_token();
      expect("]");
      if (match.range_min > match.range_max)
        error("range minimum is greater than range maximum");
      return match;
    }

    if (peek() == "{") {
      get();
      match.mode = field_match_t::mode_t::set;
      if (peek() == "}")
        error("empty set is not allowed");

      while (true) {
        match.values.push_back(parse_int_token());
        if (peek() == "}")
          break;
        expect(",");
        if (peek() == "}")
          error("trailing comma in set");
      }
      expect("}");
      return match;
    }

    match.mode = field_match_t::mode_t::exact;
    match.exact_value = parse_int_token();
    return match;
  }

  void set_selector_field(selector_t &selector, const std::string &field)
  {
    expect("=");

    if (field == "time_offset") { selector.time_offset = parse_double_token(); return; }

    if (field != "type" && field != "device" && field != "fifo" &&
        field != "column" && field != "pixel" && field != "counter" &&
        field != "tdc")
      error("unknown selector field '" + field + "'");

    auto match = parse_field_match();
    if (field == "type") { selector.type = match; return; }
    if (field == "device") { selector.device = match; return; }
    if (field == "fifo") { selector.fifo = match; return; }
    if (field == "column") { selector.column = match; return; }
    if (field == "pixel") { selector.pixel = match; return; }
    if (field == "counter") { selector.counter = match; return; }
    if (field == "tdc") { selector.tdc = match; return; }
  }

  void parse_selector()
  {
    auto name = get();
    if (selectors.count(name))
      error("duplicate selector '" + name + "'");
    expect("{");

    selector_t selector;
    while (peek() != "}") {
      if (end()) error("missing '}' in selector");
      auto field = get();
      set_selector_field(selector, field);
    }
    expect("}");
    selectors[name] = selector;
  }

  void parse_dt(condition_t &condition)
  {
    parse_dt(condition.dtmin, condition.dtmax);
  }

  void parse_dt(double &dtmin, double &dtmax)
  {
    expect("=");
    expect("[");
    dtmin = parse_double_token();
    expect(",");
    dtmax = parse_double_token();
    expect("]");
    if (dtmin > dtmax)
      error("dtmin > dtmax");
  }

  coincidence_t parse_coincidence()
  {
    auto selector_name = get();
    auto found = selectors.find(selector_name);
    if (found == selectors.end())
      error("unknown coincidence selector '" + selector_name + "'");

    coincidence_t coincidence;
    coincidence.selector = found->second;
    coincidence.selector_name = selector_name;

    expect("{");
    bool have_dt = false;
    bool have_min = false;
    while (peek() != "}") {
      if (end()) error("missing '}' in coincidence block");
      auto key = get();
      if (key == "dt") {
        parse_dt(coincidence.dtmin, coincidence.dtmax);
        have_dt = true;
      } else if (key == "min") {
        expect("=");
        coincidence.min_count = parse_int_token();
        if (coincidence.min_count <= 0)
          error("coincidence min must be > 0");
        have_min = true;
      } else if (key == "unique") {
        expect("=");
        auto value = get();
        if (value != "channel")
          error("unknown unique value '" + value + "'");
        coincidence.unique_channel = true;
      } else {
        error("unknown coincidence keyword '" + key + "'");
      }
    }
    expect("}");

    if (!have_dt)
      error("coincidence is missing dt interval");
    if (!have_min)
      error("coincidence is missing min");
    if (coincidence.dtmin > 0. || coincidence.dtmax < 0.)
      error("coincidence dt interval must contain 0 for earliest-hit reference convention");
    return coincidence;
  }

  condition_t parse_condition(bool veto)
  {
    auto selector_name = get();
    auto found = selectors.find(selector_name);
    if (found == selectors.end())
      error("unknown selector '" + selector_name + "'");

    condition_t condition;
    condition.selector = found->second;
    condition.veto = veto;
    if (veto)
      condition.min_count = 0;

    expect("{");
    bool have_dt = false;
    while (peek() != "}") {
      if (end()) error("missing '}' in condition");
      auto key = get();
      if (key == "dt") {
        parse_dt(condition);
        have_dt = true;
      } else if (key == "min" && !veto) {
        expect("=");
        condition.min_count = parse_int_token();
        if (condition.min_count < 0)
          error("negative min");
      } else {
        error("unknown condition keyword '" + key + "'");
      }
    }
    expect("}");

    if (!have_dt)
      error("condition is missing dt interval");
    return condition;
  }

  void parse_trigger()
  {
    if (have_trigger)
      error("multiple trigger blocks are not supported");
    have_trigger = true;
    expect("{");

    while (peek() != "}") {
      if (end()) error("missing '}' in trigger block");
      auto key = get();
      if (key == "seed") {
        if (trigger.has_coincidence)
          error("trigger block cannot contain both seed and coincidence");
        if (trigger.has_seed)
          error("duplicate trigger seed");
        auto selector_name = get();
        auto found = selectors.find(selector_name);
        if (found == selectors.end())
          error("unknown seed selector '" + selector_name + "'");
        trigger.seed = found->second;
        trigger.seed_name = selector_name;
        trigger.has_seed = true;
        trigger.mode = trigger_mode_t::seeded;
      } else if (key == "coincidence") {
        if (trigger.has_seed)
          error("trigger block cannot contain both seed and coincidence");
        if (trigger.has_coincidence)
          error("multiple coincidence blocks are not supported");
        trigger.coincidence = parse_coincidence();
        trigger.has_coincidence = true;
        trigger.mode = trigger_mode_t::coincidence;
      } else if (key == "frame_center") {
        if (trigger.has_frame_center)
          error("duplicate trigger frame_center");
        auto selector_name = get();
        auto found = selectors.find(selector_name);
        if (found == selectors.end())
          error("unknown frame_center selector '" + selector_name + "'");
        trigger.frame_center_selector = found->second;
        trigger.frame_center_name = selector_name;
        trigger.has_frame_center = true;
      } else if (key == "time_offset") {
        error("trigger-level time_offset has been removed; use selector time_offset and optional frame_center instead");
      } else if (key == "require") {
        trigger.conditions.push_back(parse_condition(false));
      } else if (key == "veto") {
        trigger.conditions.push_back(parse_condition(true));
      } else {
        error("unknown trigger keyword '" + key + "'");
      }
    }
    expect("}");
  }

  trigger_definition_t parse(const std::string &filename)
  {
    tokenize(filename);
    while (!end()) {
      auto key = get();
      if (key == "selector") parse_selector();
      else if (key == "trigger") parse_trigger();
      else error("unknown top-level keyword '" + key + "'");
    }

    if (!have_trigger)
      error("missing trigger block");
    if (trigger.mode == trigger_mode_t::none)
      error("trigger block is missing seed or coincidence");
    return trigger;
  }
};

struct frame_t {
  int nhits;

  int device[maxhits];
  int fifo[maxhits];
  int type[maxhits];
  int counter[maxhits];
  int column[maxhits];
  int pixel[maxhits];
  int tdc[maxhits];
  int rollover[maxhits];
  int coarse[maxhits];
  int fine[maxhits];
  double time[maxhits];

  bool add(const data_t &data)
  {
    if (nhits >= maxhits)
      return false;

    device[nhits] = data.device;
    fifo[nhits] = data.fifo;
    type[nhits] = data.type;
    counter[nhits] = data.counter;
    column[nhits] = data.column;
    pixel[nhits] = data.pixel;
    tdc[nhits] = data.tdc;
    rollover[nhits] = data.rollover;
    coarse[nhits] = data.coarse;
    fine[nhits] = data.fine;
    time[nhits] = data.time;
    ++nhits;
    return true;
  }
};

struct hit_store_t {
  int nhits;
  std::unique_ptr<int[]> frame_start;
  std::unique_ptr<int[]> frame_nhits;
  std::unique_ptr<int[]> device;
  std::unique_ptr<int[]> fifo;
  std::unique_ptr<int[]> type;
  std::unique_ptr<int[]> counter;
  std::unique_ptr<int[]> column;
  std::unique_ptr<int[]> pixel;
  std::unique_ptr<int[]> tdc;
  std::unique_ptr<int[]> rollover;
  std::unique_ptr<int[]> coarse;
  std::unique_ptr<int[]> fine;
  std::unique_ptr<double[]> time;

  hit_store_t()
    : frame_start(new int[maxframes]),
      frame_nhits(new int[maxframes]),
      device(new int[maxspillhits]),
      fifo(new int[maxspillhits]),
      type(new int[maxspillhits]),
      counter(new int[maxspillhits]),
      column(new int[maxspillhits]),
      pixel(new int[maxspillhits]),
      tdc(new int[maxspillhits]),
      rollover(new int[maxspillhits]),
      coarse(new int[maxspillhits]),
      fine(new int[maxspillhits]),
      time(new double[maxspillhits])
  {
    reset();
  }

  void reset()
  {
    nhits = 0;
  }

  void start_frame(int iframe)
  {
    frame_start[iframe] = nhits;
    frame_nhits[iframe] = 0;
  }

  bool add(const frame_t &frame, int ihit, int iframe)
  {
    if (nhits >= maxspillhits)
      return false;

    int j = nhits++;
    device[j] = frame.device[ihit];
    fifo[j] = frame.fifo[ihit];
    type[j] = frame.type[ihit];
    counter[j] = frame.counter[ihit];
    column[j] = frame.column[ihit];
    pixel[j] = frame.pixel[ihit];
    tdc[j] = frame.tdc[ihit];
    rollover[j] = frame.rollover[ihit];
    coarse[j] = frame.coarse[ihit];
    fine[j] = frame.fine[ihit];
    time[j] = frame.time[ihit];
    ++frame_nhits[iframe];
    return true;
  }
};

struct spill_t {
  int id;
  int nframes;
  hit_store_t trigger;
  hit_store_t timing;
  hit_store_t cherenkov;
  long long nunexpected_words;

  spill_t()
  {
    reset(0);
  }

  void reset(int _id)
  {
    id = _id;
    nframes = 0;
    trigger.reset();
    timing.reset();
    cherenkov.reset();
    nunexpected_words = 0;
  }

  bool add(const frame_t &frame)
  {
    if (nframes >= maxframes)
      return false;

    int ntrigger = 0;
    int ntiming = 0;
    int ncherenkov = 0;
    for (int i = 0; i < frame.nhits; ++i) {
      if (frame.type[i] == 9) ++ntrigger;
      else if (frame.type[i] == 1 && frame.device[i] == 200) ++ntiming;
      else if (frame.type[i] == 1) ++ncherenkov;
    }

    if (trigger.nhits + ntrigger > maxspillhits)
      return false;
    if (timing.nhits + ntiming > maxspillhits)
      return false;
    if (cherenkov.nhits + ncherenkov > maxspillhits)
      return false;

    int iframe = nframes;
    trigger.start_frame(iframe);
    timing.start_frame(iframe);
    cherenkov.start_frame(iframe);

    for (int i = 0; i < frame.nhits; ++i) {
      if (frame.type[i] == 9) {
        trigger.add(frame, i, iframe);
      } else if (frame.type[i] == 1 && frame.device[i] == 200) {
        timing.add(frame, i, iframe);
      } else if (frame.type[i] == 1) {
        cherenkov.add(frame, i, iframe);
      } else {
        ++nunexpected_words;
        std::cerr << "WARNING: unexpected word in accepted frame"
                  << " spill=" << id
                  << " frame=" << iframe
                  << " type=" << frame.type[i]
                  << " device=" << frame.device[i]
                  << " fifo=" << frame.fifo[i]
                  << " column=" << frame.column[i]
                  << " pixel=" << frame.pixel[i]
                  << " -- discarded" << std::endl;
      }
    }

    ++nframes;
    return true;
  }
};

struct active_frame_t {
  double event_time;
  double frame_center;
  double decision_time;
  int seed_entry;
  frame_t frame;
  std::vector<int> counts;
  bool decided = false;
  bool accepted = false;
};

struct coincidence_cluster_t {
  bool active = false;
  bool fired = false;
  double reference_time = 0.;
  int nhits = 0;
  std::set<std::tuple<int, int, int, int>> channels;

  void reset(double reference)
  {
    active = true;
    fired = false;
    reference_time = reference;
    nhits = 0;
    channels.clear();
  }

  void clear()
  {
    active = false;
    fired = false;
    nhits = 0;
    channels.clear();
  }

  void add(const data_t &data)
  {
    ++nhits;
    channels.insert(std::make_tuple(data.device, data.fifo, data.column, data.pixel));
  }

  int count(bool unique_channel) const
  {
    return unique_channel ? (int)channels.size() : nhits;
  }
};

static void
update_conditions(active_frame_t &active, const data_t &hit, const trigger_definition_t &definition)
{
  for (int i = 0; i < (int)definition.conditions.size(); ++i) {
    auto &condition = definition.conditions[i];
    auto delta_t = (hit.time - condition.selector.time_offset) - active.event_time;
    if (delta_t < condition.dtmin || delta_t > condition.dtmax)
      continue;
    if (condition.selector.match(hit))
      ++active.counts[i];
  }
}

static void
print_selector(const std::string &name, const selector_t &selector)
{
  std::cout << "selector " << name << ":" << std::endl;
  std::cout << "  type        = " << selector.type.str() << std::endl;
  std::cout << "  device      = " << selector.device.str() << std::endl;
  std::cout << "  fifo        = " << selector.fifo.str() << std::endl;
  std::cout << "  column      = " << selector.column.str() << std::endl;
  std::cout << "  pixel       = " << selector.pixel.str() << std::endl;
  std::cout << "  counter     = " << selector.counter.str() << std::endl;
  std::cout << "  tdc         = " << selector.tdc.str() << std::endl;
  std::cout << "  time_offset = " << selector.time_offset << std::endl;
}

static void
print_trigger(const trigger_definition_t &definition)
{
  std::cout << "trigger:" << std::endl;
  if (definition.mode == trigger_mode_t::seeded) {
    std::cout << "  mode         = seeded" << std::endl;
    std::cout << "  seed         = " << definition.seed_name << std::endl;
    if (definition.has_frame_center)
      std::cout << "  frame_center = " << definition.frame_center_name << std::endl;
    else
      std::cout << "  frame_center = event_time" << std::endl;
    std::cout << "  seed relative offset       = " << definition.seed.time_offset << std::endl;
    std::cout << "  frame relative to raw seed = " << definition.frame_raw_dt() << std::endl;
  } else {
    std::cout << "  mode         = coincidence" << std::endl;
    std::cout << "  coincidence  = " << definition.coincidence.selector_name << std::endl;
    std::cout << "  min          = " << definition.coincidence.min_count << std::endl;
    std::cout << "  unique       = " << (definition.coincidence.unique_channel ? "channel" : "hits") << std::endl;
    if (definition.has_frame_center)
      std::cout << "  frame_center = " << definition.frame_center_name << std::endl;
    else
      std::cout << "  frame_center = event_time" << std::endl;
  }
}

static bool
decide(active_frame_t &active, const trigger_definition_t &definition)
{
  for (int i = 0; i < (int)definition.conditions.size(); ++i) {
    auto &condition = definition.conditions[i];
    if (condition.veto) {
      if (active.counts[i] > 0)
        return false;
    } else {
      if (active.counts[i] < condition.min_count)
        return false;
    }
  }
  return true;
}

bool
trigger(const std::string filename,
        const std::string outfilename,
        const std::string configfilename,
        double window = 256.)
{
  trigger_definition_t definition;
  try {
    config_parser_t parser;
    definition = parser.parse(configfilename);
    std::cout << " --- trigger configuration: " << configfilename << std::endl;
    for (auto &selector : parser.selectors)
      print_selector(selector.first, selector.second);
    print_trigger(definition);
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return false;
  }

  auto fin = TFile::Open(filename.c_str());
  if (!fin || fin->IsZombie()) {
    std::cerr << " --- could not open input file: " << filename << std::endl;
    return false;
  }

  auto tin = (TTree *)fin->Get("alcor");
  if (!tin) {
    std::cerr << " --- could not find 'alcor' tree in input file" << std::endl;
    fin->Close();
    return false;
  }
  auto nev = tin->GetEntries();

  auto tmeta_in = (TTree *)fin->Get("spill_participation");
  spill_participation_t meta_in;
  Long64_t meta_entries = 0;
  Long64_t meta_entry = 0;
  if (tmeta_in) {
    meta_entries = tmeta_in->GetEntries();
    tmeta_in->SetBranchAddress("spill", &meta_in.spill);
    tmeta_in->SetBranchAddress("counter", &meta_in.counter);
    tmeta_in->SetBranchAddress("nsources", &meta_in.nsources);
    tmeta_in->SetBranchAddress("source_device", meta_in.source_device);
    tmeta_in->SetBranchAddress("source_fifo", meta_in.source_fifo);
  }

  data_t data;
  data.link_to_tree(tin);

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    fin->Close();
    return false;
  }

  static spill_t spill;
  auto tout = new TTree("frames", "triggered frames");
  tout->Branch("id", &spill.id, "id/I");
  tout->Branch("nframes", &spill.nframes, "nframes/I");

  auto branch_store = [&](const std::string &prefix, hit_store_t &store) {
    auto nname = std::string("n") + prefix + "hits";
    auto nleaf = nname + "/I";
    auto count = nname;
    tout->Branch(nname.c_str(), &store.nhits, nleaf.c_str());
    tout->Branch((prefix + "_frame_start").c_str(), store.frame_start.get(), (prefix + "_frame_start[nframes]/I").c_str());
    tout->Branch((prefix + "_frame_nhits").c_str(), store.frame_nhits.get(), (prefix + "_frame_nhits[nframes]/I").c_str());
    tout->Branch((prefix + "_device").c_str(), store.device.get(), (prefix + "_device[" + count + "]/I").c_str());
    tout->Branch((prefix + "_fifo").c_str(), store.fifo.get(), (prefix + "_fifo[" + count + "]/I").c_str());
    tout->Branch((prefix + "_type").c_str(), store.type.get(), (prefix + "_type[" + count + "]/I").c_str());
    tout->Branch((prefix + "_counter").c_str(), store.counter.get(), (prefix + "_counter[" + count + "]/I").c_str());
    tout->Branch((prefix + "_column").c_str(), store.column.get(), (prefix + "_column[" + count + "]/I").c_str());
    tout->Branch((prefix + "_pixel").c_str(), store.pixel.get(), (prefix + "_pixel[" + count + "]/I").c_str());
    tout->Branch((prefix + "_tdc").c_str(), store.tdc.get(), (prefix + "_tdc[" + count + "]/I").c_str());
    tout->Branch((prefix + "_rollover").c_str(), store.rollover.get(), (prefix + "_rollover[" + count + "]/I").c_str());
    tout->Branch((prefix + "_coarse").c_str(), store.coarse.get(), (prefix + "_coarse[" + count + "]/I").c_str());
    tout->Branch((prefix + "_fine").c_str(), store.fine.get(), (prefix + "_fine[" + count + "]/I").c_str());
    tout->Branch((prefix + "_time").c_str(), store.time.get(), (prefix + "_time[" + count + "]/D").c_str());
  };

  branch_store("trigger", spill.trigger);
  branch_store("timing", spill.timing);
  branch_store("cherenkov", spill.cherenkov);

  spill_participation_t spill_meta;
  auto tmeta_out = new TTree("spill_participation", "spill_participation");
  tmeta_out->Branch("spill", &spill_meta.spill, "spill/I");
  tmeta_out->Branch("counter", &spill_meta.counter, "counter/I");
  tmeta_out->Branch("nsources", &spill_meta.nsources, "nsources/I");
  tmeta_out->Branch("source_device", spill_meta.source_device, "source_device[nsources]/I");
  tmeta_out->Branch("source_fifo", spill_meta.source_fifo, "source_fifo[nsources]/I");

  auto load_spill_participation = [&](const data_t &start_word, int sequential_spill) {
    if (!tmeta_in) {
      spill_meta.clear(sequential_spill, start_word.counter);
      return spill_meta.add(start_word.device, start_word.fifo);
    }

    if (meta_entry >= meta_entries) {
      std::cerr << "ERROR: missing spill_participation entry for input spill counter "
                << start_word.counter << std::endl;
      return false;
    }

    auto bytes = tmeta_in->GetEntry(meta_entry++);
    if (bytes <= 0) {
      std::cerr << "ERROR: failed to read spill_participation entry "
                << (meta_entry - 1) << std::endl;
      return false;
    }

    if (meta_in.counter != start_word.counter) {
      std::cerr << "ERROR: spill_participation counter mismatch"
                << " data_counter=" << start_word.counter
                << " meta_counter=" << meta_in.counter << std::endl;
      return false;
    }
    if (meta_in.nsources < 0 || meta_in.nsources > maxsources) {
      std::cerr << "ERROR: invalid spill_participation nsources="
                << meta_in.nsources << std::endl;
      return false;
    }

    spill_meta = meta_in;
    return true;
  };

  std::vector<data_t> buffer;
  std::vector<active_frame_t> frames;

  int spill_id = -1;
  int nspills = 0;
  int ntriggers = 0;
  int nwritten = 0;
  int noverflow = 0;
  int ndropped = 0;
  long long nunexpected_words = 0;
  bool in_spill = false;
  double last_time = 0.;
  double history_window = definition.history_window(window);

  auto append_frame = [&](active_frame_t &active) {
    auto before = spill.nunexpected_words;
    if (!spill.add(active.frame)) {
      nunexpected_words += spill.nunexpected_words - before;
      std::cerr << " --- spill exceeds maxframes=" << maxframes
                << " or category maxspillhits=" << maxspillhits
                << ", dropping frame" << std::endl;
      ++ndropped;
      return;
    }
    nunexpected_words += spill.nunexpected_words - before;
    ++nwritten;
  };

  auto flush_frames = [&]() {
    for (auto &active : frames) {
      if (!active.decided && last_time >= active.decision_time) {
        active.accepted = decide(active, definition);
        active.decided = true;
      }
      if (active.decided && active.accepted)
        append_frame(active);
    }
    frames.clear();
  };

  auto fill_spill = [&]() {
    if (!in_spill)
      return;
    flush_frames();
    tout->Fill();
    tmeta_out->Fill();
  };

  coincidence_cluster_t coincidence_cluster;

  auto reset = [&]() {
    buffer.clear();
    frames.clear();
    coincidence_cluster.clear();
  };

  auto add_hit_to_active_frame = [&](active_frame_t &active, const data_t &hit) {
    auto frame_delta_t = hit.time - active.frame_center;
    if (frame_delta_t >= -window && frame_delta_t <= window && !active.frame.add(hit)) {
      std::cerr << " --- frame exceeds maxhits=" << maxhits << ", truncating" << std::endl;
      ++noverflow;
    }
    update_conditions(active, hit, definition);
  };

  auto create_candidate = [&](double event_time, int seed_entry, const data_t *seed) {
    active_frame_t active;
    active.event_time = event_time;
    active.frame_center = event_time + definition.frame_event_offset();
    active.decision_time = definition.decision_latest_raw_time(event_time);
    active.seed_entry = seed_entry;
    active.frame.nhits = 0;
    active.counts.assign(definition.conditions.size(), 0);

    if (seed) {
      if (!active.frame.add(*seed)) {
        std::cerr << " --- frame exceeds maxhits=" << maxhits << ", truncating" << std::endl;
        ++noverflow;
      }
      update_conditions(active, *seed, definition);
    }

    for (auto &hit : buffer)
      add_hit_to_active_frame(active, hit);

    if (!seed || seed_entry < 0)
      add_hit_to_active_frame(active, data);

    if (data.time >= active.decision_time) {
      active.accepted = decide(active, definition);
      active.decided = true;
    }

    frames.push_back(active);
    ++ntriggers;
  };

  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);

    /** start of spill is detected **/
    if (data.is_start_spill()) {
      fill_spill();
      reset();
      /** do any needed reset action **/
      spill_id = data.counter;
      if (!load_spill_participation(data, nspills)) {
        fout->Close();
        fin->Close();
        return false;
      }
      ++nspills;
      spill.reset(spill_id);
      in_spill = true;
      continue;
    }

    /** end of spill is detected **/
    if (data.is_end_spill()) {
      fill_spill();
      /** do any needed reset action **/
      reset();
      in_spill = false;
      continue;
    }

    /** not an ALCOR hit and not a trigger tag **/
    if (!data.is_alcor_hit() && !data.is_trigger_tag())
      continue;

    /** compute and store time **/
    if (!data.has_time) {
      data.set_nominal_time();
    }
    last_time = data.time;

    while (!buffer.empty() && data.time - buffer.front().time > history_window)
      buffer.erase(buffer.begin());

    for (int i = 0; i < (int)frames.size();) {
      auto &active = frames[i];
      if (iev != active.seed_entry)
        add_hit_to_active_frame(active, data);

      if (!active.decided && data.time > active.decision_time) {
        active.accepted = decide(active, definition);
        active.decided = true;
      }

      if (active.decided && data.time > active.frame_center + window) {
        if (active.accepted)
          append_frame(active);
        frames.erase(frames.begin() + i);
        continue;
      }
      ++i;
    }

    if (definition.mode == trigger_mode_t::coincidence && definition.coincidence.selector.match(data)) {
      auto corrected_time = data.time - definition.coincidence.selector.time_offset;
      if (!coincidence_cluster.active ||
          corrected_time > coincidence_cluster.reference_time + definition.coincidence.dtmax) {
        coincidence_cluster.reset(corrected_time);
      }

      auto delta_t = corrected_time - coincidence_cluster.reference_time;
      if (delta_t >= definition.coincidence.dtmin && delta_t <= definition.coincidence.dtmax)
        coincidence_cluster.add(data);

      if (!coincidence_cluster.fired &&
          coincidence_cluster.count(definition.coincidence.unique_channel) >= definition.coincidence.min_count) {
        create_candidate(coincidence_cluster.reference_time, -1, nullptr);
        coincidence_cluster.fired = true;
      }
    }

    if (definition.mode == trigger_mode_t::seeded && definition.seed.match(data))
      create_candidate(data.time - definition.seed.time_offset, iev, &data);

    buffer.push_back(data);
  }

  fill_spill();

  if (tmeta_in && meta_entry != meta_entries) {
    std::cerr << "ERROR: not all input spill_participation entries were consumed" << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  auto nmeta_out = tmeta_out->GetEntries();
  if (nmeta_out != tout->GetEntries()) {
    std::cerr << "ERROR: frames/spill_participation entry-count mismatch" << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << " --- spills processed: " << nspills << std::endl;
  std::cout << " --- trigger candidates found: " << ntriggers << std::endl;
  std::cout << " --- frames written: " << nwritten << std::endl;
  std::cout << " --- spill_participation entries: " << nmeta_out << std::endl;
  if (noverflow > 0)
    std::cout << " --- frames truncated: " << noverflow << std::endl;
  if (ndropped > 0)
    std::cout << " --- frames dropped: " << ndropped << std::endl;
  if (nunexpected_words > 0)
    std::cerr << "WARNING: discarded " << nunexpected_words
              << " unexpected words while writing triggered frames" << std::endl;
  return true;
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string output;
  std::string config;
  double window;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("config,c", po::value<std::string>(&config)->required(), "trigger configuration file")
    ("window,w", po::value<double>(&window)->default_value(256.), "frame window in clock cycles");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, options), vm);
    if (vm.count("help")) {
      std::cout << options << std::endl;
      return 0;
    }
    po::notify(vm);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << options << std::endl;
    return 1;
  }

  return trigger(input, output, config, window) ? 0 : 1;
}
