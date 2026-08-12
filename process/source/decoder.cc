#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct main_header_t {
  uint32_t caffe;
  uint32_t readout_version;
  uint32_t firmware_release;
  uint32_t run_number;
  uint32_t timestamp;
  uint32_t staging_size;
  uint32_t run_mode;
  uint32_t filter_mode;
  uint32_t device;
  uint32_t reserved1;
  uint32_t reserved2;
  uint32_t reserved3;
  uint32_t reserved4;
  uint32_t reserved5;
  uint32_t reserved6;
  uint32_t reserved7;
};

struct buffer_header_t {
  uint32_t caffe;
  uint32_t id;
  uint32_t counter;
  uint32_t size;
};

struct alcor_hit_t {
  uint32_t fine   : 9;
  uint32_t coarse : 15;
  uint32_t tdc    : 2;
  uint32_t pixel  : 3;
  uint32_t column : 3;
};

struct data_t {
  int device = 0;
  int fifo = 0;
  int type = 0;
  int counter = 0;
  int column = 0;
  int pixel = 0;
  int tdc = 0;
  int rollover = 0;
  int coarse = 0;
  int fine = 0;
};

struct stats_t {
  uint64_t buffers = 0;
  uint64_t spills_found = 0;
  uint64_t spills_written = 0;
  uint64_t spills_emptied = 0;
  uint64_t spills_suppressed_by_errors = 0;
  uint64_t payload_words_suppressed_by_errors = 0;
  uint64_t spills_incomplete = 0;
  uint64_t alcor_hits = 0;
  uint64_t trigger_tags = 0;
  uint64_t start_spill = 0;
  uint64_t end_spill = 0;
  uint64_t skipped_outside_spill = 0;
  uint64_t wrong_end_spill = 0;
  uint64_t unexpected_start_spill = 0;
  uint64_t invalid_column_hits = 0;
  uint64_t malformed_words = 0;
  uint64_t killed_fifo = 0;
  uint64_t daq_suppressed_records = 0;
  uint64_t incomplete_spills_discarded = 0;
  uint64_t incomplete_payload_words_discarded = 0;
};

struct fifo_state_t {
  bool in_spill = false;
  int rollover = 0;
  int spill_errors = 0;
  data_t start;
  std::vector<data_t> payload;
};

static bool verbose = false;
static data_t data;

static bool
is_start_spill(uint32_t word)
{
  return (word & 0xf0000000u) == 0x70000000u;
}

static bool
is_end_spill(uint32_t word)
{
  return (word & 0xf0000000u) == 0xf0000000u;
}

static bool
is_trigger_tag(uint32_t word)
{
  return (word & 0xf0000000u) == 0x90000000u;
}

static bool
is_deadbeef(uint32_t word)
{
  return word == 0xdeadbeefu;
}

static int
control_counter(uint32_t word)
{
  return (word & 0x0fff0000u) >> 16;
}

static void
time_from_control(const uint32_t *words, int &rollover, int &coarse)
{
  uint64_t trigger_time = (uint64_t)(words[0] & 0xffu) << 32;
  trigger_time |= words[1];
  coarse = trigger_time & 0x7fff;
  rollover = trigger_time >> 15;
}

static data_t
make_control(int device, int fifo, int type, int counter, int rollover, int coarse)
{
  data_t out;
  out.device = device;
  out.fifo = fifo;
  out.type = type;
  out.counter = counter;
  out.column = 0;
  out.pixel = 0;
  out.tdc = 0;
  out.rollover = rollover;
  out.coarse = coarse;
  out.fine = 0;
  return out;
}

static data_t
make_hit(int device, int fifo, const alcor_hit_t &hit, int rollover)
{
  data_t out;
  out.device = device;
  out.fifo = fifo;
  out.type = 1;
  out.counter = 0;
  out.column = hit.column;
  out.pixel = hit.pixel;
  out.tdc = hit.tdc;
  out.rollover = rollover;
  out.coarse = hit.coarse;
  out.fine = hit.fine;
  return out;
}

static bool
valid_column(int fifo, int column)
{
  int first = 2 * (fifo % 4);
  return column == first || column == first + 1;
}

static void
fill(TTree *tree, const data_t &word)
{
  data = word;
  tree->Fill();
}


static bool
is_daq_suppressed_record(const uint32_t *words, uint32_t nwords)
{
  if (nwords != 4)
    return false;
  if (!is_start_spill(words[0]) || !is_deadbeef(words[1]))
    return false;
  return is_end_spill(words[2]) && is_deadbeef(words[3]);
}

static void
write_spill(TTree *tree,
            const data_t &start,
            const std::vector<data_t> &payload,
            const data_t &end,
            int spill_errors,
            int allowed_spill_errors,
            stats_t &stats)
{
  bool suppress_spill = spill_errors > allowed_spill_errors;
  if (suppress_spill) {
    ++stats.spills_emptied;
    ++stats.spills_suppressed_by_errors;
    stats.payload_words_suppressed_by_errors += payload.size();
    return;
  }

  fill(tree, start);
  ++stats.start_spill;

  for (const auto &word : payload) {
    fill(tree, word);
    if (word.type == 1) ++stats.alcor_hits;
    if (word.type == 9) ++stats.trigger_tags;
  }

  fill(tree, end);
  ++stats.end_spill;
  ++stats.spills_written;
}

static bool
decode_alcor_buffer(const uint32_t *words,
                    uint32_t nwords,
                    int device,
                    int fifo,
                    int allowed_spill_errors,
                    TTree *tree,
                    stats_t &stats,
                    fifo_state_t &state)
{
  if (!state.in_spill && is_daq_suppressed_record(words, nwords)) {
    ++stats.daq_suppressed_records;
    if (verbose) {
      std::cerr << "WARNING: skipping DAQ-suppressed deadbeef record fifo=" << fifo
                << " start_counter=" << control_counter(words[0])
                << " end_counter=" << control_counter(words[2])
                << std::endl;
    }
    return true;
  }

  uint32_t pos = 0;

  while (pos < nwords) {
    if (!state.in_spill) {
      if (!is_start_spill(words[pos])) {
        ++stats.skipped_outside_spill;
        ++pos;
        continue;
      }

      if (pos + 1 >= nwords) {
        ++stats.spills_incomplete;
        std::cerr << "ERROR: START_SPILL without continuation word, fifo=" << fifo
                  << " pos=" << pos << std::endl;
        return false;
      }

      int start_counter = control_counter(words[pos]);
      int start_rollover = 0;
      int start_coarse = 0;
      time_from_control(&words[pos], start_rollover, start_coarse);
      state.start = make_control(device, fifo, 7, start_counter, start_rollover, start_coarse);
      state.payload.clear();
      state.spill_errors = 0;
      state.rollover = 0;
      state.in_spill = true;
      ++stats.spills_found;
      pos += 2;
      continue;
    }

    uint32_t word = words[pos];

    if (word == 0x5c5c5c5cu) {
      ++state.rollover;
      ++pos;
      continue;
    }

    if (word == 0x666caffeu) {
      ++stats.killed_fifo;
      ++state.spill_errors;
      data_t end = make_control(device, fifo, 15, state.start.counter, state.rollover, 0);
      write_spill(tree, state.start, state.payload, end, state.spill_errors,
                  allowed_spill_errors, stats);
      state.in_spill = false;
      state.payload.clear();
      ++pos;
      continue;
    }

    if (is_start_spill(word)) {
      ++stats.unexpected_start_spill;
      ++state.spill_errors;
      if (verbose)
        std::cerr << "WARNING: unexpected START_SPILL inside spill fifo=" << fifo
                  << " open=" << state.start.counter
                  << " candidate=" << control_counter(word)
                  << " pos=" << pos << std::endl;
      ++pos;
      continue;
    }

    if (is_end_spill(word)) {
      int end_counter = control_counter(word);
      if (end_counter != state.start.counter) {
        ++stats.wrong_end_spill;
        ++state.spill_errors;
        if (verbose)
          std::cerr << "WARNING: wrong END_SPILL candidate fifo=" << fifo
                    << " open=" << state.start.counter
                    << " candidate=" << end_counter
                    << " word=0x" << std::hex << word << std::dec
                    << " pos=" << pos << std::endl;
        ++pos;
        continue;
      }

      if (pos + 1 >= nwords) {
        ++stats.spills_incomplete;
        std::cerr << "ERROR: END_SPILL without continuation word, fifo=" << fifo
                  << " counter=" << end_counter << std::endl;
        return false;
      }

      int end_rollover = 0;
      int end_coarse = 0;
      time_from_control(&words[pos], end_rollover, end_coarse);
      data_t end = make_control(device, fifo, 15, end_counter, end_rollover, end_coarse);
      write_spill(tree, state.start, state.payload, end, state.spill_errors,
                  allowed_spill_errors, stats);
      state.in_spill = false;
      state.payload.clear();
      pos += 2;
      continue;
    }

    alcor_hit_t hit = *(const alcor_hit_t *)&word;
    if (!valid_column(fifo, hit.column)) {
      ++stats.invalid_column_hits;
      ++state.spill_errors;
      if (verbose)
        std::cerr << "WARNING: invalid column hit fifo=" << fifo
                  << " column=" << hit.column
                  << " word=0x" << std::hex << word << std::dec
                  << " pos=" << pos << std::endl;
      ++pos;
      continue;
    }

    state.payload.push_back(make_hit(device, fifo, hit, state.rollover));
    ++pos;
  }

  return true;
}

static bool
decode_trigger_buffer(const uint32_t *words,
                      uint32_t nwords,
                      int device,
                      int fifo,
                      int allowed_spill_errors,
                      TTree *tree,
                      stats_t &stats,
                      fifo_state_t &state)
{
  uint32_t pos = 0;

  while (pos < nwords) {
    if (!state.in_spill) {
      if (!is_start_spill(words[pos])) {
        ++stats.skipped_outside_spill;
        ++pos;
        continue;
      }
      if (pos + 1 >= nwords) {
        ++stats.spills_incomplete;
        return false;
      }
      int start_counter = control_counter(words[pos]);
      int start_rollover = 0;
      int start_coarse = 0;
      time_from_control(&words[pos], start_rollover, start_coarse);
      state.start = make_control(device, fifo, 7, start_counter, start_rollover, start_coarse);
      state.payload.clear();
      state.spill_errors = 0;
      state.in_spill = true;
      ++stats.spills_found;
      pos += 2;
      continue;
    }

    uint32_t word = words[pos];

    if (is_trigger_tag(word)) {
      if (pos + 1 >= nwords) {
        ++stats.malformed_words;
        ++state.spill_errors;
        ++pos;
        continue;
      }
      int trigger_rollover = 0;
      int trigger_coarse = 0;
      int counter = (word & 0x00ffff00u) >> 8;
      time_from_control(&words[pos], trigger_rollover, trigger_coarse);
      state.payload.push_back(make_control(device, fifo, 9, counter, trigger_rollover, trigger_coarse));
      pos += 2;
      continue;
    }

    if (is_start_spill(word)) {
      ++stats.unexpected_start_spill;
      ++state.spill_errors;
      ++pos;
      continue;
    }

    if (is_end_spill(word)) {
      int end_counter = control_counter(word);
      if (end_counter != state.start.counter) {
        ++stats.wrong_end_spill;
        ++state.spill_errors;
        ++pos;
        continue;
      }
      if (pos + 1 >= nwords) {
        ++stats.spills_incomplete;
        return false;
      }
      int end_rollover = 0;
      int end_coarse = 0;
      time_from_control(&words[pos], end_rollover, end_coarse);
      data_t end = make_control(device, fifo, 15, end_counter, end_rollover, end_coarse);
      write_spill(tree, state.start, state.payload, end, state.spill_errors,
                  allowed_spill_errors, stats);
      state.in_spill = false;
      state.payload.clear();
      pos += 2;
      continue;
    }

    ++stats.malformed_words;
    ++state.spill_errors;
    ++pos;
  }

  return true;
}


static std::string
summary_filename(const std::string &output)
{
  const std::string suffix = ".root";
  if (output.size() >= suffix.size() &&
      output.compare(output.size() - suffix.size(), suffix.size(), suffix) == 0)
    return output.substr(0, output.size() - suffix.size()) + ".summary";
  return output + ".summary";
}

static std::string
make_summary(const std::string &input,
             const std::string &output,
             Long64_t output_entries,
             const stats_t &stats)
{
  std::ostringstream out;
  out << "input: " << input << '\n';
  out << "output: " << output << '\n';
  out << "buffers: " << stats.buffers << '\n';
  out << "spills_found: " << stats.spills_found << '\n';
  out << "spills_written: " << stats.spills_written << '\n';
  out << "spills_emptied: " << stats.spills_emptied << '\n';
  out << "spills_suppressed_by_errors: " << stats.spills_suppressed_by_errors << '\n';
  out << "payload_words_suppressed_by_errors: " << stats.payload_words_suppressed_by_errors << '\n';
  out << "spills_incomplete: " << stats.spills_incomplete << '\n';
  out << "start_spill_type7: " << stats.start_spill << '\n';
  out << "end_spill_type15: " << stats.end_spill << '\n';
  out << "alcor_hits_type1: " << stats.alcor_hits << '\n';
  out << "trigger_tags_type9: " << stats.trigger_tags << '\n';
  out << "skipped_outside_spill: " << stats.skipped_outside_spill << '\n';
  out << "wrong_end_spill_candidates: " << stats.wrong_end_spill << '\n';
  out << "unexpected_start_spill: " << stats.unexpected_start_spill << '\n';
  out << "invalid_column_hits: " << stats.invalid_column_hits << '\n';
  out << "malformed_words: " << stats.malformed_words << '\n';
  out << "killed_fifo_markers: " << stats.killed_fifo << '\n';
  out << "daq_suppressed_records: " << stats.daq_suppressed_records << '\n';
  out << "incomplete_spills_discarded: " << stats.incomplete_spills_discarded << '\n';
  out << "incomplete_payload_words_discarded: " << stats.incomplete_payload_words_discarded << '\n';
  out << "output_entries: " << output_entries << '\n';
  return out.str();
}

static bool
write_summary_file(const std::string &filename, const std::string &summary)
{
  std::ofstream out(filename);
  if (!out) {
    std::cerr << "ERROR: could not create summary file: " << filename << std::endl;
    return false;
  }
  out << summary;
  return true;
}

static bool
decoder(const std::string &input,
        const std::string &output,
        int allowed_spill_errors)
{
  std::ifstream fin(input, std::ios::binary);
  if (!fin) {
    std::cerr << "ERROR: could not open input file: " << input << std::endl;
    return false;
  }

  main_header_t main_header;
  fin.read((char *)&main_header, sizeof(main_header));
  if (!fin || main_header.caffe != 0x000caffeu) {
    std::cerr << "ERROR: invalid main header" << std::endl;
    return false;
  }

  if (main_header.filter_mode != 0x0 && main_header.filter_mode != 0xf) {
    std::cerr << "ERROR: unsupported filter mode: 0x" << std::hex
              << main_header.filter_mode << std::dec << std::endl;
    return false;
  }

  auto fout = TFile::Open(output.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << output << std::endl;
    return false;
  }

  auto tree = new TTree("alcor", "ALCOR");
  tree->Branch("device", &data.device, "device/I");
  tree->Branch("fifo", &data.fifo, "fifo/I");
  tree->Branch("type", &data.type, "type/I");
  tree->Branch("counter", &data.counter, "counter/I");
  tree->Branch("column", &data.column, "column/I");
  tree->Branch("pixel", &data.pixel, "pixel/I");
  tree->Branch("tdc", &data.tdc, "tdc/I");
  tree->Branch("rollover", &data.rollover, "rollover/I");
  tree->Branch("coarse", &data.coarse, "coarse/I");
  tree->Branch("fine", &data.fine, "fine/I");

  stats_t stats;
  std::vector<fifo_state_t> states(100);
  std::vector<char> buffer(main_header.staging_size);

  while (true) {
    buffer_header_t header;
    fin.read((char *)&header, sizeof(header));
    if (fin.eof())
      break;
    if (!fin || header.caffe != 0x123caffeu) {
      std::cerr << "ERROR: invalid buffer header" << std::endl;
      fout->Close();
      return false;
    }
    if (header.size > buffer.size())
      buffer.resize(header.size);
    fin.read(buffer.data(), header.size);
    if (!fin) {
      std::cerr << "ERROR: could not read full buffer payload" << std::endl;
      fout->Close();
      return false;
    }

    ++stats.buffers;
    const auto *words = (const uint32_t *)buffer.data();
    uint32_t nwords = header.size / 4;

    bool ok = true;
    if (header.id < 32)
      ok = decode_alcor_buffer(words, nwords, main_header.device, header.id,
                               allowed_spill_errors, tree, stats, states.at(header.id));
    else if (header.id == 32 || header.id == 99)
      ok = decode_trigger_buffer(words, nwords, main_header.device, header.id,
                                  allowed_spill_errors, tree, stats, states.at(header.id));

    if (!ok) {
      fout->Close();
      return false;
    }
  }

  for (size_t i = 0; i < states.size(); ++i) {
    if (states[i].in_spill) {
      ++stats.spills_incomplete;
      ++stats.incomplete_spills_discarded;
      stats.incomplete_payload_words_discarded += states[i].payload.size();
      std::cerr << "WARNING: discarding incomplete spill at EOF fifo=" << i
                << " counter=" << states[i].start.counter
                << " payload_words=" << states[i].payload.size() << std::endl;
      states[i].in_spill = false;
      states[i].payload.clear();
    }
  }

  auto output_entries = tree->GetEntries();
  tree->Write();
  fout->Close();

  auto summary = make_summary(input, output, output_entries, stats);
  auto sfilename = summary_filename(output);
  if (!write_summary_file(sfilename, summary))
    return false;

  std::cout << summary;
  std::cout << "summary: " << sfilename << std::endl;

  return true;
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string output;
  int allowed_spill_errors = 0;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input DAT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("allowed-spill-errors", po::value<int>(&allowed_spill_errors)->default_value(0),
     "maximum number of decoding errors allowed before a spill payload is discarded")
    ("verbose,v", po::bool_switch(&verbose)->default_value(false), "verbose diagnostics");

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

  if (allowed_spill_errors < 0) {
    std::cerr << "ERROR: --allowed-spill-errors must be >= 0" << std::endl;
    return 1;
  }

  return decoder(input, output, allowed_spill_errors) ? 0 : 1;
}
