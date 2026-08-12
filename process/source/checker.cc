#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct check_t {
  Long64_t entries = 0;
  Long64_t start_spill = 0;
  Long64_t end_spill = 0;
  Long64_t alcor_hits = 0;
  Long64_t trigger_tags = 0;
  Long64_t unknown_words = 0;

  bool in_spill = false;
  bool consistent = true;
  int current_spill = -1;
  int last_start_counter = -1;
  int last_end_counter = -1;
  std::vector<std::string> errors;

  void error(const std::string &message)
  {
    consistent = false;
    errors.push_back(message);
  }
};

static void
check_word(check_t &check, const data_t &data, Long64_t entry)
{
  if (data.is_start_spill()) {
    ++check.start_spill;

    if (check.in_spill) {
      check.error("entry " + std::to_string(entry) +
                  ": START_SPILL before previous END_SPILL, counter=" +
                  std::to_string(data.counter));
    }

    if (check.last_start_counter >= 0 && data.counter != check.last_start_counter + 1) {
      check.error("entry " + std::to_string(entry) +
                  ": START_SPILL counter jump from " +
                  std::to_string(check.last_start_counter) + " to " +
                  std::to_string(data.counter));
    }

    check.in_spill = true;
    check.current_spill = data.counter;
    check.last_start_counter = data.counter;
    return;
  }

  if (data.is_end_spill()) {
    ++check.end_spill;

    if (!check.in_spill) {
      check.error("entry " + std::to_string(entry) +
                  ": END_SPILL without open START_SPILL, counter=" +
                  std::to_string(data.counter));
    } else if (data.counter != check.current_spill) {
      check.error("entry " + std::to_string(entry) +
                  ": END_SPILL counter " + std::to_string(data.counter) +
                  " does not match open START_SPILL counter " +
                  std::to_string(check.current_spill));
    }

    if (check.last_end_counter >= 0 && data.counter != check.last_end_counter + 1) {
      check.error("entry " + std::to_string(entry) +
                  ": END_SPILL counter jump from " +
                  std::to_string(check.last_end_counter) + " to " +
                  std::to_string(data.counter));
    }

    check.in_spill = false;
    check.current_spill = -1;
    check.last_end_counter = data.counter;
    return;
  }

  if (data.is_alcor_hit()) {
    ++check.alcor_hits;
    return;
  }

  if (data.is_trigger_tag()) {
    ++check.trigger_tags;
    return;
  }

  ++check.unknown_words;
}

static bool
write_check(const std::string &outfilename,
            const std::string &filename,
            const check_t &check)
{
  std::ofstream out(outfilename);
  if (!out) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    return false;
  }

  out << "input: " << filename << '\n';
  out << "entries: " << check.entries << '\n';
  out << "start_spill_type7: " << check.start_spill << '\n';
  out << "end_spill_type15: " << check.end_spill << '\n';
  out << "alcor_hits_type1: " << check.alcor_hits << '\n';
  out << "trigger_tags_type9: " << check.trigger_tags << '\n';
  out << "unknown_words: " << check.unknown_words << '\n';
  out << "spill_counter_consistent: " << (check.consistent ? "yes" : "no") << '\n';
  out << "open_spill_at_eof: " << (check.in_spill ? "yes" : "no") << '\n';

  if (check.start_spill != check.end_spill)
    out << "spill_count_balance: no\n";
  else
    out << "spill_count_balance: yes\n";

  out << "errors: " << check.errors.size() << '\n';
  for (const auto &error : check.errors)
    out << "error: " << error << '\n';

  return true;
}

bool
checker(const std::string &filename, const std::string &outfilename)
{
  auto fin = TFile::Open(filename.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << filename << std::endl;
    return false;
  }

  auto tin = (TTree *)fin->Get("alcor");
  if (!tin) {
    std::cerr << "ERROR: could not find 'alcor' tree in input file" << std::endl;
    fin->Close();
    return false;
  }

  data_t data;
  data.link_to_tree(tin);

  check_t check;
  check.entries = tin->GetEntries();

  for (Long64_t iev = 0; iev < check.entries; ++iev) {
    auto bytes = tin->GetEntry(iev);
    if (bytes <= 0) {
      std::cerr << "ERROR: ROOT GetEntry failed for entry " << iev << std::endl;
      fin->Close();
      return false;
    }

    check_word(check, data, iev);
  }

  if (check.in_spill) {
    check.error("EOF reached with open spill counter " +
                std::to_string(check.current_spill));
  }

  if (check.start_spill != check.end_spill) {
    check.error("START_SPILL count " + std::to_string(check.start_spill) +
                " differs from END_SPILL count " +
                std::to_string(check.end_spill));
  }

  auto ok = write_check(outfilename, filename, check);
  fin->Close();

  if (!ok)
    return false;

  std::cout << "input:                    " << filename << std::endl;
  std::cout << "output:                   " << outfilename << std::endl;
  std::cout << "entries:                  " << check.entries << std::endl;
  std::cout << "start spill words:        " << check.start_spill << std::endl;
  std::cout << "end spill words:          " << check.end_spill << std::endl;
  std::cout << "ALCOR hits:               " << check.alcor_hits << std::endl;
  std::cout << "trigger tags:             " << check.trigger_tags << std::endl;
  std::cout << "unknown words:            " << check.unknown_words << std::endl;
  std::cout << "spill counters consistent: " << (check.consistent ? "yes" : "no") << std::endl;

  return true;
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string output;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ASCII check file");

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

  return checker(input, output) ? 0 : 1;
}
