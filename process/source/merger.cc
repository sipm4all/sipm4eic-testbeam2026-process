#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>


struct stream_t {
  std::string filename;
  std::unique_ptr<TFile> file;
  TTree *tree = nullptr;
  Long64_t entries = 0;
  Long64_t entry = 0;
  data_t data;
  bool valid = false;
  bool read_error = false;
  Long64_t start_spills = 0;
  Long64_t end_spills = 0;

  bool open(const std::string &_filename)
  {
    filename = _filename;
    file.reset(TFile::Open(filename.c_str()));
    if (!file || file->IsZombie()) {
      std::cerr << "ERROR: could not open input file: " << filename << std::endl;
      return false;
    }

    tree = (TTree *)file->Get("alcor");
    if (!tree) {
      std::cerr << "ERROR: could not find 'alcor' tree in input file: " << filename << std::endl;
      return false;
    }

    entries = tree->GetEntries();
    data.link_to_tree(tree);
    return next();
  }

  bool next()
  {
    if (entry >= entries) {
      valid = false;
      return true;
    }

    auto current = entry;
    auto bytes = tree->GetEntry(entry++);
    if (bytes <= 0) {
      std::cerr << "ERROR: ROOT GetEntry failed: file=" << filename
                << " entry=" << current << std::endl;
      valid = false;
      read_error = true;
      return false;
    }

    if (data.is_data_word())
      data.set_time();
    valid = true;
    return true;
  }

  Long64_t current_entry() const { return valid ? entry - 1 : entry; }
};

struct output_t {
  TFile *file = nullptr;
  TTree *tree = nullptr;
  data_t data;
  std::string tmpfilename;
  std::string finalfilename;
  Long64_t total_entries = 0;
  int nfiles = 0;

  bool open(TTree *template_tree, const std::string &final)
  {
    finalfilename = final;
    tmpfilename = finalfilename + ".tmp";
    std::remove(tmpfilename.c_str());

    file = TFile::Open(tmpfilename.c_str(), "RECREATE");
    if (!file || file->IsZombie()) {
      std::cerr << "ERROR: could not create output file: " << tmpfilename << std::endl;
      close_failed();
      return false;
    }

    tree = template_tree->CloneTree(0);
    data.link_to_tree(tree);
    return true;
  }

  void fill(const data_t &word)
  {
    data = word;
    tree->Fill();
  }

  bool close_success()
  {
    if (!file)
      return true;

    auto entries = tree ? tree->GetEntries() : 0;
    if (file->Write() <= 0) {
      std::cerr << "ERROR: failed writing output file: " << tmpfilename << std::endl;
      close_failed();
      return false;
    }

    file->Close();
    delete file;
    file = nullptr;
    tree = nullptr;

    std::remove(finalfilename.c_str());
    if (std::rename(tmpfilename.c_str(), finalfilename.c_str()) != 0) {
      std::cerr << "ERROR: could not rename temporary output " << tmpfilename
                << " to " << finalfilename << std::endl;
      std::remove(tmpfilename.c_str());
      return false;
    }

    total_entries += entries;
    ++nfiles;
    return true;
  }

  void close_failed()
  {
    if (file) {
      file->Close();
      delete file;
    }
    file = nullptr;
    tree = nullptr;
    if (!tmpfilename.empty())
      std::remove(tmpfilename.c_str());
  }
};

static void
print_streams(const std::vector<stream_t> &streams)
{
  for (int i = 0; i < (int)streams.size(); ++i) {
    auto &s = streams[i];
    std::cerr << "stream " << i
              << ": file=" << s.filename
              << " entry=" << s.current_entry() << "/" << s.entries
              << " valid=" << s.valid;
    if (s.valid)
      std::cerr << " type=" << s.data.type
                << " device=" << s.data.device
                << " fifo=" << s.data.fifo;
    if (s.read_error)
      std::cerr << " read_error=1";
    std::cerr << std::endl;
  }
}

static std::string
spill_filename(const std::string &outfilename, Long64_t spill)
{
  std::string base = outfilename;
  if (base.size() >= 5 && base.substr(base.size() - 5) == ".root")
    base.resize(base.size() - 5);

  std::ostringstream ss;
  ss << base << ".spill_" << std::setw(4) << std::setfill('0') << spill << ".root";
  return ss.str();
}

bool
merger(const std::vector<std::string> filenames, const std::string outfilename, bool split_spills)
{
  if (filenames.empty()) {
    std::cerr << "ERROR: no input files provided" << std::endl;
    return false;
  }

  std::vector<stream_t> streams(filenames.size());
  Long64_t nsum = 0;
  for (int i = 0; i < (int)filenames.size(); ++i) {
    if (!streams[i].open(filenames[i]))
      return false;
    nsum += streams[i].entries;
  }

  output_t output;
  auto fail = [&](const std::string &message) {
    std::cerr << "ERROR: " << message << std::endl;
    print_streams(streams);
    output.close_failed();
    return false;
  };

  if (!split_spills && !output.open(streams[0].tree, outfilename))
    return false;

  auto all_eof = [&]() {
    for (auto &stream : streams)
      if (stream.valid || stream.entry != stream.entries)
        return false;
    return true;
  };

  auto any_eof = [&]() {
    for (auto &stream : streams)
      if (!stream.valid && stream.entry >= stream.entries)
        return true;
    return false;
  };

  auto check_eof_consistency = [&]() {
    if (!any_eof())
      return true;
    if (all_eof())
      return true;
    std::cerr << "ERROR: some streams reached EOF while others still contain entries" << std::endl;
    print_streams(streams);
    return false;
  };

  auto all_start_spill = [&]() {
    for (auto &stream : streams)
      if (!stream.valid || !stream.data.is_start_spill())
        return false;
    return true;
  };

  auto all_end_spill = [&]() {
    for (auto &stream : streams)
      if (!stream.valid || !stream.data.is_end_spill())
        return false;
    return true;
  };

  auto best_stream = [&]() {
    int ibest = -1;
    for (int i = 0; i < (int)streams.size(); ++i) {
      auto &stream = streams[i];
      if (!stream.valid || !stream.data.is_data_word())
        continue;
      if (ibest < 0 || stream.data.time < streams[ibest].data.time)
        ibest = i;
    }
    return ibest;
  };

  Long64_t nspill = 0;
  bool in_spill = false;

  for (;;) {
    if (!check_eof_consistency())
      return fail("inconsistent EOF state");

    if (all_eof()) {
      if (in_spill)
        return fail("all streams reached EOF inside an incomplete spill");
      break;
    }

    if (!in_spill) {
      if (!all_start_spill())
        return fail("streams are not synchronized at START_SPILL");

      if (split_spills && !output.open(streams[0].tree, spill_filename(outfilename, nspill)))
        return fail("could not open split-spill output file");

      output.fill(streams[0].data);
      for (auto &stream : streams) {
        ++stream.start_spills;
        if (!stream.next())
          return fail("failed while advancing past START_SPILL");
      }
      ++nspill;
      in_spill = true;
      continue;
    }

    int ibest = best_stream();
    if (ibest >= 0) {
      output.fill(streams[ibest].data);
      if (!streams[ibest].next())
        return fail("failed while advancing data stream");
      continue;
    }

    if (all_end_spill()) {
      output.fill(streams[0].data);
      for (auto &stream : streams) {
        ++stream.end_spills;
        if (!stream.next())
          return fail("failed while advancing past END_SPILL");
      }
      in_spill = false;

      if (split_spills && !output.close_success())
        return fail("could not close split-spill output file successfully");

      continue;
    }

    for (auto &stream : streams) {
      if (stream.valid && !stream.data.is_start_spill() && !stream.data.is_end_spill() && !stream.data.is_data_word())
        return fail("unexpected word type inside spill: type=" + std::to_string(stream.data.type));
    }

    return fail("streams are not synchronized at END_SPILL");
  }

  for (auto &stream : streams) {
    if (stream.entry != stream.entries)
      return fail("not all input entries were consumed");
  }

  for (int i = 1; i < (int)streams.size(); ++i) {
    if (streams[i].start_spills != streams[0].start_spills)
      return fail("START_SPILL count mismatch across streams");
    if (streams[i].end_spills != streams[0].end_spills)
      return fail("END_SPILL count mismatch across streams");
  }

  for (auto &stream : streams) {
    if (stream.start_spills != stream.end_spills)
      return fail("START_SPILL/END_SPILL count mismatch within stream");
  }

  Long64_t expected_output = nsum - 2 * nspill * ((Long64_t)streams.size() - 1);
  Long64_t actual_output = split_spills ? output.total_entries : output.tree->GetEntries();
  if (actual_output != expected_output) {
    std::cerr << "ERROR: output entry-count mismatch" << std::endl;
    std::cerr << "Nsum=" << nsum << std::endl;
    std::cerr << "Ninput=" << streams.size() << std::endl;
    std::cerr << "Nspill=" << nspill << std::endl;
    std::cerr << "expected_output=" << expected_output << std::endl;
    std::cerr << "actual_output=" << actual_output << std::endl;
    return fail("output entry-count mismatch");
  }

  if (!split_spills && !output.close_success())
    return fail("could not close output file successfully");

  std::cout << " --- merge successful" << std::endl;
  std::cout << " --- merged spills: " << nspill << std::endl;
  if (split_spills)
    std::cout << " --- output files: " << output.nfiles << std::endl;
  std::cout << " --- input entries: " << nsum << std::endl;
  std::cout << " --- output entries: " << actual_output << std::endl;
  for (int i = 0; i < (int)streams.size(); ++i) {
    auto &stream = streams[i];
    std::cout << "stream " << i << ": "
              << stream.entry << " / " << stream.entries << " consumed, "
              << "start_spill=" << stream.start_spills << ", "
              << "end_spill=" << stream.end_spills << std::endl;
  }

  return true;
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::vector<std::string> input;
  std::string output;
  bool split_spills;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::vector<std::string>>(&input)->multitoken()->required(), "input ROOT files")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("split-spills", po::bool_switch(&split_spills), "write one output ROOT file per spill");

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

  return merger(input, output, split_spills) ? 0 : 1;
}
