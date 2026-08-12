#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>


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

struct stream_t {
  std::string filename;
  std::unique_ptr<TFile> file;
  TTree *tree = nullptr;
  TTree *meta_tree = nullptr;
  Long64_t entries = 0;
  Long64_t entry = 0;
  Long64_t meta_entries = 0;
  Long64_t meta_entry = 0;
  data_t data;
  spill_participation_t meta;
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

    meta_tree = (TTree *)file->Get("spill_participation");
    if (meta_tree) {
      meta_entries = meta_tree->GetEntries();
      meta_tree->SetBranchAddress("spill", &meta.spill);
      meta_tree->SetBranchAddress("counter", &meta.counter);
      meta_tree->SetBranchAddress("nsources", &meta.nsources);
      meta_tree->SetBranchAddress("source_device", meta.source_device);
      meta_tree->SetBranchAddress("source_fifo", meta.source_fifo);
    }

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

  bool load_current_participation(spill_participation_t &out)
  {
    if (!valid || !data.is_start_spill())
      return false;

    if (!meta_tree)
      return out.add(data.device, data.fifo);

    if (meta_entry >= meta_entries) {
      std::cerr << "ERROR: missing spill_participation entry: file=" << filename
                << " counter=" << data.counter << std::endl;
      return false;
    }

    auto current = meta_entry;
    auto bytes = meta_tree->GetEntry(meta_entry++);
    if (bytes <= 0) {
      std::cerr << "ERROR: ROOT GetEntry failed for spill_participation: file="
                << filename << " entry=" << current << std::endl;
      return false;
    }

    if (meta.counter != data.counter) {
      std::cerr << "ERROR: spill_participation counter mismatch: file=" << filename
                << " data_counter=" << data.counter
                << " meta_counter=" << meta.counter << std::endl;
      return false;
    }

    for (int i = 0; i < meta.nsources; ++i) {
      if (!out.add(meta.source_device[i], meta.source_fifo[i])) {
        std::cerr << "ERROR: spill_participation source capacity exceeded" << std::endl;
        return false;
      }
    }
    return true;
  }
};

struct output_t {
  TFile *file = nullptr;
  TTree *tree = nullptr;
  TTree *meta_tree = nullptr;
  data_t data;
  spill_participation_t meta;
  std::string tmpfilename;
  std::string finalfilename;
  Long64_t total_entries = 0;
  Long64_t total_meta_entries = 0;
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

    meta_tree = new TTree("spill_participation", "spill_participation");
    meta_tree->Branch("spill", &meta.spill, "spill/I");
    meta_tree->Branch("counter", &meta.counter, "counter/I");
    meta_tree->Branch("nsources", &meta.nsources, "nsources/I");
    meta_tree->Branch("source_device", meta.source_device, "source_device[nsources]/I");
    meta_tree->Branch("source_fifo", meta.source_fifo, "source_fifo[nsources]/I");
    return true;
  }

  void fill(const data_t &word)
  {
    data = word;
    tree->Fill();
  }

  void fill_participation(const spill_participation_t &word)
  {
    meta = word;
    meta_tree->Fill();
  }

  bool close_success()
  {
    if (!file)
      return true;

    auto entries = tree ? tree->GetEntries() : 0;
    auto meta_entries = meta_tree ? meta_tree->GetEntries() : 0;
    if (file->Write() <= 0) {
      std::cerr << "ERROR: failed writing output file: " << tmpfilename << std::endl;
      close_failed();
      return false;
    }

    file->Close();
    delete file;
    file = nullptr;
    tree = nullptr;
    meta_tree = nullptr;

    std::remove(finalfilename.c_str());
    if (std::rename(tmpfilename.c_str(), finalfilename.c_str()) != 0) {
      std::cerr << "ERROR: could not rename temporary output " << tmpfilename
                << " to " << finalfilename << std::endl;
      std::remove(tmpfilename.c_str());
      return false;
    }

    total_entries += entries;
    total_meta_entries += meta_entries;
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
    meta_tree = nullptr;
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
                << " counter=" << s.data.counter
                << " device=" << s.data.device
                << " fifo=" << s.data.fifo;
    if (s.meta_tree)
      std::cerr << " meta_entry=" << s.meta_entry << "/" << s.meta_entries;
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

  auto next_spill_counter = [&]() {
    int counter = std::numeric_limits<int>::max();
    bool found = false;
    for (auto &stream : streams) {
      if (!stream.valid)
        continue;
      if (!stream.data.is_start_spill()) {
        std::cerr << "ERROR: stream is not positioned at START_SPILL outside spill" << std::endl;
        print_streams(streams);
        return std::numeric_limits<int>::min();
      }
      if (!found || stream.data.counter < counter)
        counter = stream.data.counter;
      found = true;
    }
    return found ? counter : std::numeric_limits<int>::max();
  };

  auto best_stream = [&](const std::vector<int> &participants) {
    int ibest = -1;
    for (auto i : participants) {
      auto &stream = streams[i];
      if (!stream.valid || !stream.data.is_data_word())
        continue;
      if (ibest < 0 || stream.data.time < streams[ibest].data.time)
        ibest = i;
    }
    return ibest;
  };

  Long64_t nspill = 0;
  Long64_t collapsed_markers = 0;

  for (;;) {
    if (all_eof())
      break;

    int current_spill_counter = next_spill_counter();
    if (current_spill_counter == std::numeric_limits<int>::min())
      return fail("unexpected stream state before spill");
    if (current_spill_counter == std::numeric_limits<int>::max())
      return fail("no next spill counter found before EOF");

    std::vector<int> participants;
    data_t start_word;
    bool have_start_word = false;
    spill_participation_t participation;
    participation.clear(nspill, current_spill_counter);

    for (int i = 0; i < (int)streams.size(); ++i) {
      auto &stream = streams[i];
      if (!stream.valid)
        continue;
      if (!stream.data.is_start_spill())
        return fail("stream is not positioned at START_SPILL before selecting participants");
      if (stream.data.counter != current_spill_counter)
        continue;

      if (!have_start_word) {
        start_word = stream.data;
        have_start_word = true;
      }

      if (!stream.load_current_participation(participation))
        return fail("could not load spill participation metadata");
      participants.push_back(i);
    }

    if (participants.empty() || !have_start_word)
      return fail("no stream participates in selected spill");

    if (split_spills && !output.open(streams[0].tree, spill_filename(outfilename, nspill)))
      return fail("could not open split-spill output file");

    output.fill(start_word);
    output.fill_participation(participation);

    for (auto i : participants) {
      auto &stream = streams[i];
      ++stream.start_spills;
      if (!stream.next())
        return fail("failed while advancing past START_SPILL");
      if (!stream.valid)
        return fail("stream reached EOF immediately after START_SPILL");
    }

    for (;;) {
      int ibest = best_stream(participants);
      if (ibest >= 0) {
        output.fill(streams[ibest].data);
        if (!streams[ibest].next())
          return fail("failed while advancing data stream");
        continue;
      }

      bool all_participants_at_end = true;
      for (auto i : participants) {
        auto &stream = streams[i];
        if (!stream.valid)
          return fail("participating stream reached EOF before END_SPILL");
        if (!stream.data.is_end_spill()) {
          all_participants_at_end = false;
          if (stream.data.is_start_spill())
            return fail("participating stream reached START_SPILL before END_SPILL");
          return fail("unexpected word type inside spill: type=" + std::to_string(stream.data.type));
        }
        if (stream.data.counter != current_spill_counter)
          return fail("END_SPILL counter does not match current START_SPILL counter");
      }

      if (all_participants_at_end)
        break;
    }

    output.fill(streams[participants[0]].data);
    for (auto i : participants) {
      auto &stream = streams[i];
      ++stream.end_spills;
      if (!stream.next())
        return fail("failed while advancing past END_SPILL");
    }

    collapsed_markers += 2 * ((Long64_t)participants.size() - 1);
    ++nspill;

    if (split_spills && !output.close_success())
      return fail("could not close split-spill output file successfully");
  }

  for (auto &stream : streams) {
    if (stream.entry != stream.entries)
      return fail("not all input entries were consumed");
    if (stream.meta_tree && stream.meta_entry != stream.meta_entries)
      return fail("not all input spill_participation entries were consumed");
    if (stream.start_spills != stream.end_spills)
      return fail("START_SPILL/END_SPILL count mismatch within stream");
  }

  Long64_t expected_output = nsum - collapsed_markers;
  Long64_t actual_output = split_spills ? output.total_entries : output.tree->GetEntries();
  Long64_t actual_meta = split_spills ? output.total_meta_entries : output.meta_tree->GetEntries();
  if (actual_output != expected_output) {
    std::cerr << "ERROR: output entry-count mismatch" << std::endl;
    std::cerr << "Nsum=" << nsum << std::endl;
    std::cerr << "Ninput=" << streams.size() << std::endl;
    std::cerr << "Nspill=" << nspill << std::endl;
    std::cerr << "collapsed_markers=" << collapsed_markers << std::endl;
    std::cerr << "expected_output=" << expected_output << std::endl;
    std::cerr << "actual_output=" << actual_output << std::endl;
    return fail("output entry-count mismatch");
  }

  if (actual_meta != nspill) {
    std::cerr << "ERROR: spill_participation entry-count mismatch" << std::endl;
    std::cerr << "Nspill=" << nspill << std::endl;
    std::cerr << "metadata_entries=" << actual_meta << std::endl;
    return fail("spill_participation entry-count mismatch");
  }

  if (!split_spills && !output.close_success())
    return fail("could not close output file successfully");

  std::cout << " --- merge successful" << std::endl;
  std::cout << " --- merged spills: " << nspill << std::endl;
  if (split_spills)
    std::cout << " --- output files: " << output.nfiles << std::endl;
  std::cout << " --- input entries: " << nsum << std::endl;
  std::cout << " --- output entries: " << actual_output << std::endl;
  std::cout << " --- spill_participation entries: " << actual_meta << std::endl;
  for (int i = 0; i < (int)streams.size(); ++i) {
    auto &stream = streams[i];
    std::cout << "stream " << i << ": "
              << stream.entry << " / " << stream.entries << " consumed, "
              << "start_spill=" << stream.start_spills << ", "
              << "end_spill=" << stream.end_spills;
    if (stream.meta_tree)
      std::cout << ", meta=" << stream.meta_entry << " / " << stream.meta_entries;
    std::cout << std::endl;
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
