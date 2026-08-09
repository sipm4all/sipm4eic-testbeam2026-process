#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include <iostream>
#include <string>
#include <vector>


void
sorter(const std::string filename, const std::string outfilename, int window)
{
  auto fin = TFile::Open(filename.c_str());
  if (!fin || fin->IsZombie()) {
    std::cerr << " --- could not open input file: " << filename << std::endl;
    return;
  }

  auto tin = (TTree *)fin->Get("alcor");
  if (!tin) {
    std::cerr << " --- could not find 'alcor' tree in input file" << std::endl;
    fin->Close();
    return;
  }
  auto nev = tin->GetEntries();

  data_t data;
  data.link_to_tree(tin);

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << " --- could not create output file: " << outfilename << std::endl;
    fin->Close();
    return;
  }
  auto tout = tin->CloneTree(0);

  std::vector<data_t> buffer;
  double max_time = 0.;

  auto fill_data = [&]() {
    tout->Fill();
  };

  auto flush_buffer = [&]() {
    while (!buffer.empty()) {
      int ibest = 0;
      for (int i = 1; i < (int)buffer.size(); ++i)
        if (buffer[i].time < buffer[ibest].time)
          ibest = i;
      data = buffer[ibest];
      fill_data();
      buffer.erase(buffer.begin() + ibest);
    }
  };

  auto release_buffer = [&](double clock) {
    auto &buf = buffer;
    for (;;) {
      int ibest = -1;
      for (int i = 0; i < (int)buf.size(); ++i) {
        if (buf[i].time > clock - window)
          continue;
        if (ibest < 0 || buf[i].time < buf[ibest].time)
          ibest = i;
      }
      if (ibest < 0)
        break;
      data = buf[ibest];
      fill_data();
      buf.erase(buf.begin() + ibest);
    }
  };

  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);

    /** start of spill is detected **/
    if (data.is_start_spill()) {
      auto spill = data;
      flush_buffer();
      data = spill;
      fill_data();
      /** do any needed reset action **/
      max_time = 0.;
      continue;
    }

    /** end of spill is detected **/
    if (data.is_end_spill()) {
      auto spill = data;
      flush_buffer();
      data = spill;
      fill_data();
      /** do any needed reset action **/
      max_time = 0.;
      continue;
    }

    /** not an ALCOR hit and not a trigger tag **/
    if (!data.is_alcor_hit() && !data.is_trigger_tag()) {
      fill_data();
      continue;
    }

    /** compute time if no calibrated time branch is available **/
    if (!data.has_time) {
      data.set_nominal_time();
    }

    if (buffer.empty() || data.time > max_time)
      max_time = data.time;

    buffer.push_back(data);
    release_buffer(max_time);
  }

  flush_buffer();

  fout->Write();
  fout->Close();
  fin->Close();
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string output;
  int window;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("window,w", po::value<int>(&window)->required(), "reorder window in clock cycles");

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

  sorter(input, output, window);
  return 0;
}
