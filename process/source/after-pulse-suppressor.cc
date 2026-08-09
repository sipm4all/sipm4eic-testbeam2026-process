#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include <array>
#include <iostream>
#include <string>


void
after_pulse_suppressor(const std::string filename, const std::string outfilename, double window)
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

  std::array<bool, 32> has_previous_hit = {false};
  std::array<double, 32> previous_time = {0.};

  auto reset = [&]() {
    has_previous_hit.fill(false);
  };

  for (int iev = 0; iev < nev; ++iev) {
    tin->GetEntry(iev);

    /** start of spill is detected **/
    if (data.is_start_spill()) {
      tout->Fill();
      /** do any needed reset action **/
      reset();
      continue;
    }

    /** end of spill is detected **/
    if (data.is_end_spill()) {
      tout->Fill();
      /** do any needed reset action **/
      reset();
      continue;
    }

    /** trigger tag is not considered for afterpulse suppression **/
    if (data.is_trigger_tag()) {
      tout->Fill();
      continue;
    }

    /** not an ALCOR hit **/
    if (!data.is_alcor_hit()) {
      tout->Fill();
      continue;
    }

    /** compute time if no calibrated time branch is available **/
    if (!data.has_time)
      data.set_nominal_time();

    auto channel = data.channel();
    if (channel < 0 || channel >= (int)has_previous_hit.size()) {
      tout->Fill();
      continue;
    }

    if (has_previous_hit[channel] && data.time - previous_time[channel] < window)
      continue;

    tout->Fill();
    previous_time[channel] = data.time;
    has_previous_hit[channel] = true;
  }

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
  double window;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("window,w", po::value<double>(&window)->required(), "afterpulse suppression window in clock cycles");

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

  after_pulse_suppressor(input, output, window);
  return 0;
}
