#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include <iostream>
#include <string>


bool
valid_fifo_column(const data_t &data)
{
  if (!data.is_alcor_hit())
    return true;

  int first_column = 2 * (data.fifo % 4);
  return data.column == first_column || data.column == first_column + 1;
}

bool
keep_data(const data_t &data)
{
  return valid_fifo_column(data);
}

bool
cleaner(const std::string &filename, const std::string &outfilename)
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
  auto nev = tin->GetEntries();

  data_t data;
  data.link_to_tree(tin);

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    fin->Close();
    return false;
  }
  auto tout = tin->CloneTree(0);

  Long64_t nkept = 0;
  Long64_t ndropped = 0;
  Long64_t ndropped_column = 0;

  for (Long64_t iev = 0; iev < nev; ++iev) {
    auto bytes = tin->GetEntry(iev);
    if (bytes <= 0) {
      std::cerr << "ERROR: ROOT GetEntry failed for entry " << iev << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    if (!keep_data(data)) {
      ++ndropped;
      if (data.is_alcor_hit() && !valid_fifo_column(data))
        ++ndropped_column;
      continue;
    }

    tout->Fill();
    ++nkept;
  }

  auto written = tout->GetEntries();
  if (written != nkept) {
    std::cerr << "ERROR: output entry count mismatch: kept=" << nkept
              << " output=" << written << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "input entries:    " << nev << std::endl;
  std::cout << "kept entries:     " << nkept << std::endl;
  std::cout << "dropped entries:  " << ndropped << std::endl;
  std::cout << "dropped columns:  " << ndropped_column << std::endl;
  std::cout << "output entries:   " << written << std::endl;

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
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file");

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

  return cleaner(input, output) ? 0 : 1;
}
