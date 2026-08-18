#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"
#include "geometry.h"

#include <iostream>
#include <string>

bool
coordinator(const std::string &filename,
            const std::string &outfilename)
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

  double x = 0.;
  double y = 0.;
  bool input_has_x = tin->GetBranch("x") != nullptr;
  bool input_has_y = tin->GetBranch("y") != nullptr;
  if (input_has_x)
    tin->SetBranchAddress("x", &x);
  if (input_has_y)
    tin->SetBranchAddress("y", &y);

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    fin->Close();
    return false;
  }

  auto tout = tin->CloneTree(0);
  if (!input_has_x)
    tout->Branch("x", &x, "x/D");
  if (!input_has_y)
    tout->Branch("y", &y, "y/D");

  Long64_t nalcor = 0;
  Long64_t ntiming = 0;
  Long64_t ncherenkov = 0;
  Long64_t nunmapped = 0;
  Long64_t ncontrol = 0;
  Long64_t nother = 0;

  for (Long64_t iev = 0; iev < nev; ++iev) {
    auto bytes = tin->GetEntry(iev);
    if (bytes <= 0) {
      std::cerr << "ERROR: ROOT GetEntry failed for entry " << iev << std::endl;
      fout->Close();
      fin->Close();
      return false;
    }

    x = 0.;
    y = 0.;

    if (data.is_alcor_hit()) {
      ++nalcor;
      if (data.device == 200)
        ++ntiming;
      else
        ++ncherenkov;

      if (!geometry_t::coordinate(data, x, y)) {
        ++nunmapped;
        if (nunmapped <= 20) {
          std::cerr << "WARNING: could not assign coordinates"
                    << " entry=" << iev
                    << " device=" << data.device
                    << " fifo=" << data.fifo
                    << " column=" << data.column
                    << " pixel=" << data.pixel
                    << std::endl;
        }
      }
    } else {
      if (data.is_start_spill() || data.is_end_spill())
        ++ncontrol;
      else
        ++nother;
    }

    tout->Fill();
  }

  auto written = tout->GetEntries();
  if (written != nev) {
    std::cerr << "ERROR: output entry count mismatch: input=" << nev
              << " output=" << written << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  if (fout->Write() <= 0) {
    std::cerr << "ERROR: failed writing output file: " << outfilename << std::endl;
    fout->Close();
    fin->Close();
    return false;
  }

  fout->Close();
  fin->Close();

  std::cout << "input entries:             " << nev << std::endl;
  std::cout << "ALCOR hits coordinated:    " << nalcor << std::endl;
  std::cout << "TIMING hits:               " << ntiming << std::endl;
  std::cout << "Cherenkov hits:            " << ncherenkov << std::endl;
  std::cout << "unmapped ALCOR hits:       " << nunmapped << std::endl;
  std::cout << "control words preserved:   " << ncontrol << std::endl;
  std::cout << "other words preserved:     " << nother << std::endl;
  std::cout << "output entries:            " << written << std::endl;

  if (nunmapped > 0)
    std::cerr << "WARNING: " << nunmapped
              << " ALCOR hits were written with x=0 y=0 because no geometry mapping was found"
              << std::endl;

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
    std::cerr << "ERROR: " << e.what() << std::endl;
    std::cerr << options << std::endl;
    return 1;
  }

  return coordinator(input, output) ? 0 : 1;
}
