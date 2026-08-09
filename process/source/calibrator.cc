#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include "data_word.h"

#include "calibration.h"

#include <iostream>
#include <string>


bool
calibrator(const std::string &filename,
           const std::string &outfilename,
           const std::string &config)
{
  calibration_t calibration;
  if (!calibration.load(config))
    return false;

  std::cout << "calibration file: " << config << std::endl;
  std::cout << "TDC rules:      " << calibration.tdc_rules() << std::endl;
  std::cout << "CHANNEL rules:  " << calibration.channel_rules() << std::endl;
  std::cout << "TRIGGER rules:  " << calibration.trigger_rules() << std::endl;

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
  bool input_has_time = tin->GetBranch("time") != nullptr;

  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create output file: " << outfilename << std::endl;
    fin->Close();
    return false;
  }

  auto tout = tin->CloneTree(0);
  if (!input_has_time)
    tout->Branch("time", &data.time, "time/D");

  Long64_t nalcor = 0;
  Long64_t ntrigger = 0;
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

    if (data.is_alcor_hit()) {
      if (data.tdc < 0 || data.tdc >= 4) {
        std::cerr << "ERROR: invalid TDC index " << data.tdc
                  << " at input entry " << iev << std::endl;
        fout->Close();
        fin->Close();
        return false;
      }

      tdc_calib_t tdc;
      channel_calib_t channel;
      std::string error;
      if (!calibration.tdc(data.device, data.fifo, data.column, data.pixel, data.tdc, tdc, error)) {
        std::cerr << "ERROR: " << error << " at input entry " << iev << std::endl;
        fout->Close();
        fin->Close();
        return false;
      }
      if (!calibration.channel(data.device, data.fifo, data.column, data.pixel, channel, error)) {
        std::cerr << "ERROR: " << error << " at input entry " << iev << std::endl;
        fout->Close();
        fin->Close();
        return false;
      }

      double phase = tdc.off + tdc.iif * data.fine;
      data.time = data.coarse + data_t::rollover_to_clock * data.rollover;
      data.time -= phase;
      data.time -= channel.offset;
      ++nalcor;
    } else if (data.is_trigger_tag()) {
      trigger_calib_t trigger;
      std::string error;
      if (!calibration.trigger(data.device, data.fifo, trigger, error)) {
        std::cerr << "ERROR: " << error << " at input entry " << iev << std::endl;
        fout->Close();
        fin->Close();
        return false;
      }

      data.time = data.coarse + data_t::rollover_to_clock * data.rollover;
      data.time -= trigger.offset;
      ++ntrigger;
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
  std::cout << "ALCOR hits calibrated:     " << nalcor << std::endl;
  std::cout << "trigger words calibrated:  " << ntrigger << std::endl;
  std::cout << "control words preserved:   " << ncontrol << std::endl;
  std::cout << "other words preserved:     " << nother << std::endl;
  std::cout << "output entries:            " << written << std::endl;
  std::cout << "cached TDC calibrations:   " << calibration.cached_tdc() << std::endl;
  std::cout << "cached channel offsets:    " << calibration.cached_channels() << std::endl;
  std::cout << "cached trigger offsets:    " << calibration.cached_triggers() << std::endl;

  return true;
}

int
main(int argc, char **argv)
{
  namespace po = boost::program_options;

  std::string input;
  std::string output;
  std::string config;

  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("config,c", po::value<std::string>(&config)->required(), "timing calibration file");

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

  return calibrator(input, output, config) ? 0 : 1;
}
