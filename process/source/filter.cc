#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

using filter_config_t = std::map<std::string, std::vector<std::string>>;

std::string trim(const std::string &input)
{
  std::size_t first = 0;
  while (first < input.size() && std::isspace(static_cast<unsigned char>(input[first])))
    ++first;
  std::size_t last = input.size();
  while (last > first && std::isspace(static_cast<unsigned char>(input[last - 1])))
    --last;
  return input.substr(first, last - first);
}

bool load_config(const std::string &filename, filter_config_t &config)
{
  std::ifstream stream(filename);
  if (!stream) {
    std::cerr << "ERROR: could not open filter configuration: " << filename << std::endl;
    return false;
  }

  std::string section, line;
  int line_number = 0;
  while (std::getline(stream, line)) {
    ++line_number;
    const auto comment = line.find('#');
    if (comment != std::string::npos)
      line.erase(comment);
    line = trim(line);
    if (line.empty())
      continue;

    if (line.front() == '[') {
      if (line.back() != ']' || line.size() <= 2) {
        std::cerr << "ERROR: malformed section at " << filename << ':' << line_number << std::endl;
        return false;
      }
      section = trim(line.substr(1, line.size() - 2));
      if (section.empty()) {
        std::cerr << "ERROR: empty section at " << filename << ':' << line_number << std::endl;
        return false;
      }
      config.try_emplace(section);
      continue;
    }

    if (section.empty()) {
      std::cerr << "ERROR: branch listed outside a section at " << filename << ':'
                << line_number << std::endl;
      return false;
    }
    for (const char character : line) {
      if (std::isspace(static_cast<unsigned char>(character))) {
        std::cerr << "ERROR: expected one branch name at " << filename << ':'
                  << line_number << std::endl;
        return false;
      }
    }
    auto &branches = config[section];
    if (std::find(branches.begin(), branches.end(), line) == branches.end())
      branches.push_back(line);
  }

  if (config.empty()) {
    std::cerr << "ERROR: filter configuration contains no trees" << std::endl;
    return false;
  }
  for (const auto &[tree, branches] : config) {
    if (branches.empty()) {
      std::cerr << "ERROR: tree section has no branches: [" << tree << ']' << std::endl;
      return false;
    }
  }
  return true;
}

bool copy_tree(TFile *input, TFile *output, const std::string &name,
               const std::vector<std::string> &branches)
{
  auto *tree = dynamic_cast<TTree *>(input->Get(name.c_str()));
  if (!tree) {
    std::cerr << "ERROR: requested tree is missing from input: " << name << std::endl;
    return false;
  }

  tree->SetBranchStatus("*", 0);
  for (const auto &branch : branches) {
    if (!tree->GetBranch(branch.c_str())) {
      std::cerr << "ERROR: requested branch is missing: tree=" << name
                << " branch=" << branch << std::endl;
      return false;
    }
    tree->SetBranchStatus(branch.c_str(), 1);
  }

  output->cd();
  auto *filtered = tree->CloneTree(-1, "fast");
  if (!filtered) {
    std::cerr << "ERROR: could not clone tree: " << name << std::endl;
    return false;
  }
  filtered->SetName(name.c_str());
  filtered->SetTitle(tree->GetTitle());
  std::cout << "tree " << name << ": " << tree->GetEntries()
            << " entries, " << branches.size() << " branches" << std::endl;
  return true;
}

} // namespace

int main(int argc, char **argv)
{
  namespace po = boost::program_options;
  std::string input, output, config_filename;
  po::options_description options("options");
  options.add_options()
    ("help,h", "show this help message")
    ("input,i", po::value<std::string>(&input)->required(), "input ROOT file")
    ("output,o", po::value<std::string>(&output)->required(), "output ROOT file")
    ("config,c", po::value<std::string>(&config_filename)->required(),
     "filter configuration file");

  try {
    po::variables_map variables;
    po::store(po::parse_command_line(argc, argv, options), variables);
    if (variables.count("help")) {
      std::cout << options << std::endl;
      return 0;
    }
    po::notify(variables);
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n' << options << std::endl;
    return 1;
  }

  filter_config_t config;
  if (!load_config(config_filename, config))
    return 1;

  std::unique_ptr<TFile> input_file(TFile::Open(input.c_str(), "READ"));
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "ERROR: could not open input ROOT file: " << input << std::endl;
    return 1;
  }
  std::unique_ptr<TFile> output_file(TFile::Open(output.c_str(), "RECREATE"));
  if (!output_file || output_file->IsZombie()) {
    std::cerr << "ERROR: could not create output ROOT file: " << output << std::endl;
    return 1;
  }
  for (const auto &[tree, branches] : config) {
    if (!copy_tree(input_file.get(), output_file.get(), tree, branches))
      return 1;
  }
  output_file->Write();
  std::cout << "output: " << output << std::endl;
  return 0;
}
