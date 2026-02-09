#include "forces.h"
#include "vdp_ensemble.hpp"
#include <CLI11.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

int main(int argc, char **argv) {
  CLI::App app{"Van der Pol Ensemble Simulation"};
  std::string config_path;
  std::string output_path;
  app.add_option("config", config_path)->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_path, "Output CSV file path")
      ->check(CLI::NonexistentPath | CLI::ExistingPath);
  CLI11_PARSE(app, argc, argv);

  try {
    YAML::Node config = YAML::LoadFile(config_path);

    int N = config["N"].as<int>();
    auto y0 = config["y0"].as<std::vector<double>>();
    auto freqs = config["freqs"].as<std::vector<double>>();
    auto lambdas = config["lambdas"].as<std::vector<double>>();
    auto coupling = config["coupling"].as<std::vector<double>>();
    auto adj = config["adj"].as<std::vector<std::vector<int>>>();

    if (y0.size() != static_cast<size_t>(2 * N)) {
      throw std::runtime_error("y0 size must be 2*N");
    }

    auto f_node = config["force"];
    auto pins = f_node["pins"].as<std::vector<int>>();
    auto alpha = f_node["alpha"].as<std::vector<double>>();
    auto fomg = f_node["fomg"].as<std::vector<double>>();

    auto s_node = config["sim"];
    double T = s_node["T"].as<double>();
    double dt = s_node["dt"].as<double>();

    SinForce sforce(pins, alpha, fomg);
    VdPEnsembleSolver<SinForce> solver(0.0, y0, freqs, lambdas, coupling, adj,
                                       sforce);

    std::ofstream out(output_path);
    if (!out.is_open())
      throw std::runtime_error("File error");

    out << "t";
    for (int i = 0; i < N; ++i)
      out << ",x" << i << ",y" << i;
    out << "\n";

    out.precision(8);
    out << std::fixed;

    for (double t = 0.0; t < T; t += dt) {
      out << solver.getTime();
      for (double val : solver.getState())
        out << "," << val;
      out << "\n";
      solver.step(dt);
    }

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
