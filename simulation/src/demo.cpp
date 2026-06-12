#include "config_array.hpp"
#include "forces.hpp"
#include "provenance.hpp"
#include "vdp_ensemble.hpp"
#include <CLI11.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

template <int N>
int run(const YAML::Node &config, const std::string &config_path,
        const std::string &output_path) {
  const auto y0 = to_array<2 * N>(config["y0"].as<std::vector<double>>());
  const auto freqs = to_array<N>(config["freqs"].as<std::vector<double>>());
  const auto lambdas = to_array<N>(config["lambdas"].as<std::vector<double>>());
  const auto coupling =
      to_array<N>(config["coupling"].as<std::vector<double>>());
  const auto adj = to_adj<N>(config["adj"].as<std::vector<std::vector<int>>>());

  auto f = config["force"];
  auto pins = f["pins"].as<std::vector<int>>();
  auto alpha = f["alpha"].as<std::vector<double>>();
  auto fomg = f["fomg"].as<std::vector<double>>();

  auto s = config["sim"];
  double T = s["T"].as<double>();
  double dt = s["dt"].as<double>();

  auto opt_force = forces::SinForce::create(pins, alpha, fomg);
  if (!opt_force.has_value()) {
    std::cerr << "Got an error on constructing sin force";
    return 1;
  }

  vdp_ensemble::VdPEnsembleSolver<N, forces::SinForce> solver(
      y0, freqs, lambdas, coupling, adj, opt_force.value());

  std::ofstream out(output_path);
  if (!out.is_open()) {
    std::cerr << "Got an error opening output file";
    return 1;
  }

  provenance::write_header(out, "vdp_sim", config_path, config);
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
  return 0;
}

int main(int argc, char **argv) {
  CLI::App app{"Van der Pol Ensemble Simulation"};
  std::string config_path, output_path;
  app.add_option("config", config_path)->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_path, "Output CSV file path")
      ->check(CLI::NonexistentPath | CLI::ExistingPath);
  CLI11_PARSE(app, argc, argv);

  YAML::Node config = YAML::LoadFile(config_path);
  const int N = config["N"].as<int>();
  switch (N) {
  case 1:
    return run<1>(config, config_path, output_path);
  case 2:
    return run<2>(config, config_path, output_path);
  case 3:
    return run<3>(config, config_path, output_path);
  case 4:
    return run<4>(config, config_path, output_path);
  case 5:
    return run<5>(config, config_path, output_path);
  case 6:
    return run<6>(config, config_path, output_path);
  case 7:
    return run<7>(config, config_path, output_path);
  case 8:
    return run<8>(config, config_path, output_path);
  default:
    std::cerr << "demo supports N in 1..8 (got " << N << ")\n";
    return 1;
  }
}
