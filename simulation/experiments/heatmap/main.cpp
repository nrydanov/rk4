#include "forces.h"
#include "vdp_ensemble.hpp"
#include <CLI11.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct HeatmapResult {
  double delta1;
  double delta2;
  size_t eps_index;
  double L;
};

enum class CalcFailure { WrongArguments };

int main(int argc, char **argv) {
  CLI::App app{"Van der Pol Ensemble Simulation"};
  std::string config_path;
  std::string output_path;
  app.add_option("config", config_path)->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_path, "Output CSV file path")
      ->check(CLI::NonexistentPath | CLI::ExistingPath);
  CLI11_PARSE(app, argc, argv);

  YAML::Node config = YAML::LoadFile(config_path);

  int N = config["N"].as<size_t>();
  auto y0 = config["y0"].as<std::vector<double>>();
  auto lambdas = config["lambdas"].as<std::vector<double>>();
  auto coupling = config["coupling"].as<std::vector<double>>();
  auto adj = config["adj"].as<std::vector<std::vector<int>>>();

  if (y0.size() != static_cast<size_t>(2 * N)) {
    std::cerr << "y0 size must be 2 * N";
    return 1;
  }

  const auto f_node = config["force"];
  const auto pins = f_node["pins"].as<std::vector<int>>();
  const auto alpha = f_node["alpha"].as<std::vector<double>>();
  const auto fomg = f_node["fomg"].as<std::vector<double>>();

  const auto s_node = config["sim"];
  const double T = s_node["T"].as<double>();
  const double dt = s_node["dt"].as<double>();
  const double d_min = s_node["d_min"].as<double>();
  const double d_max = s_node["d_max"].as<double>();
  const double t_trans = s_node["t_transition"].as<double>();
  const auto epsilons = s_node["epsilons"].as<std::vector<double>>();
  const int grid_size = static_cast<int>(std::abs(d_max - d_min) / dt + 1);
  const size_t n_tasks = grid_size * grid_size * epsilons.size();

  auto opt_force = forces::SinForce::create(pins, alpha, fomg);
  if (!opt_force.has_value()) {
    std::cerr << "Got an error on constructing sin force";
    return 1;
  }
  auto sforce = opt_force.value();

  std::ofstream out(output_path);
  if (!out.is_open()) {
    std::cerr << "Got an error opening output file";
    return 1;
  }
  out << "delta1" << ",delta2" << ",eps" << ",L" << std::endl;

  std::vector<tl::expected<HeatmapResult, CalcFailure>> results(n_tasks);
  std::atomic<size_t> write_idx{0};

#pragma omp parallel for schedule(dynamic) collapse(3)
  for (size_t e = 0; e < epsilons.size(); ++e) {
    double eps = epsilons[e];
    auto coupling = std::vector<double>{eps, eps, eps};

    for (int i1 = 0; i1 < grid_size; ++i1) {
      double delta1 = d_min + i1 * dt;

      for (int i2 = 0; i2 < grid_size; ++i2) {
        double delta2 = d_min + i2 * dt;

        auto freqs = std::vector<double>{1.0, 1.0 + delta1, 1.0 + delta2};
        auto opt_solver =
            vdp_ensemble::VdPEnsembleSolver<forces::SinForce>::create(
                N, 0.0, y0, freqs, lambdas, coupling, adj, sforce);
        if (!opt_solver.has_value()) {
          results[write_idx.fetch_add(1)] =
              tl::unexpected(CalcFailure::WrongArguments);
          continue;
        }

        auto solver = opt_solver.value();

        for (double t = 0.0; t < t_trans; t += dt) {
          solver.step(dt);
        }

        float L = 0.0f;
        double acc = 0.0;
        for (double t = t_trans; t < T; t += dt) {
          solver.step(dt);
          auto state = solver.getState();
          double sum_x = state[0] + state[2] + state[4];
          acc += sum_x * sum_x * dt;
        }
        L = static_cast<float>(2.0 / (T - t_trans) * acc);

        results[write_idx.fetch_add(1)] = {delta1, delta2, e, L};
      }
    }
  }

  for (auto &r : results) {
    if (!r.has_value()) {
      out << "nan" << ",nan" << ",nan" << ",nan" << std::endl;
    } else {
      auto value = r.value();
      out << value.delta1 << "," << value.delta2 << ","
          << epsilons[value.eps_index] << "," << value.L << std::endl;
    }
  }

  return 0;
}
