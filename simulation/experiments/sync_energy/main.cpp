#include "config_array.hpp"
#include "forces.hpp"
#include "goals.hpp"
#include "provenance.hpp"
#include "vdp_ensemble.hpp"
#include <CLI11.hpp>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct HeatmapResult {
  double delta1, delta2, eps, L;
};

template <int N>
int run(const YAML::Node &config, const std::string &config_path,
        const std::string &output_path) {
  const auto y0 = to_array<2 * N>(config["y0"].as<std::vector<double>>());
  const auto lambdas = to_array<N>(config["lambdas"].as<std::vector<double>>());
  const auto adj = to_adj<N>(config["adj"].as<std::vector<std::vector<int>>>());

  const auto s = config["sim"];
  const double T = s["T"].as<double>();
  const double dt = s["dt"].as<double>();
  const double d_step = s["d_step"].as<double>();
  const double d_min = s["d_min"].as<double>();
  const double d_max = s["d_max"].as<double>();
  const double t_trans = s["t_transition"].as<double>();
  const auto epsilons = s["epsilons"].as<std::vector<double>>();
  const int grid_size = static_cast<int>(std::abs(d_max - d_min) / d_step + 1);
  const size_t n_tasks = (size_t)grid_size * grid_size * epsilons.size();

  std::ofstream out(output_path);
  if (!out.is_open()) {
    std::cerr << "Got an error opening output file";
    return 1;
  }
  provenance::write_header(out, "vdp_sync_energy", config_path, config);
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "delta1,delta2,eps,L\n";

  std::vector<HeatmapResult> results(n_tasks);

#pragma omp parallel for schedule(static) collapse(3)
  for (size_t e = 0; e < epsilons.size(); ++e) {
    for (int i1 = 0; i1 < grid_size; ++i1) {
      for (int i2 = 0; i2 < grid_size; ++i2) {
        const size_t idx =
            e * (size_t)grid_size * grid_size + i1 * grid_size + i2;
        const double eps = epsilons[e];
        const double delta1 = d_min + i1 * d_step;
        const double delta2 = d_min + i2 * d_step;

        std::array<double, N> freqs{1.0, 1.0 + delta1, 1.0 + delta2};
        std::array<double, N> eps_coupling;
        eps_coupling.fill(eps);

        using Solver = vdp_ensemble::VdPEnsembleSolver<N>;
        thread_local std::optional<Solver> solver;
        if (!solver)
          solver.emplace(y0, freqs, lambdas, eps_coupling, adj,
                         forces::NoopForce{});
        else
          solver->reset(y0, freqs, eps_coupling);

        for (double t = 0.0; t < t_trans; t += dt)
          solver->step(dt);
        double L_acc = 0.0;
        for (double t = t_trans; t < T; t += dt) {
          solver->step(dt);
          L_acc += goals::coherence(solver->getState().data(), N);
        }
        const double L = 2.0 / (T - t_trans) * L_acc * dt;
        results[idx] = {delta1, delta2, eps, L};
      }
    }
  }

  for (const auto &r : results)
    out << r.delta1 << "," << r.delta2 << "," << r.eps << "," << r.L << "\n";
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
  case 3:
    return run<3>(config, config_path, output_path);
  default:
    std::cerr << "sync_energy supports only N=3 (got " << N << ")\n";
    return 1;
  }
}
