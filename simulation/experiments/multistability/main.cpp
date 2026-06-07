#include "coupling.hpp"
#include "forces.hpp"
#include "goals.hpp"
#include "phase_diff_slopes.hpp"
#include "vdp_ensemble.hpp"
#include <CLI11.hpp>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct Config {
  int N;
  std::vector<double> lambdas;
  std::vector<std::vector<int>> adj;
  double T, dt, d_min, d_step, t_trans;
  std::vector<double> epsilons;
  int grid_size;
  int n_ic;
  double ic_range;
  uint32_t seed;
};

struct Result {
  double delta1, delta2, eps;
  int coupling_type;
  double x0, y0, x1, y1, x2, y2;
  double L, A, P;
  double s01, s02, s12;
};

void print_progress(size_t done, size_t total,
                    std::chrono::steady_clock::time_point start,
                    const std::string &label) {
  double frac = (double)done / total;
  int bar_width = 35;
  int filled = (int)(frac * bar_width);
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - start).count();
  long eta = done > 0 ? (long)(elapsed * (double)(total - done) / done) : 0;

  std::cerr << "\r" << label << " [";
  for (int i = 0; i < bar_width; ++i)
    std::cerr << (i < filled ? '=' : (i == filled ? '>' : ' '));
  std::cerr << "] " << std::fixed << std::setprecision(1) << (frac * 100)
            << "% (" << done << "/" << total << ")";
  if (done > 0 && done < total)
    std::cerr << " ETA: " << eta / 60 << "m" << eta % 60 << "s";
  std::cerr << "   " << std::flush;
}

template <class CouplingFunc>
void sweep(int coupling_type_id, const std::string &label, const Config &cfg,
           std::vector<Result> &out) {
  const size_t n_tasks =
      (size_t)cfg.grid_size * cfg.grid_size * cfg.epsilons.size() * cfg.n_ic;
  std::vector<Result> local(n_tasks);
  std::atomic<size_t> idx{0};
  const size_t report_every = std::max((size_t)1, n_tasks / 50);
  auto start = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic) collapse(3)
  for (size_t e = 0; e < cfg.epsilons.size(); ++e) {
    for (int i1 = 0; i1 < cfg.grid_size; ++i1) {
      for (int i2 = 0; i2 < cfg.grid_size; ++i2) {
        double eps = cfg.epsilons[e];
        double delta1 = cfg.d_min + i1 * cfg.d_step;
        double delta2 = cfg.d_min + i2 * cfg.d_step;
        auto freqs = std::vector<double>{1.0, 1.0 + delta1, 1.0 + delta2};
        auto eps_coupling = std::vector<double>(cfg.N, eps);

        std::mt19937 rng(cfg.seed ^
                         std::hash<size_t>{}(e * cfg.grid_size * cfg.grid_size +
                                             i1 * cfg.grid_size + i2));
        std::uniform_real_distribution<double> dist(-cfg.ic_range, cfg.ic_range);

        for (int ic = 0; ic < cfg.n_ic; ++ic) {
          std::vector<double> y0(2 * cfg.N);
          for (auto &v : y0) v = dist(rng);

          forces::NoopForce noop;
          CouplingFunc cf{};
          auto opt_solver = vdp_ensemble::VdPEnsembleSolver<
              forces::NoopForce, CouplingFunc>::create(cfg.N, 0.0, y0, freqs,
                                                       cfg.lambdas,
                                                       eps_coupling, cfg.adj,
                                                       noop, cf);
          if (!opt_solver.has_value()) continue;
          auto solver = opt_solver.value();

          for (double t = 0.0; t < cfg.t_trans; t += cfg.dt) solver.step(cfg.dt);

          phase_diff_slopes<double> pds(cfg.N, 2, 0, 1);
          ampl_t<double, int> ampl_goal(2, cfg.N);
          phase_t<double, int> phase_goal(2, cfg.N);

          double L_acc = 0.0, A_acc = 0.0, P_acc = 0.0;
          int steps = 0;
          for (double t = cfg.t_trans; t < cfg.T; t += cfg.dt, ++steps) {
            solver.step(cfg.dt);
            const auto &state = solver.getState();
            pds.push(t, state.data());
            double sum_x = 0.0;
            for (int i = 0; i < cfg.N; ++i) sum_x += state[2 * i];
            L_acc += sum_x * sum_x;
            A_acc += ampl_goal(state.data());
            P_acc += phase_goal(state.data());
          }
          double L = 2.0 / (cfg.T - cfg.t_trans) * L_acc * cfg.dt;
          double A = A_acc / steps;
          double P = P_acc / steps;

          size_t pos = idx.fetch_add(1);
          local[pos] = {delta1,    delta2,    eps,   coupling_type_id,
                        y0[0],     y0[1],     y0[2], y0[3],
                        y0[4],     y0[5],     L,     A,
                        P,         pds(0, 1), pds(0, 2), pds(1, 2)};

          if (pos % report_every == 0) {
#pragma omp critical
            print_progress(idx.load(), n_tasks, start, label);
          }
        }
      }
    }
  }

  print_progress(n_tasks, n_tasks, start, label);
  std::cerr << "\n";
  size_t n = idx.load();
  out.insert(out.end(), local.begin(), local.begin() + n);
}

int main(int argc, char **argv) {
  CLI::App app{"Van der Pol Multistability Experiment"};
  std::string config_path, output_path;
  app.add_option("config", config_path)->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_path)
      ->check(CLI::NonexistentPath | CLI::ExistingPath);
  CLI11_PARSE(app, argc, argv);

  YAML::Node yaml = YAML::LoadFile(config_path);
  const auto &s = yaml["sim"];
  double d_min = s["d_min"].as<double>();
  double d_max = s["d_max"].as<double>();
  double d_step = s["d_step"].as<double>();

  Config cfg{
      .N = yaml["N"].as<int>(),
      .lambdas = yaml["lambdas"].as<std::vector<double>>(),
      .adj = yaml["adj"].as<std::vector<std::vector<int>>>(),
      .T = s["T"].as<double>(),
      .dt = s["dt"].as<double>(),
      .d_min = d_min,
      .d_step = d_step,
      .t_trans = s["t_transition"].as<double>(),
      .epsilons = s["epsilons"].as<std::vector<double>>(),
      .grid_size = static_cast<int>(std::abs(d_max - d_min) / d_step + 1),
      .n_ic = s["n_ic"].as<int>(),
      .ic_range = s["ic_range"].as<double>(),
      .seed = s["seed"].as<uint32_t>(),
  };

  std::cerr << "Grid: " << cfg.grid_size << "x" << cfg.grid_size
            << "  epsilons: " << cfg.epsilons.size()
            << "  n_ic: " << cfg.n_ic << "\n";

  std::vector<Result> results;
  sweep<coupling::Inertial>(0, "Inertial        ", cfg, results);
  sweep<coupling::InertialNorm>(1, "InertialNorm    ", cfg, results);
  sweep<coupling::Dissipative>(2, "Dissipative     ", cfg, results);
  sweep<coupling::DissipativeNorm>(3, "DissipativeNorm ", cfg, results);

  std::ofstream out(output_path);
  if (!out.is_open()) {
    std::cerr << "Failed to open output file\n";
    return 1;
  }
  out << "delta1,delta2,eps,coupling_type,x0,y0,x1,y1,x2,y2,L,A,P,s01,s02,s12\n";
  for (auto &r : results) {
    out << r.delta1 << "," << r.delta2 << "," << r.eps << "," << r.coupling_type
        << "," << r.x0 << "," << r.y0 << "," << r.x1 << "," << r.y1 << ","
        << r.x2 << "," << r.y2 << "," << r.L << "," << r.A << "," << r.P
        << "," << r.s01 << "," << r.s02 << "," << r.s12 << "\n";
  }

  return 0;
}
