#include "config_array.hpp"
#include "coupling.hpp"
#include "forces.hpp"
#include "goals.hpp"
#include "phase_diff_slopes.hpp"
#include "provenance.hpp"
#include "vdp_ensemble.hpp"
#include <CLI11.hpp>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <optional>
#include <random>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

template <int N> struct Config {
  std::array<double, N> lambdas;
  std::array<std::array<int, N>, N> adj;
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
  double x0, y0, x1, y1, x2, y2;       // начальные условия
  double xf0, yf0, xf1, yf1, xf2, yf2; // конечное состояние
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
                     std::chrono::steady_clock::now() - start)
                     .count();
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

template <int N, class CouplingFunc>
void sweep(int coupling_type_id, const std::string &label, const Config<N> &cfg,
           std::vector<Result> &out) {
  const size_t n_tasks =
      (size_t)cfg.grid_size * cfg.grid_size * cfg.epsilons.size() * cfg.n_ic;
  std::vector<Result> local(n_tasks);
  std::atomic<size_t> progress{0};
  const size_t report_every = std::max((size_t)1, n_tasks / 50);
  auto start = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(static) collapse(3)
  for (size_t e = 0; e < cfg.epsilons.size(); ++e) {
    for (int i1 = 0; i1 < cfg.grid_size; ++i1) {
      for (int i2 = 0; i2 < cfg.grid_size; ++i2) {
        const size_t base = (e * (size_t)cfg.grid_size * cfg.grid_size +
                             i1 * cfg.grid_size + i2) *
                            cfg.n_ic;
        const double delta1 = cfg.d_min + i1 * cfg.d_step;
        const double delta2 = cfg.d_min + i2 * cfg.d_step;

        std::array<double, N> freqs{1.0, 1.0 + delta1, 1.0 + delta2};
        std::array<double, N> eps_coupling;
        eps_coupling.fill(cfg.epsilons[e]);

        std::mt19937 rng(cfg.seed ^
                         std::hash<size_t>{}(e * cfg.grid_size * cfg.grid_size +
                                             i1 * cfg.grid_size + i2));
        std::uniform_real_distribution<double> dist(-cfg.ic_range,
                                                    cfg.ic_range);

        using Solver =
            vdp_ensemble::VdPEnsembleSolver<N, forces::NoopForce, CouplingFunc>;
        thread_local std::optional<Solver> solver;
        thread_local std::optional<phase_diff_slopes<double>> pds;
        thread_local std::optional<phase_t<double, int>> phase_goal;
        ampl_t<double, int> ampl_goal(2, N);

        std::array<double, 2 * N> y0;

        for (int ic = 0; ic < cfg.n_ic; ++ic) {
          for (auto &v : y0)
            v = dist(rng);

          if (!solver)
            solver.emplace(y0, freqs, cfg.lambdas, eps_coupling, cfg.adj,
                           forces::NoopForce{});
          else
            solver->reset(y0, freqs, eps_coupling);

          for (double t = 0.0; t < cfg.t_trans; t += cfg.dt)
            solver->step(cfg.dt);

          if (!pds)
            pds.emplace(N, 2, 0, 1);
          else
            pds->reset();
          if (!phase_goal)
            phase_goal.emplace(2, N);

          double L_acc = 0.0, A_acc = 0.0, P_acc = 0.0;
          int steps = 0;
          for (double t = cfg.t_trans; t < cfg.T; t += cfg.dt, ++steps) {
            solver->step(cfg.dt);
            const auto &state = solver->getState();
            std::array<double, N> raw_phase;
            for (int i = 0; i < N; ++i)
              raw_phase[i] = std::atan2(state[2 * i + 1], state[2 * i]);
            pds->push_raw(t, raw_phase.data());
            double sum_x = 0.0;
            for (int i = 0; i < N; ++i)
              sum_x += state[2 * i];
            L_acc += sum_x * sum_x;
            A_acc += ampl_goal(state.data());
            P_acc += phase_goal->from_raw_phases(raw_phase.data());
          }
          double L = 2.0 / (cfg.T - cfg.t_trans) * L_acc * cfg.dt;
          double A = A_acc / steps;
          double P = P_acc / steps;

          const auto &yf = solver->getState();
          local[base + ic] = {delta1,
                              delta2,
                              cfg.epsilons[e],
                              coupling_type_id,
                              y0[0],
                              y0[1],
                              y0[2],
                              y0[3],
                              y0[4],
                              y0[5],
                              yf[0],
                              yf[1],
                              yf[2],
                              yf[3],
                              yf[4],
                              yf[5],
                              L,
                              A,
                              P,
                              (*pds)(0, 1),
                              (*pds)(0, 2),
                              (*pds)(1, 2)};

          size_t done = progress.fetch_add(1, std::memory_order_relaxed) + 1;
          if (done % report_every == 0)
            print_progress(done, n_tasks, start, label);
        }
      }
    }
  }

  print_progress(n_tasks, n_tasks, start, label);
  std::cerr << "\n";
  out.insert(out.end(), local.begin(), local.end());
}

template <int N>
int run(const YAML::Node &yaml, const std::string &config_path,
        const std::string &output_path) {
  const auto &s = yaml["sim"];
  double d_min = s["d_min"].as<double>();
  double d_max = s["d_max"].as<double>();
  double d_step = s["d_step"].as<double>();

  Config<N> cfg{
      .lambdas = to_array<N>(yaml["lambdas"].as<std::vector<double>>()),
      .adj = to_adj<N>(yaml["adj"].as<std::vector<std::vector<int>>>()),
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
            << "  epsilons: " << cfg.epsilons.size() << "  n_ic: " << cfg.n_ic
            << "  OpenMP threads: " << omp_get_max_threads() << "\n";

  std::vector<Result> results;
  sweep<N, coupling::Inertial>(0, "Inertial        ", cfg, results);
  sweep<N, coupling::InertialNorm>(1, "InertialNorm    ", cfg, results);
  sweep<N, coupling::Dissipative>(2, "Dissipative     ", cfg, results);
  sweep<N, coupling::DissipativeNorm>(3, "DissipativeNorm ", cfg, results);

  std::ofstream out(output_path);
  if (!out.is_open()) {
    std::cerr << "Failed to open output file\n";
    return 1;
  }
  provenance::write_header(out, "vdp_multistability", config_path, yaml);
  // max_digits10 — чтобы по записанным НУ можно было в точности повторить
  // отдельную траекторию
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "delta1,delta2,eps,coupling_type,x0,y0,x1,y1,x2,y2,"
         "xf0,yf0,xf1,yf1,xf2,yf2,L,A,P,s01,s02,s12\n";
  for (auto &r : results) {
    out << r.delta1 << "," << r.delta2 << "," << r.eps << "," << r.coupling_type
        << "," << r.x0 << "," << r.y0 << "," << r.x1 << "," << r.y1 << ","
        << r.x2 << "," << r.y2 << "," << r.xf0 << "," << r.yf0 << "," << r.xf1
        << "," << r.yf1 << "," << r.xf2 << "," << r.yf2 << "," << r.L << ","
        << r.A << "," << r.P << "," << r.s01 << "," << r.s02 << "," << r.s12
        << "\n";
  }
  return 0;
}

int main(int argc, char **argv) {
  CLI::App app{"Van der Pol Multistability Experiment"};
  std::string config_path, output_path;
  app.add_option("config", config_path)->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_path)
      ->check(CLI::NonexistentPath | CLI::ExistingPath);
  CLI11_PARSE(app, argc, argv);

  YAML::Node yaml = YAML::LoadFile(config_path);
  const int N = yaml["N"].as<int>();
  switch (N) {
  case 3:
    return run<3>(yaml, config_path, output_path);
  default:
    std::cerr << "multistability supports only N=3 (got " << N << ")\n";
    return 1;
  }
}
