#pragma once
#include "solver.h"
#include <tl/expected.hpp>

namespace vdp_ensemble {

enum class ConstructError { WrongArgSize };

constexpr const char *to_string(ConstructError e) {
  switch (e) {
  case ConstructError::WrongArgSize:
    return "Wrong argument size";
    return "Unknown error";
  }
}

template <class ForceFunc> class VdPEnsembleSolver : public RK4Solver {
private:
  int N;
  std::vector<double> omega2;
  std::vector<double> lambda;
  std::vector<double> coupling;
  std::vector<std::vector<int>> adj;
  ForceFunc forces;

protected:
  void derivs(double t, const std::vector<double> &state,
              std::vector<double> &dydx) override;

  VdPEnsembleSolver(double t0, const std::vector<double> &y0,
                    const std::vector<double> &freqs,
                    const std::vector<double> &lambda,
                    const std::vector<double> &coupling,
                    const std::vector<std::vector<int>> &adj, ForceFunc &func);

public:
  static tl::expected<VdPEnsembleSolver, ConstructError>
  create(size_t N, double t0, const std::vector<double> &y0,
         const std::vector<double> &freqs, const std::vector<double> &lambda,
         const std::vector<double> &coupling,
         const std::vector<std::vector<int>> &adj, ForceFunc &func) {
    if (N != lambda.size() || N != coupling.size() || N != adj.size() ||
        N != adj[0].size() || 2 * N != y0.size()) {
      return tl::make_unexpected(ConstructError::WrongArgSize);
    }
    return VdPEnsembleSolver(t0, y0, freqs, lambda, coupling, adj, func);
  }
};

template <class ForceFunc>
VdPEnsembleSolver<ForceFunc>::VdPEnsembleSolver(
    double t0, const std::vector<double> &y0, const std::vector<double> &freqs,
    const std::vector<double> &lambda, const std::vector<double> &coupling,
    const std::vector<std::vector<int>> &adj, ForceFunc &func)
    : RK4Solver(t0, y0), N(freqs.size()), lambda(lambda), coupling(coupling),
      adj(adj), forces(func) {
  omega2.resize(N);
  for (int i = 0; i < N; ++i) {
    omega2[i] = freqs[i] * freqs[i];
  }
}

template <class ForceFunc>
void VdPEnsembleSolver<ForceFunc>::derivs(double t,
                                          const std::vector<double> &state,
                                          std::vector<double> &dydx) {
  const auto &impacts = forces(t);

  for (int i = 0; i < N; ++i) {
    int idx_x = 2 * i;
    int idx_y = 2 * i + 1;
    double xi = state[idx_x];
    double yi = state[idx_y];

    double coupling_sum = 0.0;
    for (int j = 0; j < N; ++j) {
      if (adj[i][j] == 1) {
        coupling_sum += (state[2 * j] - xi);
      }
    }

    dydx[idx_x] = yi;
    dydx[idx_y] = (lambda[i] - xi * xi) * yi - omega2[i] * xi +
                  coupling[i] * coupling_sum;
  }

  for (const auto &impact : impacts) {
    if (impact.first < (int)dydx.size()) {
      dydx[impact.first] += impact.second;
    }
  }
}

} // namespace vdp_ensemble
