#pragma once
#include "coupling.hpp"
#include "forces.hpp"
#include "solver.hpp"
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

template <class ForceFunc = forces::NoopForce, class CouplingFunc = coupling::Inertial>
class VdPEnsembleSolver : public RK4Solver {
private:
  int N;
  std::vector<double> omega2;
  std::vector<double> lambda;
  std::vector<double> coupling_coeff;
  std::vector<std::vector<int>> adj;
  ForceFunc forces;
  CouplingFunc coupling_func;

protected:
  void derivs(double t, const std::vector<double> &state,
              std::vector<double> &dydx) override;

  VdPEnsembleSolver(double t0, const std::vector<double> &y0,
                    const std::vector<double> &freqs,
                    const std::vector<double> &lambda,
                    const std::vector<double> &coupling,
                    const std::vector<std::vector<int>> &adj, ForceFunc &func,
                    CouplingFunc coupling_func);

public:
  static tl::expected<VdPEnsembleSolver, ConstructError>
  create(size_t N, double t0, const std::vector<double> &y0,
         const std::vector<double> &freqs, const std::vector<double> &lambda,
         const std::vector<double> &coupling,
         const std::vector<std::vector<int>> &adj, ForceFunc &func,
         CouplingFunc coupling_func = CouplingFunc{}) {
    if (N != lambda.size() || N != coupling.size() || N != adj.size() ||
        N != adj[0].size() || 2 * N != y0.size()) {
      return tl::make_unexpected(ConstructError::WrongArgSize);
    }
    return VdPEnsembleSolver(t0, y0, freqs, lambda, coupling, adj, func,
                             coupling_func);
  }
};

template <class ForceFunc, class CouplingFunc>
VdPEnsembleSolver<ForceFunc, CouplingFunc>::VdPEnsembleSolver(
    double t0, const std::vector<double> &y0, const std::vector<double> &freqs,
    const std::vector<double> &lambda, const std::vector<double> &coupling,
    const std::vector<std::vector<int>> &adj, ForceFunc &func,
    CouplingFunc coupling_func)
    : RK4Solver(t0, y0), N(freqs.size()), lambda(lambda),
      coupling_coeff(coupling), adj(adj), forces(func),
      coupling_func(coupling_func) {
  omega2.resize(N);
  for (int i = 0; i < N; ++i) {
    omega2[i] = freqs[i] * freqs[i];
  }
}

template <class ForceFunc, class CouplingFunc>
void VdPEnsembleSolver<ForceFunc, CouplingFunc>::derivs(
    double t, const std::vector<double> &state, std::vector<double> &dydx) {
  const auto &impacts = forces(t);

  for (int i = 0; i < N; ++i) {
    int idx_x = 2 * i;
    int idx_y = 2 * i + 1;
    double xi = state[idx_x];
    double yi = state[idx_y];

    double coupling_sum = coupling_func(i, state, adj);

    dydx[idx_x] = yi;
    dydx[idx_y] = (lambda[i] - xi * xi) * yi - omega2[i] * xi +
                  coupling_coeff[i] * coupling_sum;
  }

  for (const auto &impact : impacts) {
    if (impact.first < (int)dydx.size()) {
      dydx[impact.first] += impact.second;
    }
  }
}

} // namespace vdp_ensemble
