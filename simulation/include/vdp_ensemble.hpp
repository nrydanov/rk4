#pragma once
#include "coupling.hpp"
#include "forces.hpp"
#include "solver.hpp"
#include <array>

namespace vdp_ensemble {

template <int N = 3, class ForceFunc = forces::NoopForce,
          class CouplingFunc = coupling::Inertial>
class VdPEnsembleSolver
    : public RK4Solver<VdPEnsembleSolver<N, ForceFunc, CouplingFunc>, 2 * N> {
  friend RK4Solver<VdPEnsembleSolver, 2 * N>;

public:
  using State = std::array<double, 2 * N>;
  using Params = std::array<double, N>;
  using Adj = std::array<std::array<int, N>, N>;

  VdPEnsembleSolver(const State &y0, const Params &freqs, const Params &mu,
                    const Params &coupling, const Adj &adj, ForceFunc func,
                    CouplingFunc cf = {})
      : mu(mu), coupling_coeff(coupling), adj(adj), forces(func),
        coupling_func(cf) {
    this->y = y0;
    for (int i = 0; i < N; ++i) omega2[i] = freqs[i] * freqs[i];
  }

  void reset(const State &y0, const Params &freqs, const Params &coupling) {
    this->x = 0.0;
    this->y = y0;
    coupling_coeff = coupling;
    for (int i = 0; i < N; ++i) omega2[i] = freqs[i] * freqs[i];
  }

private:
  Params omega2{}, mu, coupling_coeff;
  Adj adj;
  ForceFunc forces;
  CouplingFunc coupling_func;

  void derivs(double t, const State &state, State &dydx) {
    const auto &impacts = forces(t);
    for (int i = 0; i < N; ++i) {
      const double xi = state[2 * i], yi = state[2 * i + 1];
      dydx[2 * i] = yi;
      dydx[2 * i + 1] = (mu[i] - xi * xi) * yi - omega2[i] * xi +
                        coupling_coeff[i] * coupling_func(i, state, adj);
    }
    for (const auto &im : impacts)
      if (im.first < 2 * N) dydx[im.first] += im.second;
  }
};

} // namespace vdp_ensemble
