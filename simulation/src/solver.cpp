#include "solver.h"

ODESolver::ODESolver(double t0, const std::vector<double> &y0) : x(t0), y(y0) {}

RK4Solver::RK4Solver(double t0, const std::vector<double> &y0)
    : ODESolver(t0, y0), dydx(y0.size()), dyt(y0.size()), dym(y0.size()),
      yt(y0.size()) {}

void RK4Solver::step(double h) {
  const int n = y.size();
  const double hh = 0.5 * h;
  const double h6 = h / 6.0;
  const double xh = x + hh;

  derivs(x, y, dydx);
  for (int i = 0; i < n; ++i)
    yt[i] = y[i] + hh * dydx[i];
  derivs(xh, yt, dyt);
  for (int i = 0; i < n; ++i)
    yt[i] = y[i] + hh * dyt[i];
  derivs(xh, yt, dym);
  for (int i = 0; i < n; ++i) {
    yt[i] = y[i] + h * dym[i];
    dym[i] += dyt[i];
  }
  derivs(x + h, yt, dyt);
  for (int i = 0; i < n; ++i) {
    y[i] += h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
  }
  x += h;
}

VanDerPolSolver::VanDerPolSolver(double t0, const std::vector<double> &y0,
                                 double lambda, double omega)
    : RK4Solver(t0, y0), lambda(lambda), omega2(omega * omega) {}

void VanDerPolSolver::derivs(double t, const std::vector<double> &state,
                             std::vector<double> &dydx) {
  (void)t;
  dydx[0] = state[1];
  dydx[1] = (lambda - state[0] * state[0]) * state[1] - omega2 * state[0];
}
