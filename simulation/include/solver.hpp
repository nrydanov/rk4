#pragma once
#include <vector>

template<typename Derived>
class RK4Solver {
public:
  double getTime() const { return x; }
  const std::vector<double> &getState() const { return y; }

  void step(double h) {
    const int n = y.size();
    const double hh = 0.5 * h;
    const double h6 = h / 6.0;
    const double xh = x + hh;

    derived().derivs(x, y, dydx);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * dydx[i];
    derived().derivs(xh, yt, dyt);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * dyt[i];
    derived().derivs(xh, yt, dym);
    for (int i = 0; i < n; ++i) {
      yt[i] = y[i] + h * dym[i];
      dym[i] += dyt[i];
    }
    derived().derivs(x + h, yt, dyt);
    for (int i = 0; i < n; ++i)
      y[i] += h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    x += h;
  }

protected:
  double x;
  std::vector<double> y;
  std::vector<double> dydx, dyt, dym, yt;

  RK4Solver(double t0, const std::vector<double> &y0)
      : x(t0), y(y0), dydx(y0.size()), dyt(y0.size()), dym(y0.size()),
        yt(y0.size()) {}

private:
  Derived &derived() { return static_cast<Derived &>(*this); }
};

class VanDerPolSolver : public RK4Solver<VanDerPolSolver> {
public:
  VanDerPolSolver(double t0, const std::vector<double> &y0, double lambda,
                  double omega)
      : RK4Solver(t0, y0), lambda(lambda), omega2(omega * omega) {}

  void derivs(double, const std::vector<double> &state,
              std::vector<double> &dydx) {
    dydx[0] = state[1];
    dydx[1] = (lambda - state[0] * state[0]) * state[1] - omega2 * state[0];
  }

private:
  double lambda;
  double omega2;
};
