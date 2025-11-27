#include <iostream>

class ODESolver {
public:
  virtual ~ODESolver() = default;

protected:
  double x;
  std::vector<double> y;
  ODESolver(double t0, const std::vector<double> &y0) : x(t0), y(y0) {}

  virtual void derivs(double x, const std::vector<double> &y,
                      std::vector<double> &dydx) = 0;
  virtual void step(double h) = 0;
  virtual double getTime() const { return x; }
  virtual const std::vector<double> &getState() const { return y; }
};

class RK4Solver : public ODESolver {
  void step(double h) override {
    const int n = y.size();
    const double hh = 0.5 * h;
    const double h6 = h / 6.0;
    const double xh = x + hh;

    // k1 = dydx * h
    derivs(x, y, dydx);

    for (int i = 0; i < n; ++i) {
      yt[i] = y[i] + hh * dydx[i];
    }
    // k2 = dyt * h
    derivs(xh, yt, dyt);

    for (int i = 0; i < n; ++i) {
      yt[i] = y[i] + hh * dyt[i];
    }
    // k3 = dym * h
    derivs(xh, yt, dym);
    for (int i = 0; i < n; ++i) {
      yt[i] = y[i] + h * dym[i];
      // dym = k2 + k3
      dym[i] += dyt[i];
    }
    // k4 = dyt * h
    derivs(x + h, yt, dyt);

    for (int i = 0; i < n; ++i) {
      y[i] += h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    }

    x += h;
  }

protected:
  RK4Solver(double t0, const std::vector<double> &y0)
      : ODESolver(t0, y0), dydx(y0.size()), dyt(y0.size()), dym(y0.size()),
        yt(y0.size()) {}

private:
  std::vector<double> dydx, dyt, dym, yt;
};

class VanDerPolSolver : public RK4Solver {
  void derivs(double t, const std::vector<double> &y,
              std::vector<double> &dydx) override {
    dydx[0] = y[1];
    dydx[1] = (lambda - y[0] * y[0]) * y[1] - omega2 * y[0];
  }

private:
  double lambda;
  double omega2;

public:
  VanDerPolSolver(double t0, const std::vector<double> &y0, double lambda,
                  double omega)
      : RK4Solver(t0, y0), lambda(lambda), omega2(omega * omega) {};
};
