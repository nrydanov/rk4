#pragma once
#include <vector>

class ODESolver {
public:
  virtual ~ODESolver() = default;
  virtual double getTime() const { return x; }
  virtual const std::vector<double> &getState() const { return y; }
  virtual void step(double h) = 0;

protected:
  double x;
  std::vector<double> y;
  ODESolver(double t0, const std::vector<double> &y0);
  virtual void derivs(double x, const std::vector<double> &state,
                      std::vector<double> &dydx) = 0;
};

class RK4Solver : public ODESolver {
public:
  void step(double h) override;

protected:
  RK4Solver(double t0, const std::vector<double> &y0);

private:
  std::vector<double> dydx, dyt, dym, yt;
};

class VanDerPolSolver : public RK4Solver {
public:
  VanDerPolSolver(double t0, const std::vector<double> &y0, double lambda,
                  double omega);

protected:
  void derivs(double t, const std::vector<double> &state,
              std::vector<double> &dydx) override;

private:
  double lambda;
  double omega2;
};
