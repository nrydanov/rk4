#pragma once
#include <array>
#include <cstddef>

// CRTP RK4-интегратор. Dim известен на этапе компиляции => состояние и рабочие
// буферы лежат в std::array (стек, ноль heap), циклы разворачиваются, derivs()
// инлайнится из производного класса вместо виртуального вызова.
template <typename Derived, std::size_t Dim>
class RK4Solver {
public:
  double getTime() const { return x; }
  const std::array<double, Dim> &getState() const { return y; }

  void step(double h) {
    const double hh = 0.5 * h, h6 = h / 6.0, xh = x + hh;
    derived().derivs(x, y, dydx);
    for (std::size_t i = 0; i < Dim; ++i) yt[i] = y[i] + hh * dydx[i];
    derived().derivs(xh, yt, dyt);
    for (std::size_t i = 0; i < Dim; ++i) yt[i] = y[i] + hh * dyt[i];
    derived().derivs(xh, yt, dym);
    for (std::size_t i = 0; i < Dim; ++i) {
      yt[i] = y[i] + h * dym[i];
      dym[i] += dyt[i];
    }
    derived().derivs(x + h, yt, dyt);
    for (std::size_t i = 0; i < Dim; ++i)
      y[i] += h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    x += h;
  }

protected:
  double x = 0.0;
  std::array<double, Dim> y{}, dydx{}, dyt{}, dym{}, yt{};

private:
  Derived &derived() { return static_cast<Derived &>(*this); }
};
