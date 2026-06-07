/*
  25.03.26 Восстановление непрерывной фазы. На вход подаются текущие xx и yy, на
  выходе - соответствующая фаза
*/

#ifndef unwrap_phase_hpp
#define unwrap_phase_hpp

#include <cmath>
#include <numbers>

template <typename Real> struct unwrap_phase {

  Real phi;

  unwrap_phase() : phi(0.0), first(true) {}

  Real operator()(Real xx, Real yy) { return push_raw(std::atan2(yy, xx)); }

  // Разворачивает уже вычисленную сырую фазу raw_phi = atan2(yy, xx). Позволяет
  // переиспользовать atan2, посчитанный один раз на стороне вызова.
  Real push_raw(Real raw_phi) {
    if (first) {
      phi = prev_raw_phi = raw_phi;
      first = false;
    } else {
      Real delta = raw_phi - prev_raw_phi;
      if (delta > PI)
        delta -= 2.0 * PI;
      else if (delta < -PI)
        delta += 2.0 * PI;
      phi += delta;
      prev_raw_phi = raw_phi;
    }
    return phi;
  }

private:
  static constexpr Real PI = std::numbers::pi_v<Real>;
  bool first;
  Real prev_raw_phi;
};

#endif
