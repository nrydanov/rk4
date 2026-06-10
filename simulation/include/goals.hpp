#ifndef goals_hpp
#define goals_hpp

#include <cmath>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

template<typename Real, typename Index>
struct ampl_t {

  static constexpr Real TINY = 100 * std::numeric_limits<Real>::epsilon();

  const Index dim_u;
  const Index nods;
  const Index size;
  const static std::string name;

  ampl_t(Index dim_u, Index nods) :
    dim_u(dim_u),
    nods(nods),
    size(dim_u * nods)
  {}

  Real operator()(const Real *const uu) const {
    Real A = Real(0);
    Real A2 = Real(0);
    for (Index i0 = 0; i0 < dim_u; ++i0) {
      Real sum = Real(0);
      for (Index i = i0; i < size; i += dim_u) {
        Real xx = uu[i];
        A2 += xx * xx;
        sum += xx;
      }
      A += sum * sum;
    }
    return A / (A2 * nods + TINY);
  }
};

template<typename Real, typename Index>
const std::string ampl_t<Real, Index>::name = "ampl_t";


template<typename Real, typename Index>
struct phase_t {

  static constexpr Real pi = std::numbers::pi_v<Real>;
  static constexpr Real two_pi = Real(2) * std::numbers::pi_v<Real>;

  const Index dim_u;
  const Index nods;
  const Index nods1;
  const Real inv_norm;
  std::vector<Real> phi;
  static const std::string name;

  phase_t(Index dim_u, Index nods) :
    dim_u(dim_u),
    nods(nods), nods1(nods-1),
    inv_norm( Real(2) / (pi * nods * nods1) ),
    phi(nods)
  {}

  Real operator()(const Real *const uu) {
    const Real *ux = uu;
    for (Index i = 0; i < nods; ++i, ux += dim_u) {
      phi[i] = std::atan2(ux[1], ux[0]);
    }
    return from_phases();
  }

  // То же, но на вход подаются уже вычисленные фазы raw[i] = atan2(y_i, x_i) —
  // чтобы не считать atan2 повторно.
  Real from_raw_phases(const Real *raw) {
    for (Index i = 0; i < nods; ++i) {
      phi[i] = raw[i];
    }
    return from_phases();
  }

private:

  // Агрегирует попарные разности фаз из phi[] в метрику когерентности.
  Real from_phases() const {
    Real Phi = Real(0);
    for (Index i = 0; i < nods1; ++i) {
      const Real phi_i = phi[i];
      for (Index j = i + 1; j < nods; ++j) {
        Phi += angle_abs_diff(phi_i, phi[j]);
      }
    }
    Phi *= inv_norm;
    Phi = Real(1) - Phi;
    return Phi;
  }

  static Real angle_abs_diff(Real phi1, Real phi2) {
    Real dd = std::abs(phi1 - phi2);
    if (dd > pi) {dd = two_pi - dd;}
    return dd;
  }

};

template<typename Real, typename Index>
const std::string phase_t<Real, Index>::name = "phase_t";

namespace goals {

// L = (x0 + x1 + ... + x_{N-1})^2 — мгновенное значение
inline double coherence(const double *state, int N) {
  double sum = 0.0;
  for (int i = 0; i < N; ++i) sum += state[2 * i];
  return sum * sum;
}

// ampl_t — мгновенное значение
inline double ampl(const double *state, int N) {
  ampl_t<double, int> goal(2, N);
  return goal(state);
}

// phase_t — мгновенное значение
inline double phase_coherence(const double *state, int N) {
  phase_t<double, int> goal(2, N);
  return goal(state);
}

} // namespace goals

#endif
