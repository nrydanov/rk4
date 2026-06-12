// Валидация вычислительного пайплайна против известных аналитических
// результатов. Каждая проверка печатает PASS/FAIL; ненулевой код возврата,
// если хоть одна провалена.
//
// 1. Метрики на синтетических состояниях: когерентное (A=1, P=1), противофаза
//    (P=0), восстановление известного слопа разности фаз.
// 2. Одиночный ВдП в форме dy = (mu - x^2) y - w^2 x: амплитуда предельного
//    цикла 2*sqrt(mu) (энергобаланс: mu A^2/2 = A^4/8), частота 1 - mu^2/16
//    (поправка Линдштедта—Пуанкаре; mu здесь совпадает с параметром
//    нормированной формы после замены x -> sqrt(mu) u).
// 3. Порядок сходимости RK4: ошибка падает в ~16 раз при dt/2.
// 4. Два идентичных осциллятора с диссипативной связью: полная синхронизация.
// 5. Язык Адлера: граница захвата Delta_c растёт линейно по eps.
// 6. Симметрия звезды N=3 относительно перестановки листьев 1<->2.

#include "coupling.hpp"
#include "forces.hpp"
#include "goals.hpp"
#include "linfit_stream.hpp"
#include "phase_diff_slopes.hpp"
#include "unwrap_phase.hpp"
#include "vdp_ensemble.hpp"
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>

namespace {

int n_pass = 0, n_fail = 0;

void report(const std::string &name, bool pass, const std::string &detail) {
  std::printf("[%s] %-42s %s\n", pass ? "PASS" : "FAIL", name.c_str(),
              detail.c_str());
  (pass ? n_pass : n_fail)++;
}

std::string fmt(const char *f, ...) {
  char buf[256];
  va_list ap;
  va_start(ap, f);
  std::vsnprintf(buf, sizeof buf, f, ap);
  va_end(ap);
  return buf;
}

// --- 1. Метрики на синтетике -----------------------------------------------

void check_metrics_synthetic() {
  // Полностью когерентное состояние: все осцилляторы в одной точке
  const double coh[6] = {1.3, 0.7, 1.3, 0.7, 1.3, 0.7};
  ampl_t<double, int> ampl(2, 3);
  phase_t<double, int> phase(2, 3);
  double A = ampl(coh), P = phase(coh);
  report("metrics: coherent state A=1, P=1",
         std::abs(A - 1.0) < 1e-12 && std::abs(P - 1.0) < 1e-12,
         fmt("A=%.15f P=%.15f", A, P));

  // Противофаза двух осцилляторов: |dphi| = pi => P = 0
  const double anti[4] = {1.0, 0.0, -1.0, 0.0};
  phase_t<double, int> phase2(2, 2);
  double P2 = phase2(anti);
  report("metrics: anti-phase pair P=0", std::abs(P2) < 1e-12,
         fmt("P=%.15f", P2));

  // Синтетические фазы phi1 = 0.3 t, phi2 = 0.1 t => слоп разности 0.2
  phase_diff_slopes<double> pds(2, 2, 0, 1);
  for (double t = 0.0; t < 100.0; t += 0.05) {
    double raw[2] = {std::remainder(0.3 * t, 2 * M_PI),
                     std::remainder(0.1 * t, 2 * M_PI)};
    pds.push_raw(t, raw);
  }
  double s = pds(0, 1);
  report("metrics: known phase-diff slope 0.2", std::abs(s - 0.2) < 1e-9,
         fmt("slope=%.12f", s));
}

// --- 2. Одиночный ВдП: амплитуда и частота ----------------------------------

void check_single_vdp() {
  constexpr double mu = 0.1, dt = 0.01;
  using Solver = vdp_ensemble::VdPEnsembleSolver<1>;
  Solver solver({2.5, 0.0}, {1.0}, {mu}, {0.0}, {{{0}}}, forces::NoopForce{});

  for (double t = 0.0; t < 300.0; t += dt) solver.step(dt);

  double amp = 0.0;
  unwrap_phase<double> uw;
  linfit_stream<true, false, double> fit;
  for (double t = 300.0; t < 800.0; t += dt) {
    solver.step(dt);
    const auto &s = solver.getState();
    amp = std::max(amp, std::abs(s[0]));
    fit.push(t, uw(s[0], s[1]));
  }
  fit.load();
  double freq = std::abs(fit.A);
  const double freq_theory = 1.0 - mu * mu / 16.0;
  const double amp_theory = 2.0 * std::sqrt(mu);

  report("single VdP: limit cycle amplitude 2*sqrt(mu)",
         std::abs(amp - amp_theory) < 0.003,
         fmt("amp=%.6f theory=%.6f", amp, amp_theory));
  report("single VdP: frequency ~ 1 - mu^2/16",
         std::abs(freq - freq_theory) < 2e-3,
         fmt("freq=%.6f theory=%.6f", freq, freq_theory));
}

// --- 3. Порядок сходимости RK4 ----------------------------------------------

std::array<double, 2> integrate_single(double dt, double T) {
  using Solver = vdp_ensemble::VdPEnsembleSolver<1>;
  Solver solver({2.0, 0.0}, {1.0}, {0.1}, {0.0}, {{{0}}}, forces::NoopForce{});
  long n = std::lround(T / dt);
  for (long i = 0; i < n; ++i) solver.step(dt);
  return solver.getState();
}

void check_rk4_order() {
  const double T = 20.0;
  auto ref = integrate_single(0.000625, T);
  auto u1 = integrate_single(0.02, T);
  auto u2 = integrate_single(0.01, T);
  double e1 = std::hypot(u1[0] - ref[0], u1[1] - ref[1]);
  double e2 = std::hypot(u2[0] - ref[0], u2[1] - ref[1]);
  double ratio = e1 / e2;
  report("RK4: 4th-order convergence (ratio ~ 16)",
         ratio > 12.0 && ratio < 20.0,
         fmt("err(dt)=%.3e err(dt/2)=%.3e ratio=%.2f", e1, e2, ratio));
}

// --- 4. Полная синхронизация двух идентичных осцилляторов -------------------

void check_identical_sync() {
  constexpr double dt = 0.01, eps = 0.1;
  using Solver =
      vdp_ensemble::VdPEnsembleSolver<2, forces::NoopForce,
                                      coupling::Dissipative>;
  Solver solver({2.0, 0.0, 0.5, 1.0}, {1.0, 1.0}, {0.1, 0.1}, {eps, eps},
                {{{0, 1}, {1, 0}}}, forces::NoopForce{});
  for (double t = 0.0; t < 400.0; t += dt) solver.step(dt);
  const auto &s = solver.getState();
  double diff = std::abs(s[0] - s[2]) + std::abs(s[1] - s[3]);
  report("identical pair, dissipative: full sync", diff < 1e-6,
         fmt("|x1-x2|+|y1-y2|=%.3e after T=400", diff));
}

// --- 5. Язык Адлера: Delta_c(eps) растёт линейно ----------------------------

bool locked(double eps, double delta) {
  constexpr double dt = 0.05, t_trans = 500.0, T = 3000.0;
  using Solver =
      vdp_ensemble::VdPEnsembleSolver<2, forces::NoopForce,
                                      coupling::Dissipative>;
  Solver solver({2.0, 0.0, 0.1, 1.8}, {1.0, 1.0 + delta}, {0.1, 0.1},
                {eps, eps}, {{{0, 1}, {1, 0}}}, forces::NoopForce{});
  for (double t = 0.0; t < t_trans; t += dt) solver.step(dt);
  phase_diff_slopes<double> pds(2, 2, 0, 1);
  for (double t = t_trans; t < T; t += dt) {
    solver.step(dt);
    pds.push(t, solver.getState().data());
  }
  return pds(0, 1) < 0.01;
}

double find_tongue_boundary(double eps) {
  double lo = 0.0, hi = 0.4; // locked при delta=0, beat при delta=0.4
  for (int i = 0; i < 14; ++i) {
    double mid = 0.5 * (lo + hi);
    (locked(eps, mid) ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

void check_adler_scaling() {
  double dc1 = find_tongue_boundary(0.05);
  double dc2 = find_tongue_boundary(0.10);
  double ratio = dc2 / dc1;
  bool linear = ratio > 1.7 && ratio < 2.3;
  bool prefactor_sane = dc1 / 0.05 > 0.3 && dc1 / 0.05 < 3.0;
  report("Adler tongue: Delta_c scales linearly in eps",
         linear && prefactor_sane,
         fmt("Dc(0.05)=%.4f Dc(0.10)=%.4f ratio=%.2f Dc/eps=%.2f", dc1, dc2,
             ratio, dc1 / 0.05));
}

// --- 6. Симметрия 1<->2 для звезды N=3 ---------------------------------------

struct AvgMetrics {
  double L, A, P, s01, s02, s12;
};

AvgMetrics run_star(double delta1, double delta2,
                    const std::array<double, 6> &y0) {
  constexpr int N = 3;
  constexpr double dt = 0.05, t_trans = 120.0, T = 800.0, eps = 0.05;
  using Solver =
      vdp_ensemble::VdPEnsembleSolver<N, forces::NoopForce,
                                      coupling::Dissipative>;
  Solver solver(y0, {1.0, 1.0 + delta1, 1.0 + delta2}, {0.1, 0.1, 0.1},
                {eps, eps, eps}, {{{0, 1, 1}, {1, 0, 0}, {1, 0, 0}}},
                forces::NoopForce{});
  for (double t = 0.0; t < t_trans; t += dt) solver.step(dt);

  phase_diff_slopes<double> pds(N, 2, 0, 1);
  ampl_t<double, int> ampl(2, N);
  phase_t<double, int> phase(2, N);
  double L = 0, A = 0, P = 0;
  int steps = 0;
  for (double t = t_trans; t < T; t += dt, ++steps) {
    solver.step(dt);
    const auto &s = solver.getState();
    pds.push(t, s.data());
    L += goals::coherence(s.data(), N);
    A += ampl(s.data());
    P += phase(s.data());
  }
  return {L / steps, A / steps, P / steps, pds(0, 1), pds(0, 2), pds(1, 2)};
}

void check_star_symmetry() {
  const std::array<double, 6> ic = {1.7, -0.4, 0.3, 1.1, -2.0, 0.8};
  const std::array<double, 6> ic_swapped = {1.7, -0.4, -2.0, 0.8, 0.3, 1.1};
  auto m = run_star(0.12, -0.2, ic);
  auto w = run_star(-0.2, 0.12, ic_swapped);
  double d = std::abs(m.L - w.L) + std::abs(m.A - w.A) + std::abs(m.P - w.P) +
             std::abs(m.s01 - w.s02) + std::abs(m.s02 - w.s01) +
             std::abs(m.s12 - w.s12);
  report("star N=3: relabel symmetry 1<->2", d < 1e-6,
         fmt("total metric deviation %.3e", d));
}

} // namespace

int main() {
  check_metrics_synthetic();
  check_single_vdp();
  check_rk4_order();
  check_identical_sync();
  check_adler_scaling();
  check_star_symmetry();
  std::printf("\n%d passed, %d failed\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
