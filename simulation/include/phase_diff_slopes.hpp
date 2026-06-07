/*
  25.03.26 Вычисление линейной аппроксимации разности фаз. Если
  уголовой коэффциент (слоп) равен нулю - имеем фазовую синхронизацию

  Пусть имеем три связанных двумерных систем (три ван дер Поля). Пусть в массиве uu
  записано текущее решение: uu = {x1, y1, x2, y2, x3, y3}. Тогда
  
  nods = 3
  local_dim = 2
  ix = 0
  iy = 1

  Тогда используем этот класс вот так:
  
  phase_diff_slopes<double> pds(nods, local_dim, ix, iy);
  for (int nn = 0; nn < n_store; ++nn) {
    sys.step();
    pds.push(sys.tt, sys.uu);
  }
  for (int ii = 0; ii < nods; ++ii) {
    for (int jj = ii + 1; jj < nods; ++jj) {
      cout << pds(ii, jj) << endl; // Это будут слопы для разности фаз подсистем ii и jj
    }
  }  
*/

#ifndef phase_diff_slopes_hpp
#define phase_diff_slopes_hpp

// #include <cmath>
#include <vector>

#include "linfit_stream.hpp"
#include "unwrap_phase.hpp"

template<typename Real>
struct phase_diff_slopes {
public:
  
  phase_diff_slopes(int nods, int local_dim, int ix, int iy) :
    nods(nods), local_dim(local_dim), ix(ix), iy(iy),
    n_phases((nods*nods-nods)/2), 
    phase(nods), slope(n_phases), dirty(true)
  {}

  void reset() {
    for (auto &p : phase) p = unwrap_phase<Real>{};
    for (auto &s : slope) s.reset();
    dirty = true;
  }

  void push(Real tt, const Real *uu) {
    for (int ii = 0, mm = 0; ii < nods; ++ii, mm += local_dim) {
      phase[ii](uu[mm+ix], uu[mm+iy]);
    }
    accumulate(tt);
  }

  // То же, но на вход подаются уже вычисленные сырые фазы raw[ii] =
  // atan2(y_ii, x_ii) — чтобы не считать atan2 повторно.
  void push_raw(Real tt, const Real *raw) {
    for (int ii = 0; ii < nods; ++ii) {
      phase[ii].push_raw(raw[ii]);
    }
    accumulate(tt);
  }

  Real operator()(int ii, int jj) {
    if (dirty) {
      for (auto &ss: slope) {ss.load();}
      dirty = false;
    }
    if (ii == jj) return 0.0;
    if (ii > jj) {int tmp = ii; ii = jj; jj = tmp;}
    // Пересчёт пары ii, jj в линейный индекс
    int mm = ii*(2*nods-ii-1)/2+jj-ii-1;
    return std::abs(slope[mm].A);
  }
  
private:

  // Общая часть push/push_raw: фазы узлов уже разложены в phase[], осталось
  // протолкнуть попарные разности в линейные регрессии.
  void accumulate(Real tt) {
    dirty = true;
    int mm = 0;
    for (int ii = 0; ii < nods; ++ii) {
      Real phi1 = phase[ii].phi;
      for (int jj = ii + 1; jj < nods; ++jj) {
        Real phi2 = phase[jj].phi;
        slope[mm++].push(tt, phi1 - phi2);
      }
    }
  }

  static constexpr bool LF_WITH_B = 1;
  static constexpr bool LF_NEED_ERR = 0;
  using LinFit = linfit_stream<LF_WITH_B, LF_NEED_ERR, Real>;

  int nods;
  int local_dim;
  int ix;
  int iy;
  int n_phases;
  std::vector<unwrap_phase<Real>> phase;
  std::vector<LinFit> slope;
  bool dirty;

};

#endif
