#pragma once
#include <cstddef>
#include <numeric>

// Функторы связи. operator() шаблонизирован по типам контейнеров состояния и
// матрицы смежности — работает и с std::array (compile-time N), и с std::vector.

namespace coupling {

// Инерционная связь через x, без нормировки (лапласиан)
// s_n = sum_j a_nj * (x_j - x_n)
struct Inertial {
  template <class State, class Adj>
  double operator()(int i, const State &s, const Adj &adj) const {
    double sum = 0.0, xi = s[2 * i];
    for (std::size_t j = 0; j < adj[i].size(); ++j)
      sum += adj[i][j] * (s[2 * j] - xi);
    return sum;
  }
};

// Инерционная связь через x, с нормировкой на число соседей
// s_n = sum_j (a_nj / k_n) * (x_j - x_n)
struct InertialNorm {
  template <class State, class Adj>
  double operator()(int i, const State &s, const Adj &adj) const {
    double sum = 0.0, xi = s[2 * i];
    int k = std::accumulate(adj[i].begin(), adj[i].end(), 0);
    if (k == 0) return 0.0;
    for (std::size_t j = 0; j < adj[i].size(); ++j)
      sum += adj[i][j] * (s[2 * j] - xi);
    return sum / k;
  }
};

// Диссипативная связь через y, без нормировки
// c_n = sum_j a_nj * (y_j - y_n)
struct Dissipative {
  template <class State, class Adj>
  double operator()(int i, const State &s, const Adj &adj) const {
    double sum = 0.0, yi = s[2 * i + 1];
    for (std::size_t j = 0; j < adj[i].size(); ++j)
      sum += adj[i][j] * (s[2 * j + 1] - yi);
    return sum;
  }
};

// Диссипативная связь через y, с нормировкой на число соседей
// c_n = sum_j (a_nj / k_n) * (y_j - y_n)
struct DissipativeNorm {
  template <class State, class Adj>
  double operator()(int i, const State &s, const Adj &adj) const {
    double sum = 0.0, yi = s[2 * i + 1];
    int k = std::accumulate(adj[i].begin(), adj[i].end(), 0);
    if (k == 0) return 0.0;
    for (std::size_t j = 0; j < adj[i].size(); ++j)
      sum += adj[i][j] * (s[2 * j + 1] - yi);
    return sum / k;
  }
};

} // namespace coupling
