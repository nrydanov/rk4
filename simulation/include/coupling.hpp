#pragma once
#include <numeric>
#include <vector>

namespace coupling {

// Инерционная связь через x, без нормировки (лапласиан)
// s_n = sum_j a_nj * (x_j - x_n)
struct Inertial {
  double operator()(int i, const std::vector<double> &state,
                    const std::vector<std::vector<int>> &adj) const {
    double sum = 0.0;
    double xi = state[2 * i];
    for (int j = 0; j < (int)adj[i].size(); ++j)
      sum += adj[i][j] * (state[2 * j] - xi);
    return sum;
  }
};

// Инерционная связь через x, с нормировкой на число соседей
// s_n = sum_j (a_nj / k_n) * (x_j - x_n)
struct InertialNorm {
  double operator()(int i, const std::vector<double> &state,
                    const std::vector<std::vector<int>> &adj) const {
    double sum = 0.0;
    double xi = state[2 * i];
    int k = std::accumulate(adj[i].begin(), adj[i].end(), 0);
    if (k == 0) return 0.0;
    for (int j = 0; j < (int)adj[i].size(); ++j)
      sum += adj[i][j] * (state[2 * j] - xi);
    return sum / k;
  }
};

// Диссипативная связь через y, без нормировки
// c_n = sum_j a_nj * (y_j - y_n)
struct Dissipative {
  double operator()(int i, const std::vector<double> &state,
                    const std::vector<std::vector<int>> &adj) const {
    double sum = 0.0;
    double yi = state[2 * i + 1];
    for (int j = 0; j < (int)adj[i].size(); ++j)
      sum += adj[i][j] * (state[2 * j + 1] - yi);
    return sum;
  }
};

// Диссипативная связь через y, с нормировкой на число соседей
// c_n = sum_j (a_nj / k_n) * (y_j - y_n)
struct DissipativeNorm {
  double operator()(int i, const std::vector<double> &state,
                    const std::vector<std::vector<int>> &adj) const {
    double sum = 0.0;
    double yi = state[2 * i + 1];
    int k = std::accumulate(adj[i].begin(), adj[i].end(), 0);
    if (k == 0) return 0.0;
    for (int j = 0; j < (int)adj[i].size(); ++j)
      sum += adj[i][j] * (state[2 * j + 1] - yi);
    return sum / k;
  }
};

} // namespace coupling
