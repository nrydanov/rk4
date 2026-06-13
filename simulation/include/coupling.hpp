#pragma once
#include <cstddef>
#include <numeric>

namespace coupling {

struct Inertial {
  template <class State, class Adj>
  double operator()(int i, const State &s, const Adj &adj) const {
    double sum = 0.0, xi = s[2 * i];
    for (std::size_t j = 0; j < adj[i].size(); ++j)
      sum += adj[i][j] * (s[2 * j] - xi);
    return sum;
  }
};

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

struct Dissipative {
  template <class State, class Adj>
  double operator()(int i, const State &s, const Adj &adj) const {
    double sum = 0.0, yi = s[2 * i + 1];
    for (std::size_t j = 0; j < adj[i].size(); ++j)
      sum += adj[i][j] * (s[2 * j + 1] - yi);
    return sum;
  }
};

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
