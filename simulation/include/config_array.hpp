#pragma once
#include <array>
#include <cstddef>
#include <vector>

template <std::size_t M>
std::array<double, M> to_array(const std::vector<double> &v) {
  std::array<double, M> a{};
  for (std::size_t i = 0; i < M && i < v.size(); ++i) a[i] = v[i];
  return a;
}

template <int N>
std::array<std::array<int, N>, N>
to_adj(const std::vector<std::vector<int>> &m) {
  std::array<std::array<int, N>, N> a{};
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) a[i][j] = m[i][j];
  return a;
}
