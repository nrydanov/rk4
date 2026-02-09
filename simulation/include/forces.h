#pragma once
#include <vector>

class SinForce {
public:
  using ForceResult = std::vector<std::pair<int, double>>;

  SinForce(const std::vector<int> &pins, const std::vector<double> &alpha,
           const std::vector<double> &fomg);

  const ForceResult &operator()(double t);

private:
  std::vector<double> amps;
  std::vector<double> freqs;
  ForceResult force;
};
