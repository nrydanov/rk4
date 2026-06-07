#pragma once
#include <cmath>
#include <tl/expected.hpp>
#include <vector>

namespace forces {

enum class ConstructError { WrongArgSize = 0 };

class SinForce {
public:
  using ForceResult = std::vector<std::pair<int, double>>;

  static tl::expected<SinForce, ConstructError>
  create(const std::vector<int> &pins, const std::vector<double> &alpha,
         const std::vector<double> &fomg) {
    if (pins.size() != alpha.size() || pins.size() != fomg.size())
      return tl::make_unexpected(ConstructError::WrongArgSize);
    return SinForce(pins, alpha, fomg);
  }

  const ForceResult &operator()(double t) {
    for (size_t i = 0; i < force.size(); ++i)
      force[i].second = amps[i] * std::sin(freqs[i] * t);
    return force;
  }

private:
  std::vector<double> amps;
  std::vector<double> freqs;
  ForceResult force;

  SinForce(const std::vector<int> &pins, const std::vector<double> &alpha,
           const std::vector<double> &fomg)
      : amps(alpha), freqs(fomg) {
    force.resize(pins.size());
    for (size_t i = 0; i < pins.size(); ++i) {
      force[i].first = 2 * pins[i] + 1;
      force[i].second = 0.0;
    }
  }
};

class NoopForce {
public:
  using ForceResult = std::vector<std::pair<int, double>>;
  const ForceResult &operator()(double) { return result; }

private:
  ForceResult result;
};

} // namespace forces
