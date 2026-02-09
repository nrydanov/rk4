#include "forces.h"

SinForce::SinForce(const std::vector<int> &pins,
                   const std::vector<double> &alpha,
                   const std::vector<double> &fomg)
    : amps(alpha), freqs(fomg) {
  force.resize(pins.size());
  for (size_t i = 0; i < pins.size(); ++i) {
    force[i].first = 2 * pins[i] + 1;
    force[i].second = 0.0;
  }
}

const SinForce::ForceResult &SinForce::operator()(double t) {
  for (size_t i = 0; i < force.size(); ++i) {
    force[i].second = amps[i] * std::sin(freqs[i] * t);
  }
  return force;
}
