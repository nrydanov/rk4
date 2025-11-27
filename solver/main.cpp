#include "CLI11.hpp"
#include "solver.cpp"
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
  CLI::App app{"Van der Pol Equation Solver"};
  // Начальные условия
  double lambda, omega, x0, v0, h, t_max;
  std::string output;
  app.add_option("-l,--lambda", lambda, "Lambda parameter")->default_val(1.0);
  app.add_option("-w,--omega", omega, "Omega frequency")->default_val(1.0);
  app.add_option("-x,--x0", x0, "Initial position")->default_val(2.0);
  app.add_option("-v,--v0", v0, "Initial velocity")->default_val(0.0);
  app.add_option("-s,--step", h, "Integration step size")->default_val(0.01);
  app.add_option("-t,--tmax", t_max, "Simulation time")->default_val(50.0);
  app.add_option("-o,--output", output, "Output file")
      ->default_val("output.csv");

  CLI11_PARSE(app, argc, argv);
  std::vector<double> y0 = {x0, v0};

  VanDerPolSolver solver(0.0, y0, lambda, omega);

  std::ofstream file(output);
  file << "t,x,v\n";

  file << solver.getTime() << "," << solver.getState()[0] << ","
       << solver.getState()[1] << "\n";

  while (solver.getTime() < t_max) {
    file << solver.getTime() << "," << solver.getState()[0] << ","
         << solver.getState()[1] << "\n";
    solver.step(h);
  }

  std::cout << "Results are saved to" << " " << output;
  return 0;
}
