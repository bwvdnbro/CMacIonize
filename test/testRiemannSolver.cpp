/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file testRiemannSolver.cpp
 *
 * @brief Unit test for the RiemannSolver class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "ExactRiemannSolver.hpp"
#include <fstream>
#include <string>

/**
 * @brief Run a Rieman solver test.
 *
 * @param solver RiemannSolver to test.
 * @param rhoL Density of the left state.
 * @param uL Velocity of the left state.
 * @param PL Pressure of the left state.
 * @param rhoR Density of the right state.
 * @param uR Velocity of the right state.
 * @param PR Pressure of the right state.
 * @param rhoexp Expected density solution.
 * @param uexp Expected velocity solution.
 * @param Pexp Expected pressure solution.
 */
void run_test(ExactRiemannSolver &solver, double rhoL, double uL, double PL,
              double rhoR, double uR, double PR, double rhoexp, double uexp,
              double Pexp) {

  double rhosol, usol, Psol;
  solver.solve(rhoL, uL, PL, rhoR, uR, PR, rhosol, usol, Psol);
  assert_values_equal_rel(rhosol, rhoexp, 1.e-4);
  assert_values_equal_rel(usol, uexp, 1.e-4);
  assert_values_equal_rel(Psol, Pexp, 1.e-4);

  // symmetry test: also test the problem with left and right states reversed
  solver.solve(rhoR, -uR, PR, rhoL, -uL, PL, rhosol, usol, Psol);
  assert_values_equal_rel(rhosol, rhoexp, 1.e-4);
  assert_values_equal_rel(usol, -uexp, 1.e-4);
  assert_values_equal_rel(Psol, Pexp, 1.e-4);
}

/**
 * @brief Plot the solution for the Riemann problem with the given left and
 * right states at the given time.
 *
 * @param solver RiemannSolver to use.
 * @param rhoL Density of the left state.
 * @param uL Velocity of the left state.
 * @param PL Pressure of the left state.
 * @param rhoR Density of the right state.
 * @param uR Velocity of the right state.
 * @param PR Pressure of the right state.
 * @param t Time at which to evaluate the solution.
 * @param filename Name of the file to write out.
 */
void plot_solution(ExactRiemannSolver &solver, double rhoL, double uL,
                   double PL, double rhoR, double uR, double PR, double t,
                   std::string filename) {

  std::ofstream ofile(filename);
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    double x = (i + 0.5) * 0.001 - 0.5;
    double dxdt = x / t;
    double rhosol, usol, Psol;
    int flag =
        solver.solve(rhoL, uL, PL, rhoR, uR, PR, rhosol, usol, Psol, dxdt);
    if (i == 0) {
      assert_condition(flag == -1);
    }
    if (i == 999) {
      assert_condition(flag == 1);
    }
    ofile << x << "\t" << rhosol << "\t" << usol << "\t" << Psol << "\n";
  }
}

/**
 * @brief Unit test for the RiemannSolver class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  ExactRiemannSolver solver(5. / 3.);

  // Toro tests
  run_test(solver, 1., 0., 1., 0.125, 0., 0.1, 0.47969, 0.841194, 0.293945);
  run_test(solver, 1., -2., 0.4, 1., 2., 0.4, 0.00617903, 0., 8.32249e-05);
  run_test(solver, 1., 0., 1000., 1., 0., 0.01, 0.615719, 18.2812, 445.626);
  run_test(solver, 1., 0., 0.01, 1., 0., 100., 0.61577, -5.78011, 44.5687);
  run_test(solver, 5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.0950,
           12.743, 8.56045, 1841.82);
  // vacuum generation
  run_test(solver, 1., -1., 1.e-6, 1., 1., 1.0005e-6, 0., 0., 0.);

  plot_solution(solver, 1., 0., 1., 0.125, 0., 0.1, 0.25,
                "test_riemann_test1.txt");
  plot_solution(solver, 1., -2., 0.4, 1., 2., 0.4, 0.15,
                "test_riemann_test2.txt");
  plot_solution(solver, 1., 0., 1000., 1., 0., 0.01, 0.012,
                "test_riemann_test3.txt");
  plot_solution(solver, 1., 0., 0.01, 1., 0., 100., 0.035,
                "test_riemann_test4.txt");
  plot_solution(solver, 5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.0950,
                0.035, "test_riemann_test5.txt");

  // very small values test
  ExactRiemannSolver solver2(1.);
  const double cs = 2875.;
  const double PL = 2.21538e-153;
  const double PR = 1.54009e-154;
  const double udiff = -24.3203;
  const double gamma = 1.00000001;
  const double uL = 0.;

  const double cs2 = cs * cs;
  const double rhoL = (gamma * PL) / cs2;
  const double rhoR = (gamma * PR) / cs2;
  const double uR = udiff + uL;

  assert_condition(rhoL == rhoL);
  assert_condition(rhoR == rhoR);

  plot_solution(solver2, rhoL, uL, PL, rhoR, uR, PR, 0.01,
                "test_riemann_small.txt");

  return 0;
}
