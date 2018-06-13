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
 * @file timeRiemannSolver.cpp
 *
 * @brief Timing test for the Riemann solver.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "ExactRiemannSolver.hpp"
#include "HLLCRiemannSolver.hpp"
#include "TimingTools.hpp"
#include <vector>

/**
 * @brief Timing test for the Riemann solver.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  timingtools_init("timeRiemannSolver", argc, argv);

  // set up the test arrays
  const uint_fast32_t num_test = 100000;
  std::vector< double > rho(num_test), P(num_test);
  std::vector< CoordinateVector<> > u(num_test), normals(num_test),
      vface(num_test);
  for (uint_fast32_t i = 0; i < num_test; ++i) {
    // densities in the range [0.125, 1.[
    rho[i] = 0.125 + Utilities::random_double() * 0.875;
    // velocities in the range [-1., 1.[
    u[i][0] = 2. * Utilities::random_double() - 1.;
    u[i][1] = 2. * Utilities::random_double() - 1.;
    u[i][2] = 2. * Utilities::random_double() - 1.;
    // pressures in the range [0.1, 1.[
    P[i] = 0.1 + Utilities::random_double() * 0.9;
    // normals
    normals[i][0] = Utilities::random_double();
    normals[i][1] = Utilities::random_double();
    normals[i][2] = Utilities::random_double();
    normals[i] /= normals[i].norm();
    // frame velocities in the range [-1., 1.[
    vface[i][0] = 2. * Utilities::random_double() - 1.;
    vface[i][1] = 2. * Utilities::random_double() - 1.;
    vface[i][2] = 2. * Utilities::random_double() - 1.;
  }

  ExactRiemannSolver exact_solver(5. / 3.);
  HLLCRiemannSolver hllc_solver(5. / 3.);
  double mflux = 0.;
  CoordinateVector<> pflux;
  double Eflux = 0.;

  timingtools_start_timing_block("ExactRiemannSolver") {
    timingtools_start_timing();
    for (uint_fast32_t i = 0; i < num_test; ++i) {
      const uint_fast32_t iplus = (i + 1) % num_test;
      exact_solver.solve_for_flux(rho[i], u[i], P[i], rho[iplus], u[iplus],
                                  P[iplus], mflux, pflux, Eflux, normals[i],
                                  vface[i]);
    }
    timingtools_stop_timing();
  }
  timingtools_end_timing_block("ExactRiemannSolver");

  timingtools_start_timing_block("HLLCRiemannSolver") {
    timingtools_start_timing();
    for (uint_fast32_t i = 0; i < num_test; ++i) {
      const uint_fast32_t iplus = (i + 1) % num_test;
      hllc_solver.solve_for_flux(rho[i], u[i], P[i], rho[iplus], u[iplus],
                                 P[iplus], mflux, pflux, Eflux, normals[i],
                                 vface[i]);
    }
    timingtools_stop_timing();
  }
  timingtools_end_timing_block("HLLCRiemannSolver");
}
