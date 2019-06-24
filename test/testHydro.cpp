/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testHydro.cpp
 *
 * @brief Unit test for the Hydro class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Hydro.hpp"

#include <fstream>

/**
 * @brief Unit test for the Hydro class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Hydro hydro(5. / 3., 100., 1.e4, 1.e99);

  HydroVariables variables[100];
  IonizationVariables ionvariables[100];
  double limiters[1000];

  for (uint_fast32_t i = 0; i < 100; ++i) {
    if (i < 50) {
      variables[i].set_primitives_density(1.);
      variables[i].set_primitives_pressure(1.);
    } else {
      variables[i].set_primitives_density(0.125);
      variables[i].set_primitives_pressure(0.1);
    }
    for (uint_fast32_t j = 0; j < 5; ++j) {
      limiters[10 * i + 2 * j] = DBL_MAX;
      limiters[10 * i + 2 * j + 1] = -DBL_MAX;
    }
    ionvariables[i].set_ionic_fraction(ION_H_n, 1.);
  }

  for (uint_fast32_t i = 0; i < 100; ++i) {
    hydro.set_conserved_variables(variables[i], 0.01);
  }

  const double dt = 0.001;
  const CoordinateVector<> dx(0.01, 0., 0.);
  for (uint_fast32_t istep = 0; istep < 100; ++istep) {

    for (uint_fast32_t i = 0; i < 100; ++i) {
      uint_fast32_t inext = (i + 1) % 100;
      hydro.do_gradient_calculation(0, variables[i], variables[inext], 100.,
                                    &limiters[10 * i], &limiters[10 * inext]);
    }

    for (uint_fast32_t i = 0; i < 100; ++i) {
      hydro.apply_slope_limiter(variables[i], &limiters[10 * i], dx);
    }

    for (uint_fast32_t i = 0; i < 100; ++i) {
      uint_fast32_t inext = (i + 1) % 100;
      hydro.do_flux_calculation(0, variables[i], variables[inext], 0.01, 1.);
    }

    for (uint_fast32_t i = 0; i < 100; ++i) {
      for (uint_fast32_t j = 0; j < 5; ++j) {
        variables[i].conserved(j) += variables[i].delta_conserved(j) * dt;
      }
    }

    for (uint_fast32_t i = 0; i < 100; ++i) {
      hydro.set_primitive_variables(variables[i], ionvariables[i], 100.);
    }

    for (uint_fast32_t i = 0; i < 100; ++i) {
      for (uint_fast32_t j = 0; j < 5; ++j) {
        variables[i].delta_conserved(j) = 0.;
        variables[i].primitive_gradients(j) = CoordinateVector<>(0.);
        limiters[10 * i + 2 * j] = DBL_MAX;
        limiters[10 * i + 2 * j + 1] = -DBL_MAX;
      }
    }
  }

  std::ofstream ofile("testHydro_result.txt");
  ofile << "# x\trho\tvx\tP\n";
  for (uint_fast32_t i = 0; i < 100; ++i) {
    ofile << (i + 0.5) * 0.01 << "\t" << variables[i].get_primitives_density()
          << "\t" << variables[i].get_primitives_velocity().x() << "\t"
          << variables[i].get_primitives_pressure() << "\n";
  }

  return 0;
}
