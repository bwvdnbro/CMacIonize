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
 * @file testDensityPDFCalculator.cpp
 *
 * @brief Unit test for the DensityPDFCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "DensityPDFCalculator.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Unit test for the DensityPDFCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double box[6] = {0., 0., 0., 1., 1., 1.};
  CoordinateVector< int_fast32_t > ncell(32, 32, 32);
  HydroDensitySubGrid subgrid(box, ncell);

  RandomGenerator random_generator;
  for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
       ++cellit) {
    const double u0 = random_generator.get_uniform_random_double();
    const double u1 = random_generator.get_uniform_random_double();
    const double rho =
        std::pow(10., std::sqrt(-2. * std::log(u0)) * std::cos(2. * M_PI * u1));
    cellit.get_hydro_variables().set_primitives_density(rho);
  }

  DensityPDFCalculator calculator(1, 1.e-3, 1.e3, 100);
  calculator.calculate_density_PDF(0, subgrid);
  calculator.output("test_densityPDFcalculator.txt");

  return 0;
}
