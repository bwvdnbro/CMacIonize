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
 * @file testSurfaceDensityCalculator.cpp
 *
 * @brief Unit test for the SurfaceDensityCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "SurfaceDensityCalculator.hpp"

/**
 * @brief Unit test for the SurfaceDensityCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double box[6] = {0., 0., 0., 1., 1., 1.};
  CoordinateVector< int_fast32_t > ncell(32, 32, 32);
  HydroDensitySubGrid subgrid(box, ncell);

  for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
       ++cellit) {
    const double r =
        (cellit.get_cell_midpoint() - CoordinateVector<>(0.5)).norm();
    cellit.get_hydro_variables().set_primitives_density(1. / (1. + r));
  }

  SurfaceDensityCalculator calculator(1, 32);
  calculator.calculate_surface_density(0, subgrid);
  calculator.output("test_surfacedensitycalculator.txt", Box<>(0., 1.));

  return 0;
}
