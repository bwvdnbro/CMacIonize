/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testTreeSelfGravity.cpp
 *
 * @brief Unit test for the TreeSelfGravity class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CartesianDensityGrid.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "TreeSelfGravity.hpp"
#include <fstream>

/**
 * @brief Unit test for the TreeSelfGravity class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  HomogeneousDensityFunction density_function(1.);
  density_function.initialize();

  CartesianDensityGrid grid(box, 32, false, true);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  // initialize hydro variables
  const CoordinateVector<> center(0.5);
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    const double r = (it.get_cell_midpoint() - center).norm();
    if (r < 0.5) {
      it.get_hydro_variables().set_conserved_mass(1. * it.get_volume());
    } else {
      it.get_hydro_variables().set_conserved_mass(0.);
    }
  }

  TreeSelfGravity self_gravity(grid, 0.25);

  cmac_status("Computing accelerations...");
  self_gravity.compute_accelerations(grid);
  cmac_status("Done.");

  std::ofstream ofile("test_treeselfgravity.txt");
  ofile << "# r (m)\ta (m s^-2)\tx (m)\ty (m)\tz (m)\tax (m s^-2)\tay (m "
           "s^-2)\taz (m s^-2)\n";
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    const CoordinateVector<> r = it.get_cell_midpoint() - center;
    const double rnorm = r.norm();
    const CoordinateVector<> a =
        it.get_hydro_variables().get_gravitational_acceleration();
    const double anorm = a.norm();
    ofile << rnorm << "\t" << anorm << "\t" << r.x() << "\t" << r.y() << "\t"
          << r.z() << "\t" << a.x() << "\t" << a.y() << "\t" << a.z() << "\n";
  }

  return 0;
}
