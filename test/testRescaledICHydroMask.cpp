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
 * @file testRescaledICHydroMask.cpp
 *
 * @brief Unit test for the RescaledICHydroMask class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CartesianDensityGrid.hpp"
#include "DiscICDensityFunction.hpp"
#include "HydroIntegrator.hpp"
#include "RescaledICHydroMask.hpp"
#include <fstream>

/**
 * @brief Unit test for the RescaledICHydroMask class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Box<> box(CoordinateVector<>(-1.), CoordinateVector<>(3.));
  DiscICDensityFunction density_function(1.e16, 500., 1., 1.5, 1., 0.5);

  CartesianDensityGrid grid(box, 64, false, true);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  const double gamma = 5. / 3.;

  HydroIntegrator integrator(gamma, false, false, 0.2);
  integrator.initialize_hydro_variables(grid);

  RescaledICHydroMask hydro_mask(CoordinateVector<>(0.), 0.5, 0.01, 1., 0.02);
  hydro_mask.initialize_mask(grid);
  hydro_mask.apply_mask(grid);

  std::ofstream ofile("test_rescaled_hydro_mask.txt");
  ofile << "# x\ty\tz\trho\tux\tuy\tuz\tp\n";
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    const CoordinateVector<> p = it.get_cell_midpoint();
    const double rho = it.get_hydro_variables().get_primitives_density();
    const CoordinateVector<> v =
        it.get_hydro_variables().get_primitives_velocity();
    const double P = it.get_hydro_variables().get_primitives_pressure();
    ofile << p.x() << "\t" << p.y() << "\t" << p.z() << "\t" << rho << "\t"
          << v.x() << "\t" << v.y() << "\t" << v.z() << "\t" << P << "\n";
  }

  return 0;
}
