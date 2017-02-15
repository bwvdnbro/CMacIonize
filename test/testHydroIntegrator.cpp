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
 * @file testHydroIntegrator.cpp
 *
 * @brief Unit test for the HydroIntegrator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CartesianDensityGrid.hpp"
#include "DensityFunction.hpp"
#include "HydroIntegrator.hpp"
#include <fstream>

/**
 * @brief DensityFunction implementation that sets up a basic Sod shock problem.
 */
class SodShockDensityFunction : public DensityFunction {
public:
  /**
   * @brief Get the density field at the given position.
   *
   * @param position Position.
   * @return DensityValues at that position.
   */
  virtual DensityValues operator()(CoordinateVector<> position) const {
    DensityValues values;
    double r = (position - CoordinateVector<>(0.5)).norm();
    if (r < 0.25) {
      values.set_number_density(1.);
      values.set_temperature(1.);
    } else {
      values.set_number_density(0.125);
      values.set_temperature(0.1);
    }
    return values;
  }
};

/**
 * @brief Unit test for the HydroIntegrator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  HydroIntegrator integrator(5. / 3.);

  Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CoordinateVector< int > ncell(32);
  SodShockDensityFunction density_function;
  CoordinateVector< bool > periodic(true);
  CartesianDensityGrid grid(box, ncell, density_function, periodic, true);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  integrator.initialize_hydro_variables(grid);

  // write initial snapshot
  {
    std::ofstream snapfile("hydro_snap_0.txt");
    double mtot = 0.;
    double etot = 0.;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      double r = (it.get_cell_midpoint() - CoordinateVector<>(0.5)).norm();
      snapfile << r << "\t" << it.get_hydro_primitive_density() << "\t"
               << it.get_hydro_primitive_velocity_x() << "\t"
               << it.get_hydro_primitive_pressure() << "\n";
      mtot += it.get_hydro_conserved_mass();
      etot += it.get_hydro_conserved_total_energy();
    }
    cmac_status("Total mass: %g, total energy: %g", mtot, etot);
  }

  for (unsigned int i = 0; i < 100; ++i) {
    integrator.do_hydro_step(grid, 0.001);
  }

  // write final snapshot
  {
    std::ofstream snapfile("hydro_snap_1.txt");
    double mtot = 0.;
    double etot = 0.;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      double r = (it.get_cell_midpoint() - CoordinateVector<>(0.5)).norm();
      snapfile << r << "\t" << it.get_hydro_primitive_density() << "\t"
               << it.get_hydro_primitive_velocity_x() << "\t"
               << it.get_hydro_primitive_pressure() << "\n";
      mtot += it.get_hydro_conserved_mass();
      etot += it.get_hydro_conserved_total_energy();
    }
    cmac_status("Total mass: %g, total energy: %g", mtot, etot);
  }

  return 0;
}
