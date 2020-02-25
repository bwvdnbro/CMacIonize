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
 * @file testSingleSupernovaPhotonSourceDistribution.cpp
 *
 * @brief Unit test for the SingleSupernovaPhotonSourceDistribution class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "SingleSupernovaPhotonSourceDistribution.hpp"
#include <fstream>

/**
 * @brief Unit test for the SingleSupernovaPhotonSourceDistribution class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double kpc_in_m = 3.086e19;
  const double Myr_in_s = 1.e6 * 365.25 * 24. * 3600.;

  CoordinateVector<> position;
  SingleSupernovaPhotonSourceDistribution distribution(position, 10. * Myr_in_s,
                                                       1.e49, 1.e44);

  // restart test
  {
    RestartWriter restart_writer(
        "singlesupernovaphotonsourcedistribution.dump");
    distribution.write_restart_file(restart_writer);
  }

  HomogeneousDensityFunction testfunction(1., 2000.);
  testfunction.initialize();
  CartesianDensityGrid grid(Box<>(-1.5 * kpc_in_m, 3. * kpc_in_m), 32, false,
                            true);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, testfunction);

  distribution.add_stellar_feedback(grid, 0.);
  DensityGrid::iterator cell =
      static_cast< DensityGrid * >(&grid)->get_cell(position);
  assert_condition(cell.get_hydro_variables().get_energy_term() == 0.);
  distribution.add_stellar_feedback(grid, 10. * Myr_in_s);
  assert_condition(cell.get_hydro_variables().get_energy_term() == 1.e44);

  // restart test
  {
    RestartReader restart_reader(
        "singlesupernovaphotonsourcedistribution.dump");
    SingleSupernovaPhotonSourceDistribution restart_distribution(
        restart_reader);
  }

  return 0;
}
