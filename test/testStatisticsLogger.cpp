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
 * @file testStatisticsLogger.cpp
 *
 * @brief Unit test for the StatisticsLogger class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "CartesianDensityGrid.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "HydroIntegrator.hpp"
#include "StatisticsLogger.hpp"

/**
 * @brief Unit test for the StatisticsLogger class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  HomogeneousDensityFunction density_function(1., 500., 1.);
  density_function.initialize();
  CartesianDensityGrid grid(Box<>(0., 1.), 10, false, true);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  HydroIntegrator integrator(5. / 3., true, true, 0.2);
  integrator.initialize_hydro_variables(grid);

  {
    StatisticsLogger logger;
    logger.write_statistics(0., grid);

    {
      RestartWriter restart_writer("statistics_logger.dump");
      logger.write_restart_file(restart_writer);
    }

    logger.write_statistics(1., grid);
  }

  /// restart test
  {
    RestartReader restart_reader("statistics_logger.dump");
    StatisticsLogger logger(restart_reader);
    logger.write_statistics(2., grid);
  }

  return 0;
}
