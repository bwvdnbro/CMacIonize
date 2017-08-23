/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensityGrid.cpp
 *
 * @brief DensityGrid implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DensityGrid.hpp"
#include "DensityGridTraversalJobMarket.hpp"

/**
 * @brief Initialize the cells in the grid.
 *
 * All implementations should call this method in their initialization()
 * routine.
 *
 * @param block Continuous block of indices to initialize.
 * @param function DensityFunction that sets the density.
 * @param worksize Number of parallel threads to use. If a negative number is
 * given, all available threads will be used.
 */
void DensityGrid::set_densities(
    std::pair< unsigned long, unsigned long > &block, DensityFunction &function,
    int worksize) {
  DensityGridInitializationFunction init(function, _has_hydro);
  WorkDistributor<
      DensityGridTraversalJobMarket< DensityGridInitializationFunction >,
      DensityGridTraversalJob< DensityGridInitializationFunction > >
      workers(worksize);

  if (_log) {
    _log->write_status("Initializing grid using ",
                       workers.get_worksize_string(), ".");
  }

  DensityGridTraversalJobMarket< DensityGridInitializationFunction > jobs(
      *this, init, block);
  workers.do_in_parallel(jobs);

  if (_log) {
    _log->write_status("Done initializing grid.");
  }
}
