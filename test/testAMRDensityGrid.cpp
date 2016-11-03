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
 * @file testAMRDensityGrid.cpp
 *
 * @brief Unit test for the AMRDensityGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AMRDensityGrid.hpp"
#include "Assert.hpp"
#include "TerminalLog.hpp"

/**
 * @brief Unit test for the AMRDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  TerminalLog log(LOGLEVEL_INFO);
  AMRDensityGrid grid(Box(CoordinateVector<>(0.), CoordinateVector<>(1.)), 32,
                      false, &log);

  assert_condition(grid.get_number_of_cells() == 32 * 32 * 32);
  // 32768 is the grid in the lower left front corner, which has indices 000 on
  // all levels
  // the 32768 bit is set to indicate its level: 5
  assert_condition(grid.get_cell_volume(32768) ==
                   (1. / 32) * (1. / 32) * (1. / 32));
  unsigned long key = grid.get_cell_index(CoordinateVector<>(0.01));
  assert_condition(key == 32768);
  CoordinateVector<> midpoint = grid.get_cell_midpoint(32768);
  assert_condition(midpoint.x() == 0.015625);
  assert_condition(midpoint.y() == 0.015625);
  assert_condition(midpoint.z() == 0.015625);

  //  unsigned int ncell = 0;
  //  for(auto it = grid.begin(); it != grid.end(); ++it){
  //    ++ncell;
  //  }
  //  assert_condition(ncell == grid.get_number_of_cells());

  return 0;
}
