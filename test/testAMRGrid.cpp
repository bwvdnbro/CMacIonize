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
 * @file testAMRGrid.cpp
 *
 * @brief Unit test for the AMRGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AMRGrid.hpp"
#include "Assert.hpp"
#include "Box.hpp"
#include "CoordinateVector.hpp"
#include <fstream>

/**
 * @brief Unit test for the AMRGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  AMRGrid< double > grid(
      Box(CoordinateVector<>(), CoordinateVector<>(2., 1., 1.)),
      CoordinateVector< int >(2, 1, 1));

  grid.create_cell(2, CoordinateVector<>(0.1));
  grid.create_cell(2, CoordinateVector<>(1.1, 0.1, 0.1));

  std::ofstream gfile("amrgrid.dump");
  grid.print(gfile);

  unsigned long key = grid.get_key(2, CoordinateVector<>(0.1));
  assert_condition(key == 64);
  key = grid.get_key(2, CoordinateVector<>(1.1, 0.1, 0.1));
  assert_condition(key == 0x0010000000000040);

  grid.get_cell(CoordinateVector<>(0.1)) = 42.;
  assert_condition(grid[64] == 42.);
  grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) = 3.14;
  assert_condition(grid[0x0010000000000040] == 3.14);

  grid[64] = 3.14;
  assert_condition(grid.get_cell(CoordinateVector<>(0.1)) == 3.14);
  grid[0x0010000000000040] = 42.;
  assert_condition(grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) == 42.);

  key = grid.get_key(2, CoordinateVector<>(0.7));
  assert_condition(key == 71);
  grid.create_cell(key) = 5.5;
  assert_condition(grid.get_cell(CoordinateVector<>(0.7)) == 5.5);

  grid.create_all_cells(3);
  assert_condition(grid.get_number_of_cells() == 2 * 8 * 8 * 8);

  return 0;
}
