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

  key = grid.get_key(CoordinateVector<>(0.1));
  assert_condition(key == 64);
  key = grid.get_key(CoordinateVector<>(1.1, 0.1, 0.1));
  assert_condition(key == 0x0010000000000040);

  grid.get_cell(CoordinateVector<>(0.1)) = 42.;
  assert_condition(grid[64] == 42.);
  grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) = 3.14;
  assert_condition(grid[0x0010000000000040] == 3.14);

  grid[64] = 3.14;
  assert_condition(grid.get_cell(CoordinateVector<>(0.1)) == 3.14);
  grid[0x0010000000000040] = 42.;
  assert_condition(grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) == 42.);

  assert_condition(grid.get_volume(64) == 0.015625);
  CoordinateVector<> midpoint = grid.get_midpoint(64);
  assert_condition(midpoint.x() == 0.125);
  assert_condition(midpoint.y() == 0.125);
  assert_condition(midpoint.z() == 0.125);

  key = grid.get_key(2, CoordinateVector<>(0.7));
  assert_condition(key == 71);
  grid.create_cell(key) = 5.5;
  assert_condition(grid.get_cell(CoordinateVector<>(0.7)) == 5.5);

  grid.create_all_cells(3);
  assert_condition(grid.get_number_of_cells() == 2 * 8 * 8 * 8);

  assert_condition(grid.get_first_key() == 512);
  unsigned long next_key = grid.get_next_key(512);
  // 576: 1*512 + (0*256 + 0*128 + 1*64) + (0*32 + 0*16 + 0*8) +
  //              (0*4 + 0*2 + 0*1)
  assert_condition(next_key == 576);
  // 1023: the above with all ones
  next_key = grid.get_next_key(1023);
  // the number below is 512 in the second block
  assert_condition(next_key == 0x0010000000000200);
  // the number below is 1023 in the second block
  next_key = grid.get_next_key(0x00100000000003ff);
  assert_condition(next_key == grid.get_max_key());

  key = grid.get_first_key();
  unsigned int ncell = 0;
  // we need a way to make sure all cells are traversed (exactly once)
  // to this end, we calculate the sum of all keys
  unsigned long keysum = 0;
  while (key != grid.get_max_key()) {
    keysum += key;
    key = grid.get_next_key(key);
    ++ncell;
  }
  assert_condition(ncell == grid.get_number_of_cells());
  // the sum should be the sum of all number between 512 and 1023
  // + all numbers between 0x0010000000000200 and 0x00100000000003ff
  unsigned long refsum =
      256 * (512 + 1023) + 256 * (0x0010000000000200 + 0x00100000000003ff);
  assert_condition(keysum == refsum);

  // refine a cell
  grid.refine_cell(703);
  // count again
  assert_condition(grid.get_number_of_cells() == 2 * 8 * 8 * 8 + 7);
  key = grid.get_first_key();
  ncell = 0;
  keysum = 0;
  while (key != grid.get_max_key()) {
    keysum += key;
    key = grid.get_next_key(key);
    ++ncell;
  }
  assert_condition(ncell == grid.get_number_of_cells());
  // the sum should now be the same as before, but without the cell we refined
  // (703), and with the 8 new cells added
  // the key for each of these cells consists of the part of 703 without the
  // level bit (512): 191, and a new level bit (4096)
  // to find the new part of the key, we just alternate all possible
  // combinations of 2048, 1024, and 512
  refsum = 256 * (512 + 1023) +
           256 * (0x0010000000000200 + 0x00100000000003ff) - 703 +
           8 * (4096 + 191) + 4 * 2048 + 4 * 1024 + 4 * 512;
  assert_condition(keysum == refsum);

  return 0;
}
