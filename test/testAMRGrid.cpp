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
  // create an empty grid
  AMRGrid< double > grid(
      Box(CoordinateVector<>(), CoordinateVector<>(2., 1., 1.)),
      CoordinateVector< int >(2, 1, 1));

  /// basic functionality

  // create two level 2 cells, one in each block
  grid.create_cell(2, CoordinateVector<>(0.1));
  grid.create_cell(2, CoordinateVector<>(1.1, 0.1, 0.1));

  // print the grid out for visual inspection
  std::ofstream gfile("amrgrid.dump");
  grid.print(gfile);

  // test the methods to generate a key from a level and a position
  unsigned long key = grid.get_key(2, CoordinateVector<>(0.1));
  assert_condition(key == 64);
  key = grid.get_key(2, CoordinateVector<>(1.1, 0.1, 0.1));
  assert_condition(key == 0x0010000000000040);

  // test the methods to get the key for the deepest cell containing a position
  key = grid.get_key(CoordinateVector<>(0.1));
  assert_condition(key == 64);
  key = grid.get_key(CoordinateVector<>(1.1, 0.1, 0.1));
  assert_condition(key == 0x0010000000000040);

  // test the access on position and access on key routines, part 1
  grid.get_cell(CoordinateVector<>(0.1)) = 42.;
  assert_condition(grid[64] == 42.);
  grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) = 3.14;
  assert_condition(grid[0x0010000000000040] == 3.14);

  // part 2
  grid[64] = 3.14;
  assert_condition(grid.get_cell(CoordinateVector<>(0.1)) == 3.14);
  grid[0x0010000000000040] = 42.;
  assert_condition(grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) == 42.);

  // test volume, midpoint and geometry routines
  assert_condition(grid.get_volume(64) == 0.015625);
  CoordinateVector<> midpoint = grid.get_midpoint(64);
  assert_condition(midpoint.x() == 0.125);
  assert_condition(midpoint.y() == 0.125);
  assert_condition(midpoint.z() == 0.125);
  Box geometry = grid.get_geometry(64);
  assert_condition(geometry.get_anchor().x() == 0.);
  assert_condition(geometry.get_anchor().y() == 0.);
  assert_condition(geometry.get_anchor().z() == 0.);
  assert_condition(geometry.get_sides().x() == 0.25);
  assert_condition(geometry.get_sides().y() == 0.25);
  assert_condition(geometry.get_sides().z() == 0.25);

  // create a new cell and immediately access it
  key = grid.get_key(2, CoordinateVector<>(0.7));
  assert_condition(key == 71);
  grid.create_cell(key) = 5.5;
  assert_condition(grid.get_cell(CoordinateVector<>(0.7)) == 5.5);

  /// advanced functionality

  // create all level 3 cells
  // we now have a Cartesian grid of 16x8x8
  grid.create_all_cells(3);
  assert_condition(grid.get_number_of_cells() == 2 * 8 * 8 * 8);

  // basic indexing: check first key and next key routines
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

  // advanced indexing: check neighbour finding routines

  // within a parent cell
  next_key = grid.get_neighbour(512, CoordinateVector< char >(0, 1, 0),
                                CoordinateVector<>(0.1, 0.125, 0.1));
  // 640: 1*512 + (0*256 + 1*128 + 0*64) + (0*32 + 0*16 + 0*8) +
  //              (0*4 + 0*2 + 0*1)
  assert_condition(next_key == 640);
  next_key = grid.get_neighbour(512, CoordinateVector< char >(1, 0, 0),
                                CoordinateVector<>(0.125, 0.1, 0.1));
  // 768: 1*512 + (1*256 + 0*128 + 0*64) + (0*32 + 0*16 + 0*8) +
  //              (0*4 + 0*2 + 0*1)
  assert_condition(next_key == 768);
  next_key = grid.get_neighbour(512, CoordinateVector< char >(0, 0, 1),
                                CoordinateVector<>(0.1, 0.1, 0.125));
  // 576: 1*512 + (0*256 + 0*128 + 1*64) + (0*32 + 0*16 + 0*8) +
  //              (0*4 + 0*2 + 0*1)
  assert_condition(next_key == 576);
  next_key = grid.get_neighbour(512, CoordinateVector< char >(0, 0, -1),
                                CoordinateVector<>(0.1, 0.1, 0.));
  assert_condition(next_key == AMRGRID_MAXKEY);
  next_key = grid.get_neighbour(512, CoordinateVector< char >(0, 1, 1),
                                CoordinateVector<>(0.1, 0.125, 0.125));
  // 704: 1*512 + (0*256 + 1*128 + 1*64) + (0*32 + 0*16 + 0*8) +
  //              (0*4 + 0*2 + 0*1)
  assert_condition(next_key == 704);

  // between neighbouring parent cells
  next_key = grid.get_neighbour(704, CoordinateVector< char >(0, 0, 1),
                                CoordinateVector<>(0.1, 0.2, 0.25));
  // 648: 1*512 + (0*256 + 1*128 + 0*64) + (0*32 + 0*16 + 1*8) +
  //              (0*4 + 0*2 + 0*1)
  assert_condition(next_key == 648);
  next_key = grid.get_neighbour(648, CoordinateVector< char >(0, 0, -1),
                                CoordinateVector<>(0.1, 0.2, 0.25));
  assert_condition(next_key == 704);

  // between neighbouring blocks
  next_key = grid.get_neighbour(1023, CoordinateVector< char >(1, 0, 0),
                                CoordinateVector<>(1., 0.9, 0.9));
  // 731: 1*512 + (0*256 + 1*128 + 1*64) + (0*32 + 1*16 + 1*8) +
  //              (0*4 + 1*2 + 1*1) (=0x2db)
  assert_condition(next_key == 0x00100000000002db);
  next_key =
      grid.get_neighbour(0x00100000000002db, CoordinateVector< char >(-1, 0, 0),
                         CoordinateVector<>(1., 0.9, 0.9));
  assert_condition(next_key == 1023);

  // check that we can walk through the grid in a straight line

  // x direction
  key = grid.get_first_key();
  unsigned int ncell = 0;
  while (key != AMRGRID_MAXKEY) {
    key = grid.get_neighbour(key, CoordinateVector< char >(1, 0, 0),
                             CoordinateVector<>(0.1));
    ++ncell;
  }
  assert_condition(ncell == 16);

  // y direction
  key = grid.get_first_key();
  ncell = 0;
  while (key != AMRGRID_MAXKEY) {
    key = grid.get_neighbour(key, CoordinateVector< char >(0, 1, 0),
                             CoordinateVector<>(0.1));
    ++ncell;
  }
  assert_condition(ncell == 8);

  // z direction
  key = grid.get_first_key();
  ncell = 0;
  while (key != AMRGRID_MAXKEY) {
    key = grid.get_neighbour(key, CoordinateVector< char >(0, 0, 1),
                             CoordinateVector<>(0.1));
    ++ncell;
  }
  assert_condition(ncell == 8);

  // diagonal in y and z (why not?)
  key = grid.get_first_key();
  ncell = 0;
  while (key != AMRGRID_MAXKEY) {
    key = grid.get_neighbour(key, CoordinateVector< char >(0, 1, 1),
                             CoordinateVector<>(0.1));
    ++ncell;
  }
  assert_condition(ncell == 8);

  // advanced indexing: first cell in a given direction
  key = grid.get_first_key(CoordinateVector< char >(1, 0, 0),
                           CoordinateVector<>(0.1));
  assert_condition(key == 512);
  key = grid.get_first_key(CoordinateVector< char >(0, 1, 0),
                           CoordinateVector<>(0.1));
  assert_condition(key == 512);
  key = grid.get_first_key(CoordinateVector< char >(0, -1, 0),
                           CoordinateVector<>(0.1, 0.9, 0.1));
  // 658: 1*512 + (0*256 + 1*128 + 0*64) + (0*32 + 1*16 + 0*8) +
  //              (0*4 + 1*2 + 0*1)
  assert_condition(key == 658);
  key = grid.get_first_key(CoordinateVector< char >(-1, -1, 0),
                           CoordinateVector<>(1.9, 0.9, 0.1));
  // 950: 1*512 + (1*256 + 1*128 + 0*64) + (1*32 + 1*16 + 0*8) +
  //              (1*4 + 1*2 + 0*1) = (0x3b6)
  assert_condition(key == 0x00100000000003b6);

  // check Morton iteration
  key = grid.get_first_key();
  ncell = 0;
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

  /// mesh refinement

  // refine a single cell
  grid.refine_cell(703);
  // repeat the Morton iteration test
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
