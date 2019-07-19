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
#include <cinttypes>
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
      Box<>(CoordinateVector<>(), CoordinateVector<>(2., 1., 1.)),
      CoordinateVector< uint_fast32_t >(2, 1, 1));

  /// basic functionality

  // create two level 2 cells, one in each block
  grid.create_cell(2, CoordinateVector<>(0.1));
  grid.create_cell(2, CoordinateVector<>(1.1, 0.1, 0.1));

  // print the grid out for visual inspection
  std::ofstream gfile("amrgrid.dump");
  grid.print(gfile);

  // test the methods to generate a key from a level and a position
  uint64_t key = grid.get_key(2, CoordinateVector<>(0.1));
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
  assert_condition(grid[64].value() == 42.);
  grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) = 3.14;
  assert_condition(grid[0x0010000000000040].value() == 3.14);

  // part 2
  grid[64].value() = 3.14;
  assert_condition(grid.get_cell(CoordinateVector<>(0.1)) == 3.14);
  grid[0x0010000000000040].value() = 42.;
  assert_condition(grid.get_cell(CoordinateVector<>(1.1, 0.1, 0.1)) == 42.);

  // test volume, midpoint and geometry routines
  assert_condition(grid[64].get_volume() == 0.015625);
  CoordinateVector<> midpoint = grid[64].get_midpoint();
  assert_condition(midpoint.x() == 0.125);
  assert_condition(midpoint.y() == 0.125);
  assert_condition(midpoint.z() == 0.125);
  Box<> geometry = grid[64].get_geometry();
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

  // check the level
  assert_condition(grid[71].get_level() == 2);
  assert_condition(grid[71].get_parent()->get_level() == 1);
  assert_condition(grid[71].get_parent()->get_parent()->get_level() == 0);
  assert_condition(grid[71].get_parent()->get_parent()->get_parent() ==
                   nullptr);

  /// advanced functionality

  // create all level 3 cells
  // we now have a Cartesian grid of 16x8x8
  grid.create_all_cells(3);
  assert_condition(grid.get_number_of_cells() == 2 * 8 * 8 * 8);

  // basic indexing: check first key and next key routines
  assert_condition(grid.get_first_key() == 512);
  uint64_t next_key = grid.get_next_key(512);
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

  // check Morton iteration
  key = grid.get_first_key();
  uint_fast32_t ncell = 0;
  // we need a way to make sure all cells are traversed (exactly once)
  // to this end, we calculate the sum of all keys
  uint64_t keysum = 0;
  while (key != grid.get_max_key()) {
    keysum += key;
    key = grid.get_next_key(key);
    ++ncell;
  }
  assert_condition(ncell == grid.get_number_of_cells());
  // the sum should be the sum of all numbers between 512 and 1023
  // + all numbers between 0x0010000000000200 and 0x00100000000003ff
  uint64_t refsum =
      256 * (512 + 1023) + 256 * (0x0010000000000200 + 0x00100000000003ff);
  assert_condition(keysum == refsum);

  /// mesh refinement

  // refine a single cell
  uint64_t refined_key = grid.refine_cell(703);
  assert_condition(refined_key == 4287);
  // repeat the Morton iteration test
  assert_condition(grid.get_number_of_cells() == 2 * 8 * 8 * 8 + 7);
  key = grid.get_first_key();
  ncell = 0;
  keysum = 0;
  while (key != grid.get_max_key()) {
    keysum += key;
    grid[key].value() = key;
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

  grid.set_ngbs(CoordinateVector< bool >(false));

  // first check an easy one:
  // 512 = 1 000 000 000
  assert_condition(grid[512].get_ngb(AMRNGBPOSITION_LEFT) == nullptr);
  // 768 = 1 100 000 000
  assert_condition(grid[512].get_ngb(AMRNGBPOSITION_RIGHT)->value() == 768.);
  assert_condition(grid[512].get_ngb(AMRNGBPOSITION_FRONT) == nullptr);
  // 640 = 1 010 000 000
  assert_condition(grid[512].get_ngb(AMRNGBPOSITION_BACK)->value() == 640.);
  assert_condition(grid[512].get_ngb(AMRNGBPOSITION_BOTTOM) == nullptr);
  // 576 = 1 001 000 000
  assert_condition(grid[512].get_ngb(AMRNGBPOSITION_TOP)->value() == 576.);

  // now a couple arbitrary ones in the middle of the first block
  uint_fast64_t tests[8][7] = {// 568 = 1 000 111 000
                               // 792 = 1 100 011 000
                               // 824 = 1 100 111 000
                               // 680 = 1 010 101 000
                               // 696 = 1 010 111 000
                               // 624 = 1 001 110 000
                               // 632 = 1 001 111 000
                               {568, 792, 824, 680, 696, 624, 632},
                               // 632 = 1 001 111 000
                               // 856 = 1 101 011 000
                               // 888 = 1 101 111 000
                               // 744 = 1 011 101 000
                               // 760 = 1 011 111 000
                               // 568 = 1 000 111 000
                               // 561 = 1 000 110 001
                               {632, 856, 888, 744, 760, 568, 561},
                               // 696 = 1 010 111 000
                               // 920 = 1 110 011 000
                               // 952 = 1 110 111 000
                               // 568 = 1 000 111 000
                               // 554 = 1 000 101 010
                               // 752 = 1 011 110 000
                               // 760 = 1 011 111 000
                               {696, 920, 952, 568, 554, 752, 760},
                               // 760 = 1 011 111 000
                               // 984 = 1 111 011 000
                               // 1016 = 1 111 111 000
                               // 632 = 1 001 111 000
                               // 618 = 1 001 101 010
                               // 696 = 1 010 111 000
                               // 689 = 1 010 110 001
                               {760, 984, 1016, 632, 618, 696, 689},
                               // 824 = 1 100 111 000
                               // 568 = 1 000 111 000
                               // 540 = 1 000 011 100
                               // 936 = 1 110 101 000
                               // 952 = 1 110 111 000
                               // 880 = 1 101 110 000
                               // 888 = 1 101 111 000
                               {824, 568, 540, 936, 952, 880, 888},
                               // 888 = 1 101 111 000
                               // 632 = 1 001 111 000
                               // 604 = 1 001 011 100
                               // 1000 = 1 111 101 000
                               // 1016 = 1 111 111 000
                               // 824 = 1 100 111 000
                               // 817 = 1 100 110 001
                               {888, 632, 604, 1000, 1016, 824, 817},
                               // 952 = 1 110 111 000
                               // 696 = 1 010 111 000
                               // 668 = 1 010 011 100
                               // 824 = 1 100 111 000
                               // 810 = 1 100 101 010
                               // 1008 = 1 111 110 000
                               // 1016 = 1 111 111 000
                               {952, 696, 668, 824, 810, 1008, 1016},
                               // 1016 = 1 111 111 000
                               // 760 = 1 011 111 000
                               // 732 = 1 011 011 100
                               // 888 = 1 101 111 000
                               // 874 = 1 101 101 010
                               // 952 = 1 110 111 000
                               // 945 = 1 110 110 001
                               {1016, 760, 732, 888, 874, 952, 945}};
  for (uint_fast8_t i = 0; i < 8; ++i) {
    cmac_status("Testing cell %" PRIuFAST64, tests[i][0]);

    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][0]].get_ngb(AMRNGBPOSITION_LEFT)->value(),
                tests[i][1]);
    assert_condition(grid[tests[i][0]].get_ngb(AMRNGBPOSITION_LEFT)->value() ==
                     tests[i][1]);
    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][1]].get_ngb(AMRNGBPOSITION_RIGHT)->value(),
                tests[i][0]);
    assert_condition(grid[tests[i][1]].get_ngb(AMRNGBPOSITION_RIGHT)->value() ==
                     tests[i][0]);

    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][0]].get_ngb(AMRNGBPOSITION_RIGHT)->value(),
                tests[i][2]);
    assert_condition(grid[tests[i][0]].get_ngb(AMRNGBPOSITION_RIGHT)->value() ==
                     tests[i][2]);
    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][2]].get_ngb(AMRNGBPOSITION_LEFT)->value(),
                tests[i][0]);
    assert_condition(grid[tests[i][2]].get_ngb(AMRNGBPOSITION_LEFT)->value() ==
                     tests[i][0]);

    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][0]].get_ngb(AMRNGBPOSITION_FRONT)->value(),
                tests[i][3]);
    assert_condition(grid[tests[i][0]].get_ngb(AMRNGBPOSITION_FRONT)->value() ==
                     tests[i][3]);
    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][3]].get_ngb(AMRNGBPOSITION_BACK)->value(),
                tests[i][0]);
    assert_condition(grid[tests[i][3]].get_ngb(AMRNGBPOSITION_BACK)->value() ==
                     tests[i][0]);

    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][0]].get_ngb(AMRNGBPOSITION_BACK)->value(),
                tests[i][4]);
    assert_condition(grid[tests[i][0]].get_ngb(AMRNGBPOSITION_BACK)->value() ==
                     tests[i][4]);
    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][4]].get_ngb(AMRNGBPOSITION_FRONT)->value(),
                tests[i][0]);
    assert_condition(grid[tests[i][4]].get_ngb(AMRNGBPOSITION_FRONT)->value() ==
                     tests[i][0]);

    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][0]].get_ngb(AMRNGBPOSITION_BOTTOM)->value(),
                tests[i][5]);
    assert_condition(
        grid[tests[i][0]].get_ngb(AMRNGBPOSITION_BOTTOM)->value() ==
        tests[i][5]);
    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][5]].get_ngb(AMRNGBPOSITION_TOP)->value(),
                tests[i][0]);
    assert_condition(grid[tests[i][5]].get_ngb(AMRNGBPOSITION_TOP)->value() ==
                     tests[i][0]);

    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][0]].get_ngb(AMRNGBPOSITION_TOP)->value(),
                tests[i][6]);
    assert_condition(grid[tests[i][0]].get_ngb(AMRNGBPOSITION_TOP)->value() ==
                     tests[i][6]);
    cmac_status("%g should be %" PRIuFAST64,
                grid[tests[i][6]].get_ngb(AMRNGBPOSITION_BOTTOM)->value(),
                tests[i][0]);
    assert_condition(
        grid[tests[i][6]].get_ngb(AMRNGBPOSITION_BOTTOM)->value() ==
        tests[i][0]);
  }

  // now check one at the border between the two blocks
  // 1023 = 1 111 111 111
  // 767 = 1 011 111 111
  assert_condition(grid[1023].get_ngb(AMRNGBPOSITION_LEFT)->value() == 767.);
  // 0x00100000000002db = 1 011 011 011 + 0x0010000000000000
  assert_condition(grid[1023].get_ngb(AMRNGBPOSITION_RIGHT)->value() ==
                   0x00100000000002db);
  // 895 = 1 101 111 111
  assert_condition(grid[1023].get_ngb(AMRNGBPOSITION_FRONT)->value() == 895.);
  assert_condition(grid[1023].get_ngb(AMRNGBPOSITION_BACK) == nullptr);
  // 959 = 1 110 111 111
  assert_condition(grid[1023].get_ngb(AMRNGBPOSITION_BOTTOM)->value() == 959.);
  assert_condition(grid[1023].get_ngb(AMRNGBPOSITION_TOP) == nullptr);

  // check a random cell
  Box<> testbox = grid[696].get_geometry();
  // 696 = 1 010 111 000
  // anchor_x = 0. + 0.25 + 0. = 0.25
  // anchor_y = 0. + 0.25 + 0.125 = 0.375
  // anchor_z = 0. + 0.25 + 0. = 0.25
  assert_condition(testbox.get_anchor().x() == 0.25);
  assert_condition(testbox.get_anchor().y() == 0.375);
  assert_condition(testbox.get_anchor().z() == 0.25);
  assert_condition(testbox.get_sides().x() == 0.125);
  assert_condition(testbox.get_sides().y() == 0.125);
  assert_condition(testbox.get_sides().z() == 0.125);

  return 0;
}
