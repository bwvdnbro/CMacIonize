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
 * @file testDensitySubGridCreator.cpp
 *
 * @brief Unit test for the DensitySubGridCreator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "DensitySubGridCreator.hpp"
#include "HomogeneousDensityFunction.hpp"

#include <fstream>
#include <vector>

/**
 * @brief Unit test for the DensitySubGridCreator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const CoordinateVector<> box_anchor(0., 0., 0.);
  const CoordinateVector<> box_sides(1., 1., 1.);
  const CoordinateVector< int_fast32_t > ncell(16, 16, 16);
  const CoordinateVector< int_fast32_t > nsubgrid(4, 4, 8);

  DensitySubGridCreator grid_creator(Box<>(box_anchor, box_sides), ncell,
                                     nsubgrid);
  HomogeneousDensityFunction density_function;
  grid_creator.initialize(density_function);

  std::ofstream ofile("testDensitySubGridCreator_grid.txt");
  for (auto gridit = grid_creator.begin(); gridit != grid_creator.all_end();
       ++gridit) {
    assert_condition((*gridit).get_number_of_cells() == 32);
    for (auto cellit = (*gridit).begin(); cellit != (*gridit).end(); ++cellit) {
      cellit.set_neutral_fraction(gridit.get_index());
    }
    (*gridit).print_intensities(ofile);
  }

  // we manually check two subgrids:
  //  - subgrid 0 (0, 0, 0)
  DensitySubGrid &grid0 = *grid_creator.get_subgrid(0);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 41);

  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 9);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 33);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 40);

  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_FACE_X_N) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 32);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_FACE_Y_N) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 8);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_FACE_Z_N) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(grid0.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 1);

  //  - subgrid 82 (2, 2, 2)
  DensitySubGrid &grid82 = *grid_creator.get_subgrid(82);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NNN) == 41);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NNP) == 43);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NPN) == 57);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NPP) == 59);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PNN) == 105);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PNP) == 107);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PPN) == 121);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 123);

  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) == 73);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) == 75);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) == 89);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 91);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) == 49);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) == 51);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) == 113);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 115);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) == 42);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) == 58);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) == 106);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 122);

  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_X_N) == 50);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 114);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 74);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 90);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 81);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 83);

  std::vector< uint_fast8_t > copy_levels(
      grid_creator.number_of_original_subgrids(), 0);
  copy_levels[82] = 1;
  copy_levels[83] = 2;
  grid_creator.create_copies(copy_levels);

  assert_condition(grid_creator.number_of_actual_subgrids() == 132);

  // we repeat the neighbour checks for subgrid 82
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NNN) == 41);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NNP) == 43);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NPN) == 57);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_NPP) == 59);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PNN) == 105);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PNP) == 107);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PPN) == 121);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 123);

  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) == 73);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) == 75);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) == 89);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 91);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) == 49);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) == 51);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) == 113);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 115);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) == 42);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) == 58);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) == 106);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 122);

  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_X_N) == 50);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 114);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 74);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 90);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 81);
  assert_condition(grid82.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 83);

  // also check the copy of 82
  DensitySubGrid &grid128 = *grid_creator.get_subgrid(128);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_NNN) == 41);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_NNP) == 43);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_NPN) == 57);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_NPP) == 59);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_PNN) == 105);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_PNP) == 107);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_PPN) == 121);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 123);

  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) == 73);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) == 75);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) == 89);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 91);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) == 49);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) == 51);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) == 113);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 115);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) == 42);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) == 58);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) == 106);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 122);

  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_FACE_X_N) == 50);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 114);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 74);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 90);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 81);
  assert_condition(grid128.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 129);

  // check subgrid 83 (2, 2, 3)
  DensitySubGrid &grid83 = *grid_creator.get_subgrid(83);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_NNN) == 42);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_NNP) == 44);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_NPN) == 58);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_NPP) == 60);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_PNN) == 106);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_PNP) == 108);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_PPN) == 122);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 124);

  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) == 74);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) == 76);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) == 90);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 92);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) == 50);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) == 52);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) == 114);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 116);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) == 43);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) == 59);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) == 107);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 123);

  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_FACE_X_N) == 51);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 115);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 75);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 91);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 82);
  assert_condition(grid83.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 84);

  // check the 3 copies of subgrid 83:
  //  - copy 1
  DensitySubGrid &grid129 = *grid_creator.get_subgrid(129);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_NNN) == 42);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_NNP) == 44);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_NPN) == 58);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_NPP) == 60);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_PNN) == 106);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_PNP) == 108);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_PPN) == 122);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 124);

  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) == 74);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) == 76);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) == 90);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 92);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) == 50);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) == 52);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) == 114);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 116);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) == 43);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) == 59);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) == 107);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 123);

  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_FACE_X_N) == 51);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 115);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 75);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 91);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 82);
  assert_condition(grid129.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 84);

  //  - copy 2
  DensitySubGrid &grid130 = *grid_creator.get_subgrid(130);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_NNN) == 42);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_NNP) == 44);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_NPN) == 58);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_NPP) == 60);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_PNN) == 106);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_PNP) == 108);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_PPN) == 122);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 124);

  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) == 74);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) == 76);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) == 90);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 92);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) == 50);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) == 52);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) == 114);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 116);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) == 43);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) == 59);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) == 107);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 123);

  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_FACE_X_N) == 51);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 115);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 75);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 91);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 128);
  assert_condition(grid130.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 84);

  //  - copy 3
  DensitySubGrid &grid131 = *grid_creator.get_subgrid(131);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_NNN) == 42);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_NNP) == 44);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_NPN) == 58);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_NPP) == 60);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_PNN) == 106);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_PNP) == 108);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_PPN) == 122);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_CORNER_PPP) == 124);

  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_X_NN) == 74);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_X_NP) == 76);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_X_PN) == 90);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 92);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) == 50);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) == 52);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) == 114);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 116);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) == 43);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) == 59);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) == 107);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 123);

  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_FACE_X_N) == 51);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_FACE_X_P) == 115);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 75);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 91);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 128);
  assert_condition(grid131.get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 84);

  return 0;
}
