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

  const double box_anchor[3] = {0., 0., 0.};
  const double box_sides[3] = {1., 1., 1.};
  const int_fast32_t ncell[3] = {16, 16, 16};
  const int_fast32_t nsubgrid[3] = {4, 4, 8};

  DensitySubGridCreator grid_creator(box_anchor, box_sides, ncell, nsubgrid);
  const int_fast32_t ngrid = nsubgrid[0] * nsubgrid[1] * nsubgrid[2];
  std::vector< DensitySubGrid * > subgrids(ngrid, nullptr);

  std::ofstream ofile("testDensitySubGridCreator_grid.txt");
  for (int_fast32_t igrid = 0; igrid < ngrid; ++igrid) {
    DensitySubGrid *subgrid = grid_creator.create_subgrid(igrid);
    const size_t nsubcells = subgrid->get_number_of_cells();
    assert_condition(nsubcells == 32);
    for (size_t i = 0; i < nsubcells; ++i) {
      subgrid->set_neutral_fraction(i, igrid);
    }
    subgrid->print_intensities(ofile);
    subgrids[igrid] = subgrid;
  }

  // we manually check two subgrids:
  //  - subgrid 0 (0, 0, 0)
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   41);

  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) == 9);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) == 33);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) == 40);

  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_FACE_X_N) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_FACE_X_P) == 32);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 8);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) ==
                   NEIGHBOUR_OUTSIDE);
  assert_condition(subgrids[0]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 1);

  //  - subgrid 82 (2, 2, 2)
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   41);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   43);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   57);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   59);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   105);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   107);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   121);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   123);

  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   73);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   75);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   89);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) ==
                   91);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   49);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   51);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   113);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) ==
                   115);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   42);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   58);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   106);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) ==
                   122);

  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_X_N) == 50);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
                   114);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 74);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 90);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 81);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 83);

  std::vector< uint_fast8_t > copy_levels(ngrid, 0);
  copy_levels[82] = 1;
  copy_levels[83] = 2;
  std::vector< uint_fast32_t > originals, copies;
  grid_creator.create_copies(subgrids, copy_levels, originals, copies);

  assert_condition(subgrids.size() == 132);
  assert_condition(originals.size() == 4);
  assert_condition(copies.size() == ngrid);

  assert_condition(originals[0] == 82);
  assert_condition(copies[82] == ngrid);
  assert_condition(originals[1] == 83);
  assert_condition(copies[83] == ngrid + 1);
  assert_condition(originals[2] == 83);
  assert_condition(originals[3] == 83);

  // we repeat the neighbour checks for subgrid 82
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   41);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   43);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   57);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   59);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   105);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   107);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   121);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   123);

  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   73);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   75);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   89);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) ==
                   91);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   49);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   51);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   113);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) ==
                   115);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   42);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   58);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   106);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) ==
                   122);

  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_X_N) == 50);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
                   114);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 74);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 90);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 81);
  assert_condition(subgrids[82]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 83);

  // also check the copy of 82
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   41);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   43);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   57);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   59);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   105);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   107);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   121);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   123);

  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   73);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   75);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   89);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) ==
                   91);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   49);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   51);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   113);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) ==
                   115);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   42);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   58);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   106);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) ==
                   122);

  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_FACE_X_N) ==
                   50);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
                   114);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) ==
                   74);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) ==
                   90);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) ==
                   81);
  assert_condition(subgrids[128]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) ==
                   129);

  // check subgrid 83 (2, 2, 3)
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   42);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   44);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   58);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   60);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   106);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   108);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   122);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   124);

  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   74);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   76);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   90);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) ==
                   92);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   50);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   52);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   114);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) ==
                   116);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   43);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   59);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   107);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) ==
                   123);

  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_FACE_X_N) == 51);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
                   115);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) == 75);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) == 91);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) == 82);
  assert_condition(subgrids[83]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) == 84);

  // check the 3 copies of subgrid 83:
  //  - copy 1
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   42);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   44);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   58);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   60);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   106);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   108);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   122);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   124);

  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   74);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   76);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   90);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) ==
                   92);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   50);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   52);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   114);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) ==
                   116);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   43);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   59);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   107);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) ==
                   123);

  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_FACE_X_N) ==
                   51);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
                   115);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) ==
                   75);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) ==
                   91);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) ==
                   82);
  assert_condition(subgrids[129]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) ==
                   84);

  //  - copy 2
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   42);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   44);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   58);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   60);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   106);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   108);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   122);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   124);

  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   74);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   76);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   90);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) ==
                   92);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   50);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   52);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   114);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) ==
                   116);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   43);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   59);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   107);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) ==
                   123);

  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_FACE_X_N) ==
                   51);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
                   115);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) ==
                   75);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) ==
                   91);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) ==
                   128);
  assert_condition(subgrids[130]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) ==
                   84);

  //  - copy 3
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_NNN) ==
                   42);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_NNP) ==
                   44);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_NPN) ==
                   58);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_NPP) ==
                   60);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_PNN) ==
                   106);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_PNP) ==
                   108);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_PPN) ==
                   122);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_CORNER_PPP) ==
                   124);

  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_X_NN) ==
                   74);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_X_NP) ==
                   76);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_X_PN) ==
                   90);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_X_PP) ==
                   92);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NN) ==
                   50);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Y_NP) ==
                   52);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PN) ==
                   114);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Y_PP) ==
                   116);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NN) ==
                   43);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Z_NP) ==
                   59);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PN) ==
                   107);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_EDGE_Z_PP) ==
                   123);

  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_FACE_X_N) ==
                   51);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_FACE_X_P) ==
                   115);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_FACE_Y_N) ==
                   75);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_FACE_Y_P) ==
                   91);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_FACE_Z_N) ==
                   128);
  assert_condition(subgrids[131]->get_neighbour(TRAVELDIRECTION_FACE_Z_P) ==
                   84);

  for (size_t igrid = 0; igrid < subgrids.size(); ++igrid) {
    delete subgrids[igrid];
  }

  return 0;
}
