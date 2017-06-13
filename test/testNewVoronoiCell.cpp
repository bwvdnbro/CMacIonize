/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testNewVoronoiCell.cpp
 *
 * @brief Unit test for the NewVoronoiCell class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "NewVoronoiCell.hpp"

#include <fstream>

/**
 * @brief Names for the tests that require specific cell setups.
 */
enum NewVoronoiCellTestName {
  /*! @brief Test for the 1 to 4 flip routine. */
  NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP,
  /*! @brief Test for the 2 to 3 flip routine. */
  NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP,
  /*! @brief Test for the 3 to 2 flip routine. */
  NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP,
  /*! @brief Test for the 4 to 4 flip routine. */
  NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP
};

/**
 * @brief Set up the internal variables for a specific unit test.
 *
 * @param test NewVoronoiCellTestName.
 */
void NewVoronoiCell::setup_test(int test) {
  switch (test) {
  case NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _ngbs.resize(5);

    // convenient names used in documentation figure
    const unsigned int v[5] = {0, 1, 2, 3, 4};
    const unsigned int t[1] = {0};
    const unsigned int ngb[4] = {1, 2, 3, 4};

    // we need 1 tetrahedron + 4 dummy neighbours to check the neighbour
    // handling
    _tetrahedra.resize(5);
    _tetrahedra[t[0]] = VoronoiTetrahedron(v[0], v[1], v[2], v[3], ngb[0],
                                           ngb[1], ngb[2], ngb[3], 0, 0, 0, 0);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _ngbs.resize(5);

    // convenient names used in documentation figure
    const unsigned int v[5] = {0, 1, 2, 3, 4};
    const unsigned int t[2] = {0, 1};
    const unsigned int ngb[6] = {2, 3, 4, 5, 6, 7};

    // we need 2 tetrahedra + 6 dummy neighbours to check the neighbour handling
    _tetrahedra.resize(8);
    _tetrahedra[t[0]] = VoronoiTetrahedron(v[0], v[1], v[2], v[3], ngb[0],
                                           ngb[3], t[1], ngb[4], 0, 0, 3, 0);
    _tetrahedra[t[1]] = VoronoiTetrahedron(v[0], v[1], v[3], v[4], ngb[1],
                                           ngb[2], ngb[5], t[0], 0, 0, 0, 2);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _ngbs.resize(5);

    // convenient names used in documentation figure
    const unsigned int v[5] = {0, 1, 2, 3, 4};
    const unsigned int t[3] = {0, 1, 2};
    const unsigned int ngb[6] = {3, 4, 5, 6, 7, 8};

    // we need 3 tetrahedra + 6 dummy neighbours to check the neighbour handling
    _tetrahedra.resize(9);
    _tetrahedra[t[0]] = VoronoiTetrahedron(v[0], v[1], v[2], v[4], t[2], t[1],
                                           ngb[5], ngb[4], 3, 3, 0, 0);
    _tetrahedra[t[1]] = VoronoiTetrahedron(v[0], v[4], v[2], v[3], t[2], ngb[3],
                                           ngb[2], t[0], 1, 0, 0, 1);
    _tetrahedra[t[2]] = VoronoiTetrahedron(v[4], v[1], v[2], v[3], ngb[0], t[1],
                                           ngb[1], t[0], 0, 0, 0, 0);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _ngbs.resize(6);

    // convenient names used in documentation figure
    const unsigned int v[6] = {0, 1, 2, 3, 4, 5};
    const unsigned int t[4] = {0, 1, 2, 3};
    const unsigned int ngb[8] = {4, 5, 6, 7, 8, 9, 10, 11};

    // we need 4 tetrahedra + 8 dummy neighbours to check the neighbour handling
    _tetrahedra.resize(12);
    _tetrahedra[t[0]] = VoronoiTetrahedron(v[0], v[1], v[2], v[3], ngb[0],
                                           ngb[3], t[1], t[2], 0, 0, 3, 2);
    _tetrahedra[t[1]] = VoronoiTetrahedron(v[0], v[1], v[3], v[4], ngb[1],
                                           ngb[2], t[3], t[0], 0, 0, 1, 2);
    _tetrahedra[t[2]] = VoronoiTetrahedron(v[0], v[1], v[5], v[2], ngb[4],
                                           ngb[7], t[0], t[3], 0, 0, 3, 3);
    _tetrahedra[t[3]] = VoronoiTetrahedron(v[0], v[5], v[1], v[4], ngb[5], t[1],
                                           ngb[6], t[2], 0, 2, 0, 3);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  }
}

/**
 * @brief Check the internal variables after a specific unit test.
 *
 * @param test NewVoronoiCellTestName.
 */
void NewVoronoiCell::check_test(int test) {
  switch (test) {
  case NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP: {

    // three new tetrahedra should have been created by this flip
    assert_condition(_tetrahedra.size() == 8);

    // convenient names used in documentation figure
    const unsigned int v[5] = {0, 1, 2, 3, 4};
    const unsigned int tn[4] = {0, 5, 6, 7};
    const unsigned int ngb[4] = {1, 2, 3, 4};

    // check the individual tetrahedra
    assert_condition(_tetrahedra[tn[0]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(1) == v[1]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(2) == v[2]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(3) == v[4]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(0) == tn[3]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(1) == tn[2]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(2) == tn[1]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(3) == ngb[3]);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(0) == 3);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(1) == 3);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(2) == 3);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(3) == 0);
    assert_condition(_tetrahedra[ngb[3]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[3]].get_ngb_index(0) == 3);

    assert_condition(_tetrahedra[tn[1]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(1) == v[1]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(2) == v[4]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(3) == v[3]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(0) == tn[3]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(1) == tn[2]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(2) == ngb[2]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(3) == tn[0]);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(0) == 2);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(1) == 2);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(3) == 2);
    assert_condition(_tetrahedra[ngb[2]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[2]].get_ngb_index(0) == 2);

    assert_condition(_tetrahedra[tn[2]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(1) == v[4]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(2) == v[2]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(3) == v[3]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(0) == tn[3]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(1) == ngb[1]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(2) == tn[1]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(3) == tn[0]);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(2) == 1);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(3) == 1);
    assert_condition(_tetrahedra[ngb[1]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[ngb[1]].get_ngb_index(0) == 1);

    assert_condition(_tetrahedra[tn[3]].get_vertex(0) == v[4]);
    assert_condition(_tetrahedra[tn[3]].get_vertex(1) == v[1]);
    assert_condition(_tetrahedra[tn[3]].get_vertex(2) == v[2]);
    assert_condition(_tetrahedra[tn[3]].get_vertex(3) == v[3]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(0) == ngb[0]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(1) == tn[2]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(2) == tn[1]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(3) == tn[0]);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(3) == 0);
    assert_condition(_tetrahedra[ngb[0]].get_neighbour(0) == tn[3]);
    assert_condition(_tetrahedra[ngb[0]].get_ngb_index(0) == 0);

    return;
  }
  case NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP: {

    // one extra tetrahedron should have been created by this flip
    assert_condition(_tetrahedra.size() == 9);

    // convenient names used in documentation figure
    const unsigned int v[5] = {0, 1, 2, 3, 4};
    const unsigned int tn[3] = {0, 1, 8};
    const unsigned int ngb[6] = {2, 3, 4, 5, 6, 7};

    // check the individual tetrahedra
    assert_condition(_tetrahedra[tn[0]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(1) == v[1]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(2) == v[2]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(3) == v[4]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(1) == tn[1]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(2) == ngb[5]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(3) == ngb[4]);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(0) == 3);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(1) == 3);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(3) == 0);
    assert_condition(_tetrahedra[ngb[5]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[5]].get_ngb_index(0) == 2);
    assert_condition(_tetrahedra[ngb[4]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[4]].get_ngb_index(0) == 3);

    assert_condition(_tetrahedra[tn[1]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(1) == v[4]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(2) == v[2]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(3) == v[3]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(1) == ngb[3]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(2) == ngb[2]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(3) == tn[0]);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(3) == 1);
    assert_condition(_tetrahedra[ngb[3]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[3]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[ngb[2]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[2]].get_ngb_index(0) == 2);

    assert_condition(_tetrahedra[tn[2]].get_vertex(0) == v[4]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(1) == v[1]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(2) == v[2]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(3) == v[3]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(0) == ngb[0]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(1) == tn[1]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(2) == ngb[1]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(3) == tn[0]);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(3) == 0);
    assert_condition(_tetrahedra[ngb[0]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[ngb[0]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[ngb[1]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[ngb[1]].get_ngb_index(0) == 2);

    return;
  }
  case NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP: {

    // no new tetrahedra should have been created by this flip
    assert_condition(_tetrahedra.size() == 9);
    // one free spot should have opened up
    assert_condition(_free_tetrahedra.size() == 1);
    assert_condition(_free_tetrahedra[0] == 2);

    // convenient names used in documentation figure
    const unsigned int v[5] = {0, 1, 2, 3, 4};
    const unsigned int tn[2] = {0, 1};
    const unsigned int ngb[6] = {3, 4, 5, 6, 7, 8};

    // check the individual tetrahedra
    assert_condition(_tetrahedra[tn[0]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(1) == v[1]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(2) == v[2]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(3) == v[3]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(0) == ngb[0]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(1) == ngb[3]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(2) == tn[1]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(3) == ngb[4]);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(2) == 3);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(3) == 0);
    assert_condition(_tetrahedra[ngb[0]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[0]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[ngb[3]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[3]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[ngb[4]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[4]].get_ngb_index(0) == 3);

    assert_condition(_tetrahedra[tn[1]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(1) == v[1]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(2) == v[3]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(3) == v[4]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(0) == ngb[1]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(1) == ngb[2]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(2) == ngb[5]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(3) == tn[0]);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(3) == 2);
    assert_condition(_tetrahedra[ngb[1]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[1]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[ngb[2]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[2]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[ngb[5]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[5]].get_ngb_index(0) == 2);

    return;
  }
  case NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP: {

    // no new tetrahedra should have been created by this flip
    assert_condition(_tetrahedra.size() == 12);

    // convenient names used in documentation figure
    const unsigned int v[6] = {0, 1, 2, 3, 4, 5};
    const unsigned int tn[4] = {0, 1, 2, 3};
    const unsigned int ngb[8] = {4, 5, 6, 7, 8, 9, 10, 11};

    // check the individual tetrahedra
    assert_condition(_tetrahedra[tn[0]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(1) == v[3]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(2) == v[5]);
    assert_condition(_tetrahedra[tn[0]].get_vertex(3) == v[2]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(1) == ngb[7]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(2) == ngb[3]);
    assert_condition(_tetrahedra[tn[0]].get_neighbour(3) == tn[2]);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[0]].get_ngb_index(3) == 3);
    assert_condition(_tetrahedra[ngb[7]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[7]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[ngb[3]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[ngb[3]].get_ngb_index(0) == 2);

    assert_condition(_tetrahedra[tn[1]].get_vertex(0) == v[1]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(1) == v[5]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(2) == v[3]);
    assert_condition(_tetrahedra[tn[1]].get_vertex(3) == v[2]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(0) == tn[0]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(1) == ngb[0]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(2) == ngb[4]);
    assert_condition(_tetrahedra[tn[1]].get_neighbour(3) == tn[3]);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[1]].get_ngb_index(3) == 3);
    assert_condition(_tetrahedra[ngb[0]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[0]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[ngb[4]].get_neighbour(0) == tn[1]);
    assert_condition(_tetrahedra[ngb[4]].get_ngb_index(0) == 2);

    assert_condition(_tetrahedra[tn[2]].get_vertex(0) == v[0]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(1) == v[5]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(2) == v[3]);
    assert_condition(_tetrahedra[tn[2]].get_vertex(3) == v[4]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(0) == tn[3]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(1) == ngb[2]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(2) == ngb[6]);
    assert_condition(_tetrahedra[tn[2]].get_neighbour(3) == tn[0]);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[2]].get_ngb_index(3) == 3);
    assert_condition(_tetrahedra[ngb[2]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[ngb[2]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[ngb[6]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[ngb[6]].get_ngb_index(0) == 2);

    assert_condition(_tetrahedra[tn[3]].get_vertex(0) == v[1]);
    assert_condition(_tetrahedra[tn[3]].get_vertex(1) == v[3]);
    assert_condition(_tetrahedra[tn[3]].get_vertex(2) == v[5]);
    assert_condition(_tetrahedra[tn[3]].get_vertex(3) == v[4]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(0) == tn[2]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(1) == ngb[5]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(2) == ngb[1]);
    assert_condition(_tetrahedra[tn[3]].get_neighbour(3) == tn[1]);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(0) == 0);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(1) == 0);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(2) == 0);
    assert_condition(_tetrahedra[tn[3]].get_ngb_index(3) == 3);
    assert_condition(_tetrahedra[ngb[5]].get_neighbour(0) == tn[3]);
    assert_condition(_tetrahedra[ngb[5]].get_ngb_index(0) == 1);
    assert_condition(_tetrahedra[ngb[1]].get_neighbour(0) == tn[3]);
    assert_condition(_tetrahedra[ngb[1]].get_ngb_index(0) == 2);

    return;
  }
  }
}

/**
 * @brief Check if NewVoronoiCell::get_positive_permutation can complete the
 * given permutation of 0123.
 *
 * @param a First index.
 * @param b Second index.
 * @param c Third index.
 * @param d Fourth index.
 */
inline void check_permutation(unsigned char a, unsigned char b, unsigned char c,
                              unsigned char d) {
  unsigned char v[4];
  v[0] = a;
  v[1] = b;
  NewVoronoiCell::get_positive_permutation(v);
  assert_condition(v[2] == c);
  assert_condition(v[3] == d);
  assert_condition(
      NewVoronoiCell::positive_permutation(v[0], v[1], v[2], v[3]) == true);
}

/**
 * @brief Unit test for the NewVoronoiCell class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  /// test VoronoiTetrahedron
  {
    std::vector< CoordinateVector<> > positions(4);
    positions[0] = CoordinateVector<>(1., 1., 1.);
    positions[1] = CoordinateVector<>(-1., -1., 1.);
    positions[2] = CoordinateVector<>(-1., 1., -1.);
    positions[3] = CoordinateVector<>(1., -1., -1.);
    VoronoiTetrahedron tetrahedron(0, 1, 2, 3);
    CoordinateVector<> midpoint =
        tetrahedron.get_midpoint_circumsphere(positions);
    assert_condition(midpoint.x() == 0.);
    assert_condition(midpoint.y() == 0.);
    assert_condition(midpoint.z() == 0.);
    cmac_status("Midpoint tetrahedron works!");

    CoordinateVector<> centroid = tetrahedron.get_centroid(positions);
    assert_condition(centroid.x() == 0.);
    assert_condition(centroid.y() == 0.);
    assert_condition(centroid.z() == 0.);
    cmac_status("Centroid tetrahedron works!");

    double volume = tetrahedron.get_volume(positions);
    assert_condition(volume == 8. / 3.);
    cmac_status("Volume tetrahedron works!");
  }

  /// test geometrical routines
  {
    Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    std::vector< CoordinateVector<> > positions(1);
    positions[0] = CoordinateVector<>(0.25, 0.25, 0.25);
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > long_positions(1);
    long_positions[0] = CoordinateVector< unsigned long >(1250);
    VoronoiBox< unsigned long > voronoi_box(long_positions[0], box_anchor,
                                            box_sides);
    NewVoronoiCell cell(0, voronoi_box, long_positions);
    cell.finalize(box, positions);

    std::vector< CoordinateVector<> > full_volume_positions(4);
    full_volume_positions[0] = CoordinateVector<>(6.26786);
    full_volume_positions[1] = CoordinateVector<>(-8.125, 3.5, 3.5);
    full_volume_positions[2] = CoordinateVector<>(3.5, -8.125, 3.5);
    full_volume_positions[3] = CoordinateVector<>(3.5, 3.5, -8.125);
    VoronoiTetrahedron full_volume_tetrahedron(0, 1, 2, 3);
    const double full_volume =
        full_volume_tetrahedron.get_volume(full_volume_positions);
    const CoordinateVector<> full_centroid =
        full_volume_tetrahedron.get_centroid(full_volume_positions);
    const double full_area_023 =
        0.5 *
        CoordinateVector<>::cross_product(
            full_volume_positions[2] - full_volume_positions[0],
            full_volume_positions[3] - full_volume_positions[0])
            .norm();
    const CoordinateVector<> full_midpoint_023 =
        (full_volume_positions[0] + full_volume_positions[2] +
         full_volume_positions[3]) /
        3.;
    const double full_area_123 =
        0.5 *
        CoordinateVector<>::cross_product(
            full_volume_positions[2] - full_volume_positions[1],
            full_volume_positions[3] - full_volume_positions[1])
            .norm();
    const CoordinateVector<> full_midpoint_123 =
        (full_volume_positions[1] + full_volume_positions[2] +
         full_volume_positions[3]) /
        3.;
    const double full_area_013 =
        0.5 *
        CoordinateVector<>::cross_product(
            full_volume_positions[1] - full_volume_positions[0],
            full_volume_positions[3] - full_volume_positions[0])
            .norm();
    const CoordinateVector<> full_midpoint_013 =
        (full_volume_positions[0] + full_volume_positions[1] +
         full_volume_positions[3]) /
        3.;
    const double full_area_012 =
        0.5 *
        CoordinateVector<>::cross_product(
            full_volume_positions[1] - full_volume_positions[0],
            full_volume_positions[2] - full_volume_positions[0])
            .norm();
    const CoordinateVector<> full_midpoint_012 =
        (full_volume_positions[0] + full_volume_positions[1] +
         full_volume_positions[2]) /
        3.;

    const double tolerance = 1.e-6;

    double volume = cell.get_volume();
    assert_values_equal_rel(volume, full_volume, tolerance);
    cmac_status("Cell volume computation works!");

    CoordinateVector<> centroid = cell.get_centroid();
    assert_values_equal_rel(centroid.x(), full_centroid.x(), tolerance);
    assert_values_equal_rel(centroid.y(), full_centroid.y(), tolerance);
    assert_values_equal_rel(centroid.z(), full_centroid.z(), tolerance);
    cmac_status("Cell centroid computation works!");

    std::vector< VoronoiFace > faces = cell.get_faces();
    assert_condition(faces.size() == 4);

    assert_values_equal_rel(faces[0].get_surface_area(), full_area_023,
                            tolerance);
    CoordinateVector<> midpoint = faces[0].get_midpoint();
    assert_values_equal_rel(midpoint.x(), full_midpoint_023.x(), tolerance);
    assert_values_equal_rel(midpoint.y(), full_midpoint_023.y(), tolerance);
    assert_values_equal_rel(midpoint.z(), full_midpoint_023.z(), tolerance);
    assert_condition(faces[0].get_neighbour() == NEWVORONOICELL_BOX_CORNER1);

    assert_values_equal_rel(faces[1].get_surface_area(), full_area_013,
                            tolerance);
    midpoint = faces[1].get_midpoint();
    assert_values_equal_rel(midpoint.x(), full_midpoint_013.x(), tolerance);
    assert_values_equal_rel(midpoint.y(), full_midpoint_013.y(), tolerance);
    assert_values_equal_rel(midpoint.z(), full_midpoint_013.z(), tolerance);
    assert_condition(faces[1].get_neighbour() == NEWVORONOICELL_BOX_CORNER2);

    assert_values_equal_rel(faces[2].get_surface_area(), full_area_012,
                            tolerance);
    midpoint = faces[2].get_midpoint();
    assert_values_equal_rel(midpoint.x(), full_midpoint_012.x(), tolerance);
    assert_values_equal_rel(midpoint.y(), full_midpoint_012.y(), tolerance);
    assert_values_equal_rel(midpoint.z(), full_midpoint_012.z(), tolerance);
    assert_condition(faces[2].get_neighbour() == NEWVORONOICELL_BOX_CORNER3);

    assert_values_equal_rel(faces[3].get_surface_area(), full_area_123,
                            tolerance);
    midpoint = faces[3].get_midpoint();
    assert_values_equal_rel(midpoint.x(), full_midpoint_123.x(), tolerance);
    assert_values_equal_rel(midpoint.y(), full_midpoint_123.y(), tolerance);
    assert_values_equal_rel(midpoint.z(), full_midpoint_123.z(), tolerance);
    assert_condition(faces[3].get_neighbour() == NEWVORONOICELL_BOX_CORNER0);

    cmac_status("Cell face computation works!");
  }

  /// test for NewVoronoiCell::find_tetrahedron
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(4);
    positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    positions[1] = CoordinateVector< unsigned long >(1250, 1500, 1500);
    // point on axis
    positions[2] = CoordinateVector< unsigned long >(1250);
    // point on face
    positions[3] = CoordinateVector< unsigned long >(1250, 1250, 1500);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    NewVoronoiCell cell(0, box, positions);

    unsigned int tetrahedra[4];
    assert_condition(cell.find_tetrahedron(1, box, positions, tetrahedra) == 1);
    assert_condition(tetrahedra[0] == 1);

    assert_condition(cell.find_tetrahedron(2, box, positions, tetrahedra) == 3);
    assert_condition(tetrahedra[0] == 3);
    assert_condition(tetrahedra[1] == 1);
    assert_condition(tetrahedra[2] == 2);

    assert_condition(cell.find_tetrahedron(3, box, positions, tetrahedra) == 2);
    assert_condition(tetrahedra[0] == 1);
    assert_condition(tetrahedra[1] == 2);

    cmac_status("Find tetrahedron works!");
  }

  /// test permutation routine
  {
    assert_condition(NewVoronoiCell::positive_permutation(0, 1, 2, 3) == true);
    assert_condition(NewVoronoiCell::positive_permutation(0, 2, 3, 1) == true);
    assert_condition(NewVoronoiCell::positive_permutation(0, 3, 1, 2) == true);

    assert_condition(NewVoronoiCell::positive_permutation(1, 0, 3, 2) == true);
    assert_condition(NewVoronoiCell::positive_permutation(1, 2, 0, 3) == true);
    assert_condition(NewVoronoiCell::positive_permutation(1, 3, 2, 0) == true);

    assert_condition(NewVoronoiCell::positive_permutation(2, 0, 1, 3) == true);
    assert_condition(NewVoronoiCell::positive_permutation(2, 1, 3, 0) == true);
    assert_condition(NewVoronoiCell::positive_permutation(2, 3, 0, 1) == true);

    assert_condition(NewVoronoiCell::positive_permutation(3, 0, 2, 1) == true);
    assert_condition(NewVoronoiCell::positive_permutation(3, 1, 0, 2) == true);
    assert_condition(NewVoronoiCell::positive_permutation(3, 2, 1, 0) == true);

    assert_condition(NewVoronoiCell::positive_permutation(0, 1, 3, 2) == false);
    assert_condition(NewVoronoiCell::positive_permutation(0, 2, 1, 3) == false);
    assert_condition(NewVoronoiCell::positive_permutation(0, 3, 2, 1) == false);

    assert_condition(NewVoronoiCell::positive_permutation(1, 0, 2, 3) == false);
    assert_condition(NewVoronoiCell::positive_permutation(1, 2, 3, 0) == false);
    assert_condition(NewVoronoiCell::positive_permutation(1, 3, 0, 2) == false);

    assert_condition(NewVoronoiCell::positive_permutation(2, 0, 3, 1) == false);
    assert_condition(NewVoronoiCell::positive_permutation(2, 1, 0, 3) == false);
    assert_condition(NewVoronoiCell::positive_permutation(2, 3, 1, 0) == false);

    assert_condition(NewVoronoiCell::positive_permutation(3, 0, 1, 2) == false);
    assert_condition(NewVoronoiCell::positive_permutation(3, 1, 2, 0) == false);
    assert_condition(NewVoronoiCell::positive_permutation(3, 2, 0, 1) == false);

    cmac_status("Permutation checks work!");
  }

  /// test permutation completion routine
  {
    check_permutation(0, 1, 2, 3);
    check_permutation(0, 2, 3, 1);
    check_permutation(0, 3, 1, 2);

    check_permutation(1, 0, 3, 2);
    check_permutation(1, 2, 0, 3);
    check_permutation(1, 3, 2, 0);

    check_permutation(2, 0, 1, 3);
    check_permutation(2, 1, 3, 0);
    check_permutation(2, 3, 0, 1);

    check_permutation(3, 0, 2, 1);
    check_permutation(3, 1, 0, 2);
    check_permutation(3, 2, 1, 0);

    cmac_status("Permutation completion works!");
  }

  /// test 1 to 4 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(1);
    positions[0] = CoordinateVector< unsigned long >(1500);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    NewVoronoiCell cell(0, box, positions);
    cell.setup_test(NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP);
    std::vector< bool > queue(5, false);
    cell.one_to_four_flip(4, 0, queue);
    assert_condition(queue.size() == 8);
    assert_condition(queue[0] == true);
    assert_condition(queue[1] == false);
    assert_condition(queue[2] == false);
    assert_condition(queue[3] == false);
    assert_condition(queue[4] == false);
    assert_condition(queue[5] == true);
    assert_condition(queue[6] == true);
    assert_condition(queue[7] == true);
    cell.check_test(NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP);

    cmac_status("1 to 4 flip works!");
  }

  /// test 2 to 3 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(1);
    positions[0] = CoordinateVector< unsigned long >(1500);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    NewVoronoiCell cell(0, box, positions);
    cell.setup_test(NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP);
    std::vector< bool > queue(8, false);
    unsigned int next_check = cell.two_to_three_flip(0, 1, 2, 3, queue, 8);
    assert_condition(next_check == 0);
    assert_condition(queue.size() == 9);
    assert_condition(queue[0] == true);
    assert_condition(queue[1] == true);
    assert_condition(queue[2] == false);
    assert_condition(queue[3] == false);
    assert_condition(queue[4] == false);
    assert_condition(queue[5] == false);
    assert_condition(queue[6] == false);
    assert_condition(queue[7] == false);
    assert_condition(queue[8] == true);
    cell.check_test(NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP);

    cmac_status("2 to 3 flip works!");
  }

  /// test 3 to 2 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(1);
    positions[0] = CoordinateVector< unsigned long >(1500);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    NewVoronoiCell cell(0, box, positions);
    cell.setup_test(NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP);
    std::vector< bool > queue(9, false);
    unsigned int next_check = cell.three_to_two_flip(0, 1, 2, queue, 8);
    assert_condition(next_check == 0);
    assert_condition(queue.size() == 9);
    assert_condition(queue[0] == true);
    assert_condition(queue[1] == true);
    assert_condition(queue[3] == false);
    assert_condition(queue[4] == false);
    assert_condition(queue[5] == false);
    assert_condition(queue[6] == false);
    assert_condition(queue[7] == false);
    assert_condition(queue[8] == false);
    cell.check_test(NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP);

    cmac_status("3 to 2 flip works!");
  }

  /// test 4 to 4 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(1);
    positions[0] = CoordinateVector< unsigned long >(1500);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    NewVoronoiCell cell(0, box, positions);
    cell.setup_test(NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP);
    std::vector< bool > queue(12, false);
    unsigned int next_check = cell.four_to_four_flip(0, 1, 2, 3, queue, 11);
    assert_condition(next_check == 0);
    assert_condition(queue.size() == 12);
    assert_condition(queue[0] == true);
    assert_condition(queue[1] == true);
    assert_condition(queue[2] == true);
    assert_condition(queue[3] == true);
    assert_condition(queue[4] == false);
    assert_condition(queue[5] == false);
    assert_condition(queue[6] == false);
    assert_condition(queue[7] == false);
    assert_condition(queue[8] == false);
    assert_condition(queue[9] == false);
    assert_condition(queue[10] == false);
    assert_condition(queue[11] == false);
    cell.check_test(NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP);

    cmac_status("4 to 4 flip works!");
  }

  /// tests for NewVoronoiCell::intersect
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(2);
    positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    positions[1] = CoordinateVector< unsigned long >(1250, 1500, 1500);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    NewVoronoiCell cell(0, box, positions);

    cell.intersect(1, box, positions);

    cell.check_empty_circumsphere(box, positions);

    std::ofstream ofile("new_voronoi_cell.txt");
    cell.print_tetrahedra(ofile, box, positions);

    cmac_status("First intersection worked!");
  }
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(2);
    positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    positions[1] = CoordinateVector< unsigned long >(1005, 1500, 1250);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    NewVoronoiCell cell(0, box, positions);

    cell.intersect(1, box, positions);

    cell.check_empty_circumsphere(box, positions);

    std::ofstream ofile("new_voronoi_cell.txt");
    cell.print_tetrahedra(ofile, box, positions);
    ofile << positions[1].x() << "\t" << positions[1].y() << "\t"
          << positions[1].z() << "\n";
    ofile.close();

    cmac_status("Second intersection worked!");
  }

  return 0;
}
