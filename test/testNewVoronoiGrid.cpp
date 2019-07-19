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
 * @file testNewVoronoiGrid.cpp
 *
 * @brief Unit test for the NewVoronoiGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "NewVoronoiCellConstructor.hpp"
#include "NewVoronoiGrid.hpp"
#include "Timer.hpp"
#include "Utilities.hpp"

#include <fstream>

/**
 * @brief Names for the tests that require specific cell setups.
 */
enum NewVoronoiCellTestName {
  /*! @brief Test for the 1 to 4 flip routine. */
  NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP,
  /*! @brief Test for the degenerate 2 to 6 flip routine. */
  NEWVORONOICELL_TEST_TWO_TO_SIX_FLIP,
  /*! @brief Test for the degenerate n to 2n flip routine. */
  NEWVORONOICELL_TEST_N_TO_2N_FLIP,
  /*! @brief Test for the 2 to 3 flip routine. */
  NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP,
  /*! @brief Test for the 3 to 2 flip routine. */
  NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP,
  /*! @brief Test for the degenerate 4 to 4 flip routine. */
  NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP
};

/**
 * @brief Set up the internal variables for a specific unit test.
 *
 * @param test NewVoronoiCellTestName.
 */
void NewVoronoiCellConstructor::setup_test(int_fast32_t test) {
  switch (test) {
  case NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _vertices_size = 5;

    // convenient names used in documentation figure
    const uint_fast32_t v[5] = {0, 1, 2, 3, 4};
    const uint_fast32_t t[1] = {0};
    const uint_fast32_t ngb[4] = {1, 2, 3, 4};

    // we need 1 tetrahedron + 4 dummy neighbours to check the neighbour
    // handling
    _tetrahedra_size = 5;
    _tetrahedra[t[0]] = NewVoronoiTetrahedron(
        v[0], v[1], v[2], v[3], ngb[0], ngb[1], ngb[2], ngb[3], 0, 0, 0, 0);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_TWO_TO_SIX_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _vertices_size = 6;

    // convenient names used in documentation figure
    const uint_fast32_t v[6] = {0, 1, 2, 3, 4, 5};
    const uint_fast32_t t[2] = {0, 1};
    const uint_fast32_t ngb[6] = {2, 3, 4, 5, 6, 7};

    // we need 2 tetrahedra + 6 dummy neighbours to check the neighbour handling
    _tetrahedra_size = 8;
    _tetrahedra[t[0]] = NewVoronoiTetrahedron(v[0], v[1], v[2], v[3], ngb[0],
                                              ngb[3], t[1], ngb[4], 0, 0, 3, 0);
    _tetrahedra[t[1]] = NewVoronoiTetrahedron(v[0], v[1], v[3], v[4], ngb[1],
                                              ngb[2], ngb[5], t[0], 0, 0, 0, 2);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_N_TO_2N_FLIP: {
    // we test for n = 5

    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _vertices_size = 8;

    // convenient names used in documentation figure
    const uint_fast32_t v[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    const uint_fast32_t t[5] = {0, 1, 2, 3, 4};
    const uint_fast32_t ngb[10] = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

    // we need 5 tetrahedra + 10 dummy neighbours to check the neighbour
    // handling
    _tetrahedra_size = 15;
    _tetrahedra[t[0]] = NewVoronoiTetrahedron(v[0], v[6], v[1], v[5], t[1],
                                              ngb[0], t[4], ngb[1], 2, 0, 0, 0);
    _tetrahedra[t[1]] = NewVoronoiTetrahedron(v[1], v[6], v[2], v[5], t[2],
                                              ngb[2], t[0], ngb[3], 2, 0, 0, 0);
    _tetrahedra[t[2]] = NewVoronoiTetrahedron(v[2], v[6], v[3], v[5], t[3],
                                              ngb[4], t[1], ngb[5], 2, 0, 0, 0);
    _tetrahedra[t[3]] = NewVoronoiTetrahedron(v[3], v[6], v[4], v[5], t[4],
                                              ngb[6], t[2], ngb[7], 2, 0, 0, 0);
    _tetrahedra[t[4]] = NewVoronoiTetrahedron(v[4], v[6], v[0], v[5], t[0],
                                              ngb[8], t[3], ngb[9], 2, 0, 0, 0);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _vertices_size = 5;

    // convenient names used in documentation figure
    const uint_fast32_t v[5] = {0, 1, 2, 3, 4};
    const uint_fast32_t t[2] = {0, 1};
    const uint_fast32_t ngb[6] = {2, 3, 4, 5, 6, 7};

    // we need 2 tetrahedra + 6 dummy neighbours to check the neighbour handling
    _tetrahedra_size = 8;
    _tetrahedra[t[0]] = NewVoronoiTetrahedron(v[0], v[1], v[2], v[3], ngb[0],
                                              ngb[3], t[1], ngb[4], 0, 0, 3, 0);
    _tetrahedra[t[1]] = NewVoronoiTetrahedron(v[0], v[1], v[3], v[4], ngb[1],
                                              ngb[2], ngb[5], t[0], 0, 0, 0, 2);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _vertices_size = 5;

    // convenient names used in documentation figure
    const uint_fast32_t v[5] = {0, 1, 2, 3, 4};
    const uint_fast32_t t[3] = {0, 1, 2};
    const uint_fast32_t ngb[6] = {3, 4, 5, 6, 7, 8};

    // we need 3 tetrahedra + 6 dummy neighbours to check the neighbour handling
    _tetrahedra_size = 9;
    _tetrahedra[t[0]] = NewVoronoiTetrahedron(v[0], v[1], v[2], v[4], t[2],
                                              t[1], ngb[5], ngb[4], 3, 3, 0, 0);
    _tetrahedra[t[1]] = NewVoronoiTetrahedron(v[0], v[4], v[2], v[3], t[2],
                                              ngb[3], ngb[2], t[0], 1, 0, 0, 1);
    _tetrahedra[t[2]] = NewVoronoiTetrahedron(v[4], v[1], v[2], v[3], ngb[0],
                                              t[1], ngb[1], t[0], 0, 0, 0, 0);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  case NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP: {
    // the exact positions of the vertices are irrelevant for this test, so we
    // just use the initial values (0)
    _vertices_size = 6;

    // convenient names used in documentation figure
    const uint_fast32_t v[6] = {0, 1, 2, 3, 4, 5};
    const uint_fast32_t t[4] = {0, 1, 2, 3};
    const uint_fast32_t ngb[8] = {4, 5, 6, 7, 8, 9, 10, 11};

    // we need 4 tetrahedra + 8 dummy neighbours to check the neighbour handling
    _tetrahedra_size = 12;
    _tetrahedra[t[0]] = NewVoronoiTetrahedron(v[0], v[1], v[2], v[3], ngb[0],
                                              ngb[3], t[1], t[2], 0, 0, 3, 2);
    _tetrahedra[t[1]] = NewVoronoiTetrahedron(v[0], v[1], v[3], v[4], ngb[1],
                                              ngb[2], t[3], t[0], 0, 0, 1, 2);
    _tetrahedra[t[2]] = NewVoronoiTetrahedron(v[0], v[1], v[5], v[2], ngb[4],
                                              ngb[7], t[0], t[3], 0, 0, 3, 3);
    _tetrahedra[t[3]] = NewVoronoiTetrahedron(v[0], v[5], v[1], v[4], ngb[5],
                                              t[1], ngb[6], t[2], 0, 2, 0, 3);
    // we do not initialize the ngb tetrahedra, as their exact setup is
    // irrelevant for this test

    return;
  }
  }
}

/**
 * @brief Check that the elements of the given tetrahedron are what they should
 * be.
 *
 * @param tetrahedron Index of the tetrahedron to check.
 * @param v0 Expected first vertex.
 * @param v1 Expected second vertex.
 * @param v2 Expected third vertex.
 * @param v3 Expected fourth vertex.
 * @param ngb0 Expected first neighbour.
 * @param ngb1 Expected second neighbour.
 * @param ngb2 Expected third neighbour.
 * @param ngb3 Expected fourth neighbour.
 * @param ngbi0 Expected first neighbour index.
 * @param ngbi1 Expected second neighbour index.
 * @param ngbi2 Expected third neighbour index.
 * @param ngbi3 Expected fourth neighbour index.
 */
#define assert_tetrahedron_right(tetrahedron, v0, v1, v2, v3, ngb0, ngb1,      \
                                 ngb2, ngb3, ngbi0, ngbi1, ngbi2, ngbi3)       \
  assert_condition(_tetrahedra[tetrahedron].get_vertex(0) == v0);              \
  assert_condition(_tetrahedra[tetrahedron].get_vertex(1) == v1);              \
  assert_condition(_tetrahedra[tetrahedron].get_vertex(2) == v2);              \
  assert_condition(_tetrahedra[tetrahedron].get_vertex(3) == v3);              \
  assert_condition(_tetrahedra[tetrahedron].get_neighbour(0) == ngb0);         \
  assert_condition(_tetrahedra[tetrahedron].get_neighbour(1) == ngb1);         \
  assert_condition(_tetrahedra[tetrahedron].get_neighbour(2) == ngb2);         \
  assert_condition(_tetrahedra[tetrahedron].get_neighbour(3) == ngb3);         \
  assert_condition(_tetrahedra[tetrahedron].get_ngb_index(0) == ngbi0);        \
  assert_condition(_tetrahedra[tetrahedron].get_ngb_index(1) == ngbi1);        \
  assert_condition(_tetrahedra[tetrahedron].get_ngb_index(2) == ngbi2);        \
  assert_condition(_tetrahedra[tetrahedron].get_ngb_index(3) == ngbi3);

/**
 * @brief Check that the first neighbour of the given tetrahedron matches the
 * given neighbour and neighbour index.
 *
 * @param tetrahedron Index of the tetrahedron to check.
 * @param neighbour Expected first neighbour.
 * @param ngb_index Expected first neighbour index.
 */
#define assert_neighbour_right(tetrahedron, neighbour, ngb_index)              \
  assert_condition(_tetrahedra[tetrahedron].get_neighbour(0) == neighbour);    \
  assert_condition(_tetrahedra[tetrahedron].get_ngb_index(0) == ngb_index);

/**
 * @brief Check the internal variables after a specific unit test.
 *
 * @param test NewVoronoiCellTestName.
 */
void NewVoronoiCellConstructor::check_test(int_fast32_t test) {
  switch (test) {
  case NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP: {

    // three new tetrahedra should have been created by this flip
    assert_condition(_tetrahedra_size == 8);

    // convenient names used in documentation figure
    const uint_fast32_t v[5] = {0, 1, 2, 3, 4};
    const uint_fast32_t tn[4] = {0, 5, 6, 7};
    const uint_fast32_t ngb[4] = {1, 2, 3, 4};

    // check the individual tetrahedra
    assert_tetrahedron_right(tn[0], v[0], v[1], v[2], v[4], tn[3], tn[2], tn[1],
                             ngb[3], 3, 3, 3, 0);
    assert_neighbour_right(ngb[3], tn[0], 3);

    assert_tetrahedron_right(tn[1], v[0], v[1], v[4], v[3], tn[3], tn[2],
                             ngb[2], tn[0], 2, 2, 0, 2);
    assert_neighbour_right(ngb[2], tn[1], 2);

    assert_tetrahedron_right(tn[2], v[0], v[4], v[2], v[3], tn[3], ngb[1],
                             tn[1], tn[0], 1, 0, 1, 1);
    assert_neighbour_right(ngb[1], tn[2], 1);

    assert_tetrahedron_right(tn[3], v[4], v[1], v[2], v[3], ngb[0], tn[2],
                             tn[1], tn[0], 0, 0, 0, 0);
    assert_neighbour_right(ngb[0], tn[3], 0);

    return;
  }
  case NEWVORONOICELL_TEST_TWO_TO_SIX_FLIP: {

    // four extra tetrahedra should have been created by this flip
    assert_condition(_tetrahedra_size == 12);

    // convenient names used in documentation figure
    const uint_fast32_t v[6] = {0, 1, 2, 3, 4, 5};
    const uint_fast32_t tn[6] = {0, 1, 8, 9, 10, 11};
    const uint_fast32_t ngb[6] = {2, 3, 4, 5, 6, 7};

    // check the individual tetrahedra
    assert_tetrahedron_right(tn[0], v[0], v[1], v[2], v[5], tn[2], tn[1], tn[3],
                             ngb[4], 3, 3, 3, 0);
    assert_neighbour_right(ngb[4], tn[0], 3);

    assert_tetrahedron_right(tn[1], v[0], v[5], v[2], v[3], tn[2], ngb[3],
                             tn[4], tn[0], 1, 0, 3, 1);
    assert_neighbour_right(ngb[3], tn[1], 1);

    assert_tetrahedron_right(tn[2], v[5], v[1], v[2], v[3], ngb[0], tn[1],
                             tn[5], tn[0], 0, 0, 3, 0);
    assert_neighbour_right(ngb[0], tn[2], 0);

    assert_tetrahedron_right(tn[3], v[0], v[1], v[5], v[4], tn[5], tn[4],
                             ngb[5], tn[0], 2, 2, 0, 2);
    assert_neighbour_right(ngb[5], tn[3], 2);

    assert_tetrahedron_right(tn[4], v[0], v[5], v[3], v[4], tn[5], ngb[2],
                             tn[3], tn[1], 1, 0, 1, 2);
    assert_neighbour_right(ngb[2], tn[4], 1);

    assert_tetrahedron_right(tn[5], v[5], v[1], v[3], v[4], ngb[1], tn[4],
                             tn[3], tn[2], 0, 0, 0, 2);
    assert_neighbour_right(ngb[1], tn[5], 0);

    return;
  }
  case NEWVORONOICELL_TEST_N_TO_2N_FLIP: {

    // we test for n = 5

    // n extra tetrahedra should have been created by this flip
    assert_condition(_tetrahedra_size == 20);

    // convenient names used in documentation figure
    const uint_fast32_t v[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    const uint_fast32_t tn[10] = {0, 1, 2, 3, 4, 15, 16, 17, 18, 19};
    const uint_fast32_t ngb[10] = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

    // check the individual tetrahedra
    assert_tetrahedron_right(tn[0], v[0], v[7], v[1], v[5], tn[2], ngb[0],
                             tn[8], tn[1], 2, 0, 0, 1);
    assert_neighbour_right(ngb[0], tn[0], 1);

    assert_tetrahedron_right(tn[1], v[0], v[6], v[1], v[7], tn[3], tn[0], tn[9],
                             ngb[1], 2, 3, 0, 0);
    assert_neighbour_right(ngb[1], tn[1], 3);

    assert_tetrahedron_right(tn[2], v[1], v[7], v[2], v[5], tn[4], ngb[2],
                             tn[0], tn[3], 2, 0, 0, 1);
    assert_neighbour_right(ngb[2], tn[2], 1);

    assert_tetrahedron_right(tn[3], v[1], v[6], v[2], v[7], tn[5], tn[2], tn[1],
                             ngb[3], 2, 3, 0, 0);
    assert_neighbour_right(ngb[3], tn[3], 3);

    assert_tetrahedron_right(tn[4], v[2], v[7], v[3], v[5], tn[6], ngb[4],
                             tn[2], tn[5], 2, 0, 0, 1);
    assert_neighbour_right(ngb[4], tn[4], 1);

    assert_tetrahedron_right(tn[5], v[2], v[6], v[3], v[7], tn[7], tn[4], tn[3],
                             ngb[5], 2, 3, 0, 0);
    assert_neighbour_right(ngb[5], tn[5], 3);

    assert_tetrahedron_right(tn[6], v[3], v[7], v[4], v[5], tn[8], ngb[6],
                             tn[4], tn[7], 2, 0, 0, 1);
    assert_neighbour_right(ngb[6], tn[6], 1);

    assert_tetrahedron_right(tn[7], v[3], v[6], v[4], v[7], tn[9], tn[6], tn[5],
                             ngb[7], 2, 3, 0, 0);
    assert_neighbour_right(ngb[7], tn[7], 3);

    assert_tetrahedron_right(tn[8], v[4], v[7], v[0], v[5], tn[0], ngb[8],
                             tn[6], tn[9], 2, 0, 0, 1);
    assert_neighbour_right(ngb[8], tn[8], 1);

    assert_tetrahedron_right(tn[9], v[4], v[6], v[0], v[7], tn[1], tn[8], tn[7],
                             ngb[9], 2, 3, 0, 0);
    assert_neighbour_right(ngb[9], tn[9], 3);

    return;
  }
  case NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP: {

    // one extra tetrahedron should have been created by this flip
    assert_condition(_tetrahedra_size == 9);

    // convenient names used in documentation figure
    const uint_fast32_t v[5] = {0, 1, 2, 3, 4};
    const uint_fast32_t tn[3] = {0, 1, 8};
    const uint_fast32_t ngb[6] = {2, 3, 4, 5, 6, 7};

    // check the individual tetrahedra
    assert_tetrahedron_right(tn[0], v[0], v[1], v[2], v[4], tn[2], tn[1],
                             ngb[5], ngb[4], 3, 3, 0, 0);
    assert_neighbour_right(ngb[5], tn[0], 2);
    assert_neighbour_right(ngb[4], tn[0], 3);

    assert_tetrahedron_right(tn[1], v[0], v[4], v[2], v[3], tn[2], ngb[3],
                             ngb[2], tn[0], 1, 0, 0, 1);
    assert_neighbour_right(ngb[3], tn[1], 1);
    assert_neighbour_right(ngb[2], tn[1], 2);

    assert_tetrahedron_right(tn[2], v[4], v[1], v[2], v[3], ngb[0], tn[1],
                             ngb[1], tn[0], 0, 0, 0, 0);
    assert_neighbour_right(ngb[0], tn[2], 0);
    assert_neighbour_right(ngb[1], tn[2], 2);

    return;
  }
  case NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP: {

    // no new tetrahedra should have been created by this flip
    assert_condition(_tetrahedra_size == 9);
    // one free spot should have opened up
    assert_condition(_free_size == 1);
    assert_condition(_free_tetrahedra[0] == 2);

    // convenient names used in documentation figure
    const uint_fast32_t v[5] = {0, 1, 2, 3, 4};
    const uint_fast32_t tn[2] = {0, 1};
    const uint_fast32_t ngb[6] = {3, 4, 5, 6, 7, 8};

    // check the individual tetrahedra
    assert_tetrahedron_right(tn[0], v[0], v[1], v[2], v[3], ngb[0], ngb[3],
                             tn[1], ngb[4], 0, 0, 3, 0);
    assert_neighbour_right(ngb[0], tn[0], 0);
    assert_neighbour_right(ngb[3], tn[0], 1);
    assert_neighbour_right(ngb[4], tn[0], 3);

    assert_tetrahedron_right(tn[1], v[0], v[1], v[3], v[4], ngb[1], ngb[2],
                             ngb[5], tn[0], 0, 0, 0, 2);
    assert_neighbour_right(ngb[1], tn[1], 0);
    assert_neighbour_right(ngb[2], tn[1], 1);
    assert_neighbour_right(ngb[5], tn[1], 2);

    return;
  }
  case NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP: {

    // no new tetrahedra should have been created by this flip
    assert_condition(_tetrahedra_size == 12);

    // convenient names used in documentation figure
    const uint_fast32_t v[6] = {0, 1, 2, 3, 4, 5};
    const uint_fast32_t tn[4] = {0, 1, 2, 3};
    const uint_fast32_t ngb[8] = {4, 5, 6, 7, 8, 9, 10, 11};

    // check the individual tetrahedra
    assert_tetrahedron_right(tn[0], v[0], v[3], v[5], v[2], tn[1], ngb[7],
                             ngb[3], tn[2], 0, 0, 0, 3);
    assert_neighbour_right(ngb[7], tn[0], 1);
    assert_neighbour_right(ngb[3], tn[0], 2);

    assert_tetrahedron_right(tn[1], v[1], v[5], v[3], v[2], tn[0], ngb[0],
                             ngb[4], tn[3], 0, 0, 0, 3);
    assert_neighbour_right(ngb[0], tn[1], 1);
    assert_neighbour_right(ngb[4], tn[1], 2);

    assert_tetrahedron_right(tn[2], v[0], v[5], v[3], v[4], tn[3], ngb[2],
                             ngb[6], tn[0], 0, 0, 0, 3);
    assert_neighbour_right(ngb[2], tn[2], 1);
    assert_neighbour_right(ngb[6], tn[2], 2);

    assert_tetrahedron_right(tn[3], v[1], v[3], v[5], v[4], tn[2], ngb[5],
                             ngb[1], tn[1], 0, 0, 0, 3);
    assert_neighbour_right(ngb[5], tn[3], 1);
    assert_neighbour_right(ngb[1], tn[3], 2);

    return;
  }
  }
}

/**
 * @brief Check if NewVoronoiCellConstructor::get_positive_permutation can
 * complete the given permutation of 0123.
 *
 * We define this as a macro so that we get useful line numbers in our error
 * messages (i.e. the line where this macro is called, and not the line in this
 * macro where it goes wrong).
 *
 * @param a First index.
 * @param b Second index.
 * @param c Third index.
 * @param d Fourth index.
 */
#define check_permutation(a, b, c, d)                                          \
  {                                                                            \
    unsigned char v[4] = {a, b, 4, 4};                                         \
    NewVoronoiCellConstructor::get_positive_permutation(v);                    \
    assert_condition(v[2] == c);                                               \
    assert_condition(v[3] == d);                                               \
    assert_condition(NewVoronoiCellConstructor::positive_permutation(          \
                         v[0], v[1], v[2], v[3]) == true);                     \
  }

/**
 * @brief Get the rescaled coordinate in the range [1,2[ that corresponds to the
 * given 53-bit integer mantissa.
 *
 * @param integer_coordinate Integer mantissa.
 * @return Rescaled coordinate.
 */
inline static double get_rescaled_coordinate(unsigned long integer_coordinate) {
  ExactGeometricTests::binary_double real_value;
  // initialize as a value in the correct range to get the correct sign bit
  // and exponent
  real_value.dvalue = 1.1;
  // now overwrite the mantissa
  real_value.parts.mantissa = integer_coordinate;
  // return the new value
  return real_value.dvalue;
}

/**
 * @brief Get the rescaled coordinates in the range [1,2[ that correspond to the
 * given 53-bit integer coordinates.
 *
 * @param integer_coordinates Integer coordinates.
 * @return Rescaled coordinates.
 */
inline static CoordinateVector<> get_rescaled_coordinates(
    const CoordinateVector< unsigned long > &integer_coordinates) {
  return CoordinateVector<>(get_rescaled_coordinate(integer_coordinates.x()),
                            get_rescaled_coordinate(integer_coordinates.y()),
                            get_rescaled_coordinate(integer_coordinates.z()));
}

/**
 * @brief Unit test for the NewVoronoiGrid class.
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
    NewVoronoiTetrahedron tetrahedron(0, 1, 2, 3);
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
  /// part 1: initial tetrahedron without any other insertions
  {
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    std::vector< CoordinateVector<> > positions(1);
    positions[0] = CoordinateVector<>(0.25, 0.25, 0.25);
    NewVoronoiCellConstructor cell_constructor;
    cell_constructor.setup(0, positions, box, positions, box, false);
    NewVoronoiCell cell = cell_constructor.get_cell(box, positions);

    std::vector< CoordinateVector<> > full_volume_positions(4);
    full_volume_positions[0] = CoordinateVector<>(6.26786);
    full_volume_positions[1] = CoordinateVector<>(-8.125, 3.5, 3.5);
    full_volume_positions[2] = CoordinateVector<>(3.5, -8.125, 3.5);
    full_volume_positions[3] = CoordinateVector<>(3.5, 3.5, -8.125);
    NewVoronoiTetrahedron full_volume_tetrahedron(0, 1, 2, 3);
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
    cmac_status("Geometry, part 1: Cell volume computation works!");

    CoordinateVector<> centroid = cell.get_centroid();
    assert_values_equal_rel(centroid.x(), full_centroid.x(), tolerance);
    assert_values_equal_rel(centroid.y(), full_centroid.y(), tolerance);
    assert_values_equal_rel(centroid.z(), full_centroid.z(), tolerance);
    cmac_status("Geometry, part 1: Cell centroid computation works!");

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

    cmac_status("Geometry, part 1: Cell face computation works!");
  }

  /// test for NewVoronoiCellConstructor::find_tetrahedron
  {
    std::vector< CoordinateVector< unsigned long > > positions(4);
    positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    positions[1] = CoordinateVector< unsigned long >(1250, 1500, 1500);
    // point on axis
    positions[2] = CoordinateVector< unsigned long >(1250);
    // point on face
    positions[3] = CoordinateVector< unsigned long >(1250, 1250, 1500);

    std::vector< CoordinateVector<> > rescaled_positions(4);
    rescaled_positions[0] = get_rescaled_coordinates(positions[0]);
    rescaled_positions[1] = get_rescaled_coordinates(positions[1]);
    rescaled_positions[2] = get_rescaled_coordinates(positions[2]);
    rescaled_positions[3] = get_rescaled_coordinates(positions[3]);
    NewVoronoiBox rescaled_box(
        Box<>(CoordinateVector<>(0.), CoordinateVector<>(1.)));

    NewVoronoiCellConstructor cell;
    cell.setup(0, rescaled_positions, rescaled_box, rescaled_positions,
               rescaled_box, false);

    uint_fast32_t tetrahedra[4];
    assert_condition(cell.find_tetrahedron(1, rescaled_box, rescaled_positions,
                                           tetrahedra) == 1);
    assert_condition(tetrahedra[0] == 1);

    assert_condition(cell.find_tetrahedron(2, rescaled_box, rescaled_positions,
                                           tetrahedra) == 3);
    assert_condition(tetrahedra[0] == 3);
    assert_condition(tetrahedra[1] == 2);
    assert_condition(tetrahedra[2] == 1);

    assert_condition(cell.find_tetrahedron(3, rescaled_box, rescaled_positions,
                                           tetrahedra) == 2);
    assert_condition(tetrahedra[0] == 1);
    assert_condition(tetrahedra[1] == 2);

    cmac_status("Find tetrahedron works!");
  }

  /// test permutation routine
  {
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(0, 1, 2, 3) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(0, 2, 3, 1) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(0, 3, 1, 2) == true);

    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(1, 0, 3, 2) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(1, 2, 0, 3) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(1, 3, 2, 0) == true);

    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(2, 0, 1, 3) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(2, 1, 3, 0) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(2, 3, 0, 1) == true);

    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(3, 0, 2, 1) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(3, 1, 0, 2) == true);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(3, 2, 1, 0) == true);

    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(0, 1, 3, 2) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(0, 2, 1, 3) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(0, 3, 2, 1) == false);

    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(1, 0, 2, 3) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(1, 2, 3, 0) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(1, 3, 0, 2) == false);

    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(2, 0, 3, 1) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(2, 1, 0, 3) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(2, 3, 1, 0) == false);

    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(3, 0, 1, 2) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(3, 1, 2, 0) == false);
    assert_condition(
        NewVoronoiCellConstructor::positive_permutation(3, 2, 0, 1) == false);

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
    NewVoronoiCellConstructor cell;

    cell.setup_test(NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP);
    uint_fast32_t tn[4];
    cell.one_to_four_flip(4, 0, tn);
    assert_condition(tn[0] == 0);
    assert_condition(tn[1] == 5);
    assert_condition(tn[2] == 6);
    assert_condition(tn[3] == 7);
    cell.check_test(NEWVORONOICELL_TEST_ONE_TO_FOUR_FLIP);

    cmac_status("1 to 4 flip works!");
  }

  /// test 2 to 6 flip
  {
    NewVoronoiCellConstructor cell;

    cell.setup_test(NEWVORONOICELL_TEST_TWO_TO_SIX_FLIP);
    uint_fast32_t tetrahedra[2] = {0, 1};
    uint_fast32_t tn[6];
    cell.two_to_six_flip(5, tetrahedra, tn);
    assert_condition(tn[0] == 0);
    assert_condition(tn[1] == 1);
    assert_condition(tn[2] == 8);
    assert_condition(tn[3] == 9);
    assert_condition(tn[4] == 10);
    assert_condition(tn[5] == 11);
    cell.check_test(NEWVORONOICELL_TEST_TWO_TO_SIX_FLIP);

    cmac_status("2 to 6 flip works!");
  }

  /// test n to 2n flip
  {
    NewVoronoiCellConstructor cell;

    cell.setup_test(NEWVORONOICELL_TEST_N_TO_2N_FLIP);
    uint_fast32_t tetrahedra[5] = {0, 1, 2, 3, 4};
    uint_fast32_t tn[2 * UCHAR_MAX];
    cell.n_to_2n_flip(7, tetrahedra, 5, tn);
    assert_condition(tn[0] == 0);
    assert_condition(tn[1] == 1);
    assert_condition(tn[2] == 2);
    assert_condition(tn[3] == 3);
    assert_condition(tn[4] == 4);
    assert_condition(tn[5] == 15);
    assert_condition(tn[6] == 16);
    assert_condition(tn[7] == 17);
    assert_condition(tn[8] == 18);
    assert_condition(tn[9] == 19);
    cell.check_test(NEWVORONOICELL_TEST_N_TO_2N_FLIP);

    cmac_status("n to 2n flip works!");
  }

  /// test 2 to 3 flip
  {
    NewVoronoiCellConstructor cell;

    cell.setup_test(NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP);
    uint_fast32_t tn[3];
    cell.two_to_three_flip(0, 1, 2, 3, tn);
    assert_condition(tn[0] == 0);
    assert_condition(tn[1] == 1);
    assert_condition(tn[2] == 8);
    cell.check_test(NEWVORONOICELL_TEST_TWO_TO_THREE_FLIP);

    cmac_status("2 to 3 flip works!");
  }

  /// test 3 to 2 flip
  {
    NewVoronoiCellConstructor cell;

    cell.setup_test(NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP);
    uint_fast32_t tn[2];
    cell.three_to_two_flip(0, 1, 2, tn);
    assert_condition(tn[0] == 0);
    assert_condition(tn[1] == 1);
    cell.check_test(NEWVORONOICELL_TEST_THREE_TO_TWO_FLIP);

    cmac_status("3 to 2 flip works!");
  }

  /// test 4 to 4 flip
  {
    NewVoronoiCellConstructor cell;

    cell.setup_test(NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP);
    uint_fast32_t tn[4];
    cell.four_to_four_flip(0, 1, 2, 3, tn);
    assert_condition(tn[0] == 0);
    assert_condition(tn[1] == 1);
    assert_condition(tn[2] == 2);
    assert_condition(tn[3] == 3);
    cell.check_test(NEWVORONOICELL_TEST_FOUR_TO_FOUR_FLIP);

    cmac_status("4 to 4 flip works!");
  }

  /// tests for NewVoronoiCellConstructor::intersect
  /// simple insertion with a single 1 to 4 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > integer_positions(2);
    integer_positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    integer_positions[1] = CoordinateVector< unsigned long >(1250, 1500, 1500);
    std::vector< CoordinateVector<> > real_positions(2);
    real_positions[0] = get_rescaled_coordinates(integer_positions[0]);
    real_positions[1] = get_rescaled_coordinates(integer_positions[1]);
    NewVoronoiBox real_voronoi_box(
        Box<>(get_rescaled_coordinates(box_anchor),
              get_rescaled_coordinates(box_anchor + box_sides) -
                  get_rescaled_coordinates(box_anchor)));

    NewVoronoiCellConstructor cell;
    cell.setup(0, real_positions, real_voronoi_box, real_positions,
               real_voronoi_box, false);

    cell.intersect(1, real_voronoi_box, real_positions, real_voronoi_box,
                   real_positions);
    cell.check_empty_circumsphere(real_voronoi_box, real_positions);

    std::ofstream ofile("new_voronoi_cell_1_to_4.txt");
    cell.print_tetrahedra(ofile, real_voronoi_box, real_positions);

    cmac_status("Simple 1 to 4 flip insertion worked!");
  }
  /// simple insertion with a single 2 to 6 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > integer_positions(2);
    integer_positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    integer_positions[1] = CoordinateVector< unsigned long >(1250, 1250, 1500);

    std::vector< CoordinateVector<> > real_positions(2);
    real_positions[0] = get_rescaled_coordinates(integer_positions[0]);
    real_positions[1] = get_rescaled_coordinates(integer_positions[1]);
    NewVoronoiBox real_voronoi_box(
        Box<>(get_rescaled_coordinates(box_anchor),
              get_rescaled_coordinates(box_anchor + box_sides) -
                  get_rescaled_coordinates(box_anchor)));

    NewVoronoiCellConstructor cell;
    cell.setup(0, real_positions, real_voronoi_box, real_positions,
               real_voronoi_box, false);

    cell.intersect(1, real_voronoi_box, real_positions, real_voronoi_box,
                   real_positions);
    cell.check_empty_circumsphere(real_voronoi_box, real_positions);

    std::ofstream ofile("new_voronoi_cell_2_to_6.txt");
    cell.print_tetrahedra(ofile, real_voronoi_box, real_positions);

    cmac_status("Simple 2 to 6 flip insertion worked!");
  }
  /// simple insertion with a single n to 2n flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > integer_positions(2);
    integer_positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    integer_positions[1] = CoordinateVector< unsigned long >(1250, 1250, 1250);

    std::vector< CoordinateVector<> > real_positions(2);
    real_positions[0] = get_rescaled_coordinates(integer_positions[0]);
    real_positions[1] = get_rescaled_coordinates(integer_positions[1]);
    NewVoronoiBox real_voronoi_box(
        Box<>(get_rescaled_coordinates(box_anchor),
              get_rescaled_coordinates(box_anchor + box_sides) -
                  get_rescaled_coordinates(box_anchor)));

    NewVoronoiCellConstructor cell;
    cell.setup(0, real_positions, real_voronoi_box, real_positions,
               real_voronoi_box, false);

    cell.intersect(1, real_voronoi_box, real_positions, real_voronoi_box,
                   real_positions);
    cell.check_empty_circumsphere(real_voronoi_box, real_positions);

    std::ofstream ofile("new_voronoi_cell_n_to_2n.txt");
    cell.print_tetrahedra(ofile, real_voronoi_box, real_positions);

    cmac_status("Simple n to 2n flip insertion worked!");
  }
  /// 1 to 4 flip insertion with a single 2 to 3 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > integer_positions(2);
    integer_positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    integer_positions[1] = CoordinateVector< unsigned long >(1005, 1500, 1250);

    std::vector< CoordinateVector<> > real_positions(2);
    real_positions[0] = get_rescaled_coordinates(integer_positions[0]);
    real_positions[1] = get_rescaled_coordinates(integer_positions[1]);
    NewVoronoiBox real_voronoi_box(
        Box<>(get_rescaled_coordinates(box_anchor),
              get_rescaled_coordinates(box_anchor + box_sides) -
                  get_rescaled_coordinates(box_anchor)));

    NewVoronoiCellConstructor cell;
    cell.setup(0, real_positions, real_voronoi_box, real_positions,
               real_voronoi_box, false);

    cell.intersect(1, real_voronoi_box, real_positions, real_voronoi_box,
                   real_positions);
    cell.check_empty_circumsphere(real_voronoi_box, real_positions);

    std::ofstream ofile("new_voronoi_cell_2_to_3.txt");
    cell.print_tetrahedra(ofile, real_voronoi_box, real_positions);
    ofile << real_positions[1].x() << "\t" << real_positions[1].y() << "\t"
          << real_positions[1].z() << "\n";
    ofile.close();

    cmac_status("Simple insertion + 2 to 3 flip worked!");
  }
  /// 2 insertions with a 3 to 2 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > integer_positions(3);
    integer_positions[0] = CoordinateVector< unsigned long >(1500);
    integer_positions[1] = CoordinateVector< unsigned long >(1782, 1393, 1839);
    integer_positions[2] = CoordinateVector< unsigned long >(1197, 1910, 1797);

    std::vector< CoordinateVector<> > real_positions(3);
    real_positions[0] = get_rescaled_coordinates(integer_positions[0]);
    real_positions[1] = get_rescaled_coordinates(integer_positions[1]);
    real_positions[2] = get_rescaled_coordinates(integer_positions[2]);
    NewVoronoiBox real_voronoi_box(
        Box<>(get_rescaled_coordinates(box_anchor),
              get_rescaled_coordinates(box_anchor + box_sides) -
                  get_rescaled_coordinates(box_anchor)));

    NewVoronoiCellConstructor cell;
    cell.setup(0, real_positions, real_voronoi_box, real_positions,
               real_voronoi_box, false);

    for (uint_fast8_t i = 1; i < 3; ++i) {
      cell.intersect(i, real_voronoi_box, real_positions, real_voronoi_box,
                     real_positions);
      cell.check_empty_circumsphere(real_voronoi_box, real_positions);
    }

    std::ofstream ofile("new_voronoi_cell_3_to_2.txt");
    cell.print_tetrahedra(ofile, real_voronoi_box, real_positions);
    ofile << real_positions[1].x() << "\t" << real_positions[1].y() << "\t"
          << real_positions[1].z() << "\n";
    ofile.close();

    cmac_status("Test with 3 to 2 flip worked!");
  }
  /// 1 to 4 flip insertion with a 4 to 4 flip
  {
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > integer_positions(2);
    integer_positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    integer_positions[1] = CoordinateVector< unsigned long >(1500, 1750, 1500);

    std::vector< CoordinateVector<> > real_positions(2);
    real_positions[0] = get_rescaled_coordinates(integer_positions[0]);
    real_positions[1] = get_rescaled_coordinates(integer_positions[1]);
    NewVoronoiBox real_voronoi_box(
        Box<>(get_rescaled_coordinates(box_anchor),
              get_rescaled_coordinates(box_anchor + box_sides) -
                  get_rescaled_coordinates(box_anchor)));

    NewVoronoiCellConstructor cell;
    cell.setup(0, real_positions, real_voronoi_box, real_positions,
               real_voronoi_box, false);

    cell.intersect(1, real_voronoi_box, real_positions, real_voronoi_box,
                   real_positions);
    cell.check_empty_circumsphere(real_voronoi_box, real_positions);

    std::ofstream ofile("new_voronoi_cell_4_to_4.txt");
    cell.print_tetrahedra(ofile, real_voronoi_box, real_positions);
    ofile << real_positions[1].x() << "\t" << real_positions[1].y() << "\t"
          << real_positions[1].z() << "\n";
    ofile.close();

    cmac_status("Simple insertion + 4 to 4 flip worked!");
  }

  /// test geometrical routines
  /// part 2: initial cell with reflective copies
  {
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    std::vector< CoordinateVector<> > positions(1);
    positions[0] = CoordinateVector<>(0.25, 0.25, 0.25);
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > long_positions(1);
    long_positions[0] = CoordinateVector< unsigned long >(1250);

    std::vector< CoordinateVector<> > real_positions(1);
    real_positions[0] = get_rescaled_coordinates(long_positions[0]);
    NewVoronoiBox real_voronoi_box(
        Box<>(get_rescaled_coordinates(box_anchor),
              get_rescaled_coordinates(box_anchor + box_sides) -
                  get_rescaled_coordinates(box_anchor)));

    NewVoronoiCellConstructor cell_constructor;
    cell_constructor.setup(0, positions, box, real_positions, real_voronoi_box,
                           true);

    NewVoronoiCell cell = cell_constructor.get_cell(box, positions);

    const double tolerance = 1.e-15;

    double volume = cell.get_volume();
    assert_values_equal_rel(volume, 1., tolerance);
    cmac_status("Geometry, part 2: Cell volume computation works!");

    CoordinateVector<> centroid = cell.get_centroid();
    assert_values_equal_rel(centroid.x(), 0.5, tolerance);
    assert_values_equal_rel(centroid.y(), 0.5, tolerance);
    assert_values_equal_rel(centroid.z(), 0.5, tolerance);
    cmac_status("Geometry, part 2: Cell centroid computation works!");

    std::vector< VoronoiFace > faces = cell.get_faces();
    assert_condition(faces.size() == 6);

    assert_values_equal_rel(faces[3].get_surface_area(), 1., tolerance);
    CoordinateVector<> midpoint = faces[3].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0., tolerance);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[3].get_neighbour() == NEWVORONOICELL_BOX_LEFT);
    std::vector< CoordinateVector<> > vertices = faces[3].get_vertices();
    assert_condition(vertices.size() == 4);
    assert_condition(vertices[0] == CoordinateVector<>(0., 0., 1.));
    assert_condition(vertices[1] == CoordinateVector<>(0., 0., 0.));
    assert_condition(vertices[2] == CoordinateVector<>(0., 1., 0.));
    assert_condition(vertices[3] == CoordinateVector<>(0., 1., 1.));

    assert_values_equal_rel(faces[0].get_surface_area(), 1., tolerance);
    midpoint = faces[0].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 1., tolerance);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[0].get_neighbour() == NEWVORONOICELL_BOX_RIGHT);
    vertices = faces[0].get_vertices();
    assert_condition(vertices.size() == 4);
    assert_condition(vertices[0] == CoordinateVector<>(1., 0., 1.));
    assert_condition(vertices[1] == CoordinateVector<>(1., 0., 0.));
    assert_condition(vertices[2] == CoordinateVector<>(1., 1., 0.));
    assert_condition(vertices[3] == CoordinateVector<>(1., 1., 1.));

    assert_values_equal_rel(faces[2].get_surface_area(), 1., tolerance);
    midpoint = faces[2].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.y(), 0., tolerance);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[2].get_neighbour() == NEWVORONOICELL_BOX_FRONT);
    vertices = faces[2].get_vertices();
    assert_condition(vertices.size() == 4);
    assert_condition(vertices[0] == CoordinateVector<>(1., 0., 1.));
    assert_condition(vertices[1] == CoordinateVector<>(0., 0., 1.));
    assert_condition(vertices[2] == CoordinateVector<>(0., 0., 0.));
    assert_condition(vertices[3] == CoordinateVector<>(1., 0., 0.));

    assert_values_equal_rel(faces[5].get_surface_area(), 1., tolerance);
    midpoint = faces[5].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.y(), 1., tolerance);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[5].get_neighbour() == NEWVORONOICELL_BOX_BACK);
    vertices = faces[5].get_vertices();
    assert_condition(vertices.size() == 4);
    assert_condition(vertices[0] == CoordinateVector<>(1., 1., 1.));
    assert_condition(vertices[1] == CoordinateVector<>(1., 1., 0.));
    assert_condition(vertices[2] == CoordinateVector<>(0., 1., 0.));
    assert_condition(vertices[3] == CoordinateVector<>(0., 1., 1.));

    assert_values_equal_rel(faces[4].get_surface_area(), 1., tolerance);
    midpoint = faces[4].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.z(), 0., tolerance);
    assert_condition(faces[4].get_neighbour() == NEWVORONOICELL_BOX_BOTTOM);
    vertices = faces[4].get_vertices();
    assert_condition(vertices.size() == 4);
    assert_condition(vertices[0] == CoordinateVector<>(0., 0., 0.));
    assert_condition(vertices[1] == CoordinateVector<>(1., 0., 0.));
    assert_condition(vertices[2] == CoordinateVector<>(1., 1., 0.));
    assert_condition(vertices[3] == CoordinateVector<>(0., 1., 0.));

    assert_values_equal_rel(faces[1].get_surface_area(), 1., tolerance);
    midpoint = faces[1].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.z(), 1., tolerance);
    assert_condition(faces[1].get_neighbour() == NEWVORONOICELL_BOX_TOP);
    vertices = faces[1].get_vertices();
    assert_condition(vertices.size() == 4);
    assert_condition(vertices[0] == CoordinateVector<>(1., 0., 1.));
    assert_condition(vertices[1] == CoordinateVector<>(1., 1., 1.));
    assert_condition(vertices[2] == CoordinateVector<>(0., 1., 1.));
    assert_condition(vertices[3] == CoordinateVector<>(0., 0., 1.));

    cmac_status("Geometry, part 2: Cell face computation works!");
  }

  /// test NewVoronoiGrid construction: random generators
  {
    const uint_fast32_t ncell = 100;
    std::vector< CoordinateVector<> > positions(ncell);
    for (uint_fast32_t i = 0; i < ncell; ++i) {
      positions[i] = Utilities::random_position();
    }

    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    NewVoronoiGrid grid(positions, box);

    Timer timer;
    timer.start();
    grid.compute_grid();
    timer.stop();

    double time_per_cell = timer.value() / ncell;
    cmac_status("Standard grid construction works (%g s, %g s/cell)!",
                timer.value(), time_per_cell);
  }

  /// test NewVoronoiGrid construction: regular generators
  {
    const uint_fast32_t ncell_1D = 5;
    const uint_fast32_t ncell_2D = ncell_1D * ncell_1D;
    const uint_fast32_t ncell_3D = ncell_2D * ncell_1D;
    std::vector< CoordinateVector<> > positions(ncell_3D);
    for (uint_fast32_t ix = 0; ix < ncell_1D; ++ix) {
      for (uint_fast32_t iy = 0; iy < ncell_1D; ++iy) {
        for (uint_fast32_t iz = 0; iz < ncell_1D; ++iz) {
          positions[ncell_2D * ix + ncell_1D * iy + iz] =
              CoordinateVector<>((ix + 0.5) / ncell_1D, (iy + 0.5) / ncell_1D,
                                 (iz + 0.5) / ncell_1D);
        }
      }
    }

    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    NewVoronoiGrid grid(positions, box);

    Timer timer;
    timer.start();
    grid.compute_grid();
    timer.stop();

    double time_per_cell = timer.value() / ncell_3D;
    cmac_status(
        "Regular (degenerate) grid construction works (%g s, %g s/cell)!",
        timer.value(), time_per_cell);
  }

  return 0;
}
