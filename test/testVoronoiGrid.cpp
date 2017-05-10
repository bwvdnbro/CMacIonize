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
 * @file testVoronoiGrid.cpp
 *
 * @brief Unit test for the VoronoiGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Utilities.hpp"
#include "VoronoiCell.hpp"
#include "VoronoiGrid.hpp"

#include <fstream>

/**
 * @brief Aliases for testcase setups.
 *
 * Most of these tests were copied from the 3D Voronoi grid implementation in
 * SWIFT (https://gitlab.cosma.dur.ac.uk/swift/swiftsim), written by Bert
 * Vandenbroucke.
 * The SWIFT code contains some comments that explain the origin of the PATH
 * naming convention.
 *
 * If you ever wonder why someone would go through all the trouble of
 * implementing all these tests: VORONOITEST_PATH_2_1 actually helped finding a
 * very stupid but nonetheless extremely hard to find bug.
 */
enum {
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, and its first edge intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_0 = 0,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, and its second edge intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_1,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, and its second (and also last) edge
   * intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_2,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, and is also the closest vertex ot the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly returns an error code in this case.
   */
  VORONOITEST_PATH_1_3,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, the second vertex is closer to the plane, and
   * the first edge of the second vertex intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_4_0,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, the second vertex is closer to the plane, and
   * the second edge of the second vertex intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_4_1,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, the second vertex is closer to the plane, and
   * the second edge of the second vertex intersects the plane.
   *
   * This case is identical to VORONOITEST_PATH_1_4_1 above, but we reorder some
   * of the edge connections to force the algorithm to enter an extra loop.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_4_2,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, the second vertex is closer to the plane, and
   * the second edge of the second vertex intersects the plane.
   *
   * The difference between this case and cases VORONOITEST_PATH_1_4_1 and
   * VORONOITEST_PATH_1_4_2 is that the first edge of the second vertex in this
   * case links to the first vertex, which means we end up in a different part
   * of the algorithm.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_4_3,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, the second vertex is closer to the plane, and
   * the third (but not last) edge of the second vertex intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_4_4,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, the second vertex is closer to the plane, and
   * the third (and last) edge of the second vertex intersects the plane.
   *
   * The only difference between this test and VORONOITEST_PATH_1_4_4 is the
   * fact that now the third edge is the last edge of the second vertex.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_1_4_5,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, the second vertex is closer and is the
   * closest vertex to the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly returns an error code in this case.
   */
  VORONOITEST_PATH_1_4_6,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * above the intersecting plane, and the second vertex is on the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects a complicated setup in this case.
   */
  VORONOITEST_PATH_1_5,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, and its first edge intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects a complicated setup in this case.
   */
  VORONOITEST_PATH_2_0,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, and its second edge intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_1,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, and its second (and last) edge intersects the
   * plane.
   *
   * The only difference between this test and VORONOITEST_PATH_2_1 is that in
   * this case the second edge is also the last edge.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_2,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, and is the closest vertex to the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly returns a zero status. Although in this case this would
   * be the entire algorithm.
   */
  VORONOITEST_PATH_2_3,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, the second vertex is closer, and its first
   * edge intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_4_0,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, the second vertex is closer, and its second
   * (but not last) edge intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_4_1,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, the second vertex is closer, and its second
   * (and last) edge intersects the plane.
   *
   * The only difference between this test and VORONOITEST_PATH_2_4_1 is that in
   * this case the second edge is also the last edge.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_4_2,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, the second vertex is closer, and its second
   * edge intersects the plane.
   *
   * The difference between this case and cases VORONOITEST_PATH_2_4_1 and
   * VORONOITEST_PATH_2_4_2 is that the first edge of the second vertex in this
   * case links to the first vertex, which means we end up in a different part
   * of the algorithm.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_4_3,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, the second vertex is closer, and its third
   * (but not last) edge intersects the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_4_4,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, the second vertex is closer, and its third
   * (and last) edge intersects the plane.
   *
   * The only difference between this test and VORONOITEST_PATH_2_4_4 is that in
   * this case the third edge is also the last edge.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects an intersected edge in this case.
   */
  VORONOITEST_PATH_2_4_5,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, and the second vertex is closer and is the
   * closest vertex to the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly returns a zero status. Although in this case this would
   * be the entire algorithm.
   */
  VORONOITEST_PATH_2_4_6,
  /**
   * @brief Test the intersection algorithm when the first vertex of the cell is
   * below the intersecting plane, and the second vertex is on the plane.
   *
   * We only test the first part of the algorithm, and check if the intersection
   * routine correctly detects a complicated setup in this case.
   */
  VORONOITEST_PATH_2_5,
  /**
   * @brief Test the vertex deletion algorithm.
   *
   * This test sets up a cell with six vertices, two of which need to be
   * deleted. We check if the two vertices are deleted, and if the indices in
   * the edge connections between the remaining vertices are updated correctly.
   */
  VORONOITEST_VERTEX_DELETION,
  /**
   * @brief Test the standard constructor by checking that the initial vertex
   * positions are correct.
   */
  VORONOITEST_CONSTRUCTOR,
  /**
   * @brief Test the vertex intersection algorithm, by intersecting a cubic cell
   * with a single neighbouring generator.
   *
   * We use the standard setup for this test, and check if the resulting cell
   * has the correct vertices and connections.
   */
  VORONOITEST_INTERSECTION
};

/**
 * @brief Overwrites the default Voronoi cell initialization and sets the cell
 * variables to values that trigger a certain behaviour that is needed for a
 * unit test.
 *
 * @param testcase ID of a testcase.
 */
void VoronoiCell::setup_variables_for_test(int testcase) {
  switch (testcase) {

  case VORONOITEST_PATH_1_0: {
    _vertices.resize(2);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(2);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint_index(1, 0, 0);
    break;
  }

  case VORONOITEST_PATH_1_1: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(2., 0., 0.);
    _vertices[2] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint(0, 1, 2);
    set_edge_endpoint_index(0, 1, 0);
    _edges[1].resize(3);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 0);
    set_edge_endpoint_index(2, 0, 1);
    break;
  }

  case VORONOITEST_PATH_1_2: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(2., 0., 0.);
    _vertices[2] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(2);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint(0, 1, 2);
    set_edge_endpoint_index(0, 1, 0);
    _edges[1].resize(3);
    _edges[2].resize(3);
    set_edge_endpoint(2, 1, 0);
    set_edge_endpoint_index(2, 0, 1);
    break;
  }

  case VORONOITEST_PATH_1_3: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(2., 0., 0.);
    _vertices[2] = CoordinateVector<>(2., 0., 0.);
    _vertices[3] = CoordinateVector<>(2., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint(0, 1, 2);
    set_edge_endpoint(0, 2, 3);
    break;
  }

  case VORONOITEST_PATH_1_4_0: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.75, 0., 0.);
    _vertices[2] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(2);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 2);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 2);
    set_edge_endpoint_index(1, 0, 0);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 1);
    set_edge_endpoint_index(2, 0, 0);
    break;
  }

  case VORONOITEST_PATH_1_4_1: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.75, 0., 0.);
    _vertices[2] = CoordinateVector<>(2., 0., 0.);
    _vertices[3] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 3);
    _edges[1].resize(4);
    set_edge_endpoint(1, 0, 2);
    set_edge_endpoint(1, 1, 3);
    set_edge_endpoint(1, 3, 0);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 1, 0);
    set_edge_endpoint_index(1, 3, 0);
    _edges[2].resize(3);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 1);
    break;
  }

  case VORONOITEST_PATH_1_4_2: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.75, 0., 0.);
    _vertices[2] = CoordinateVector<>(2., 0., 0.);
    _vertices[3] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 3);
    // the line below is the only difference between this test and
    // VORONOITEST_PATH_1_4_2 (+ the indices of the third edge)
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 2);
    set_edge_endpoint(1, 1, 3);
    set_edge_endpoint(1, 2, 0);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 1, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 1);
    break;
  }

  case VORONOITEST_PATH_1_4_3: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.75, 0., 0.);
    _vertices[2] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(2);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(3);
    // the line below is the only difference between this test and
    // VORONOITEST_PATH_1_4_0
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 1);
    set_edge_endpoint_index(2, 0, 1);
    break;
  }

  case VORONOITEST_PATH_1_4_4: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.75, 0., 0.);
    _vertices[2] = CoordinateVector<>(2., 0., 0.);
    _vertices[3] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(4);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    set_edge_endpoint(1, 2, 3);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 1);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 0);
    break;
  }

  case VORONOITEST_PATH_1_4_5: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.75, 0., 0.);
    _vertices[2] = CoordinateVector<>(2., 0., 0.);
    _vertices[3] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    // the line below is the only difference between this test and
    // VORONOITEST_PATH_1_4_4
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    set_edge_endpoint(1, 2, 3);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 1);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 0);
    break;
  }

  case VORONOITEST_PATH_1_4_6: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.75, 0., 0.);
    _vertices[2] = CoordinateVector<>(2., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(2);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    set_edge_endpoint_index(1, 0, 0);
    _edges[2].resize(3);
    break;
  }

  case VORONOITEST_PATH_1_5: {
    _vertices.resize(2);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.5, 0., 0.);
    _edges.resize(2);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint_index(1, 0, 0);
    break;
  }

  case VORONOITEST_PATH_2_0: {
    _vertices.resize(2);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(2);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint_index(1, 0, 0);
    break;
  }

  case VORONOITEST_PATH_2_1: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-2., 0., 0.);
    _vertices[2] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint(0, 1, 2);
    set_edge_endpoint_index(0, 1, 0);
    _edges[1].resize(3);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 0);
    set_edge_endpoint_index(2, 0, 1);
    break;
  }

  case VORONOITEST_PATH_2_2: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-2., 0., 0.);
    _vertices[2] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(3);
    // the line below is the only difference between this test and
    // VORONOITEST_PATH_2_1
    _edges[0].resize(2);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint(0, 1, 2);
    set_edge_endpoint_index(0, 1, 0);
    _edges[1].resize(3);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 0);
    set_edge_endpoint_index(2, 0, 1);
    break;
  }

  case VORONOITEST_PATH_2_3: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-2., 0., 0.);
    _vertices[2] = CoordinateVector<>(-2., 0., 0.);
    _vertices[3] = CoordinateVector<>(-2., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint(0, 1, 2);
    set_edge_endpoint(0, 2, 3);
    _edges[1].resize(3);
    _edges[2].resize(3);
    _edges[3].resize(3);
    break;
  }

  case VORONOITEST_PATH_2_4_0: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-0.5, 0., 0.);
    _vertices[2] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 2);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 2);
    set_edge_endpoint(1, 2, 0);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 1);
    set_edge_endpoint_index(2, 0, 0);
    break;
  }

  case VORONOITEST_PATH_2_4_1: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-0.5, 0., 0.);
    _vertices[2] = CoordinateVector<>(-2., 0., 0.);
    _vertices[3] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 3);
    _edges[1].resize(4);
    set_edge_endpoint(1, 0, 2);
    set_edge_endpoint(1, 1, 3);
    set_edge_endpoint(1, 3, 0);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 1, 0);
    set_edge_endpoint_index(1, 3, 0);
    _edges[2].resize(3);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 1);
    break;
  }

  case VORONOITEST_PATH_2_4_2: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-0.5, 0., 0.);
    _vertices[2] = CoordinateVector<>(-2., 0., 0.);
    _vertices[3] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 3);
    // the line below is the only difference between this test and
    // VORONOITEST_PATH_2_4_2 (+ the indices of the third edge)
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 2);
    set_edge_endpoint(1, 1, 3);
    set_edge_endpoint(1, 2, 0);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 1, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 1);
    break;
  }

  case VORONOITEST_PATH_2_4_3: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-0.5, 0., 0.);
    _vertices[2] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 1, 0);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 1);
    set_edge_endpoint_index(2, 0, 1);
    break;
  }

  case VORONOITEST_PATH_2_4_4: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-0.5, 0., 0.);
    _vertices[2] = CoordinateVector<>(-2., 0., 0.);
    _vertices[3] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(4);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    set_edge_endpoint(1, 2, 3);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 2);
    break;
  }

  case VORONOITEST_PATH_2_4_5: {
    _vertices.resize(4);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-0.5, 0., 0.);
    _vertices[2] = CoordinateVector<>(-2., 0., 0.);
    _vertices[3] = CoordinateVector<>(1., 0., 0.);
    _edges.resize(4);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    // the line below is the only difference between this test and
    // VORONOITEST_PATH_2_4_4
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    set_edge_endpoint(1, 2, 3);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 1);
    set_edge_endpoint_index(3, 0, 2);
    break;
  }

  case VORONOITEST_PATH_2_4_6: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(-0.5, 0., 0.);
    _vertices[2] = CoordinateVector<>(-2., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(2);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint(1, 1, 2);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 1, 0);
    _edges[2].resize(3);
    break;
  }

  case VORONOITEST_PATH_2_5: {
    _vertices.resize(2);
    _vertices[0] = CoordinateVector<>(-1., 0., 0.);
    _vertices[1] = CoordinateVector<>(0.5, 0., 0.);
    _edges.resize(2);
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 0);
    set_edge_endpoint_index(1, 0, 0);
    break;
  }

  case VORONOITEST_VERTEX_DELETION: {
    _vertices.resize(6);
    // the actual positions of the vertices are completely irrelevant for this
    // test, so we don't set them
    _edges.resize(6);
    // the first vertex needs to be deleted, and is connected to the fourth
    // vertex
    _edges[0].resize(3);
    set_edge_endpoint(0, 0, 3);
    set_edge_endpoint(0, 1, -1);
    set_edge_endpoint(0, 2, -1);
    set_edge_endpoint_index(0, 0, 0);
    _edges[1].resize(3);
    set_edge_endpoint(1, 0, 2);
    set_edge_endpoint(1, 1, 4);
    set_edge_endpoint(1, 2, 5);
    set_edge_endpoint_index(1, 0, 0);
    set_edge_endpoint_index(1, 1, 0);
    set_edge_endpoint_index(1, 2, 0);
    _edges[2].resize(3);
    set_edge_endpoint(2, 0, 1);
    set_edge_endpoint(2, 1, 5);
    set_edge_endpoint(2, 2, 4);
    set_edge_endpoint_index(2, 0, 0);
    set_edge_endpoint_index(2, 1, 1);
    set_edge_endpoint_index(2, 2, 1);
    _edges[3].resize(3);
    set_edge_endpoint(3, 0, 0);
    set_edge_endpoint(3, 1, -1);
    set_edge_endpoint(3, 2, -1);
    set_edge_endpoint_index(3, 0, 1);
    _edges[4].resize(3);
    set_edge_endpoint(4, 0, 1);
    set_edge_endpoint(4, 1, 2);
    set_edge_endpoint(4, 2, 5);
    set_edge_endpoint_index(4, 0, 1);
    set_edge_endpoint_index(4, 1, 2);
    set_edge_endpoint_index(4, 2, 2);
    _edges[5].resize(3);
    set_edge_endpoint(5, 0, 1);
    set_edge_endpoint(5, 1, 2);
    set_edge_endpoint(5, 2, 4);
    set_edge_endpoint_index(5, 0, 2);
    set_edge_endpoint_index(5, 1, 1);
    set_edge_endpoint_index(5, 2, 2);
    break;
  }
  }
}

/**
 * @brief Check that the Voronoi cell contains the right vertices and edge
 * connections after a routine has changed its contents.
 *
 * @param testcase ID of a testcase.
 */
void VoronoiCell::check_variables_after_test(int testcase) {

  switch (testcase) {
  case VORONOITEST_VERTEX_DELETION: {
    assert_condition(_vertices.size() == 4);
    assert_condition(_edges.size() == 4);

    assert_condition(_edges[0].size() == 3);
    assert_condition(get_edge_endpoint(0, 0) == 1);
    assert_condition(get_edge_endpoint(0, 1) == 2);
    assert_condition(get_edge_endpoint(0, 2) == 3);
    assert_condition(get_edge_endpoint_index(0, 0) == 0);
    assert_condition(get_edge_endpoint_index(0, 1) == 0);
    assert_condition(get_edge_endpoint_index(0, 2) == 0);

    assert_condition(_edges[1].size() == 3);
    assert_condition(get_edge_endpoint(1, 0) == 0);
    assert_condition(get_edge_endpoint(1, 1) == 3);
    assert_condition(get_edge_endpoint(1, 2) == 2);
    assert_condition(get_edge_endpoint_index(1, 0) == 0);
    assert_condition(get_edge_endpoint_index(1, 1) == 1);
    assert_condition(get_edge_endpoint_index(1, 2) == 1);

    assert_condition(_edges[2].size() == 3);
    assert_condition(get_edge_endpoint(2, 0) == 0);
    assert_condition(get_edge_endpoint(2, 1) == 1);
    assert_condition(get_edge_endpoint(2, 2) == 3);
    assert_condition(get_edge_endpoint_index(2, 0) == 1);
    assert_condition(get_edge_endpoint_index(2, 1) == 2);
    assert_condition(get_edge_endpoint_index(2, 2) == 2);

    assert_condition(_edges[3].size() == 3);
    assert_condition(get_edge_endpoint(3, 0) == 0);
    assert_condition(get_edge_endpoint(3, 1) == 1);
    assert_condition(get_edge_endpoint(3, 2) == 2);
    assert_condition(get_edge_endpoint_index(3, 0) == 2);
    assert_condition(get_edge_endpoint_index(3, 1) == 1);
    assert_condition(get_edge_endpoint_index(3, 2) == 2);

    break;
  }

  case VORONOITEST_CONSTRUCTOR: {
    assert_condition(_vertices.size() == 8);
    assert_condition(_edges.size() == 8);

    assert_condition(_vertices[0] == CoordinateVector<>(-0.5, -0.5, -0.5));
    assert_condition(_vertices[1] == CoordinateVector<>(-0.5, -0.5, 0.5));
    assert_condition(_vertices[2] == CoordinateVector<>(-0.5, 0.5, -0.5));
    assert_condition(_vertices[3] == CoordinateVector<>(-0.5, 0.5, 0.5));
    assert_condition(_vertices[4] == CoordinateVector<>(0.5, -0.5, -0.5));
    assert_condition(_vertices[5] == CoordinateVector<>(0.5, -0.5, 0.5));
    assert_condition(_vertices[6] == CoordinateVector<>(0.5, 0.5, -0.5));
    assert_condition(_vertices[7] == CoordinateVector<>(0.5, 0.5, 0.5));

    break;
  }

  case VORONOITEST_INTERSECTION: {
    assert_condition(_vertices.size() == 10);
    assert_condition(_edges.size() == 10);

    assert_condition(_vertices[0] == CoordinateVector<>(-0.5, -0.5, -0.5));
    assert_condition(_vertices[1] == CoordinateVector<>(-0.5, -0.5, 0.5));
    assert_condition(_vertices[2] == CoordinateVector<>(-0.5, 0.5, -0.5));
    assert_condition(_vertices[3] == CoordinateVector<>(-0.5, 0.5, 0.5));
    assert_condition(_vertices[4] == CoordinateVector<>(0.5, 0.5, -0.5));
    assert_condition(_vertices[5] == CoordinateVector<>(0.5, 0.5, 0.5));
    assert_condition(_vertices[6] == CoordinateVector<>(0., -0.5, -0.5));
    assert_condition(_vertices[7] == CoordinateVector<>(0.5, 0., -0.5));
    assert_condition(_vertices[8] == CoordinateVector<>(0.5, 0., 0.5));
    assert_condition(_vertices[9] == CoordinateVector<>(0., -0.5, 0.5));

    assert_condition(_edges[0].size() == 3);
    assert_condition(get_edge_endpoint(0, 0) == 1);
    assert_condition(get_edge_endpoint(0, 1) == 2);
    assert_condition(get_edge_endpoint(0, 2) == 6);
    assert_condition(get_edge_endpoint_index(0, 0) == 0);
    assert_condition(get_edge_endpoint_index(0, 1) == 2);
    assert_condition(get_edge_endpoint_index(0, 2) == 1);
    assert_condition(get_edge_neighbour(0, 0) == VORONOI_BOX_FRONT);
    assert_condition(get_edge_neighbour(0, 1) == VORONOI_BOX_LEFT);
    assert_condition(get_edge_neighbour(0, 2) == VORONOI_BOX_BOTTOM);

    assert_condition(_edges[1].size() == 3);
    assert_condition(get_edge_endpoint(1, 0) == 0);
    assert_condition(get_edge_endpoint(1, 1) == 9);
    assert_condition(get_edge_endpoint(1, 2) == 3);
    assert_condition(get_edge_endpoint(1, 0) == 0);
    assert_condition(get_edge_endpoint_index(1, 1) == 1);
    assert_condition(get_edge_endpoint_index(1, 2) == 1);
    assert_condition(get_edge_neighbour(1, 0) == VORONOI_BOX_LEFT);
    assert_condition(get_edge_neighbour(1, 1) == VORONOI_BOX_FRONT);
    assert_condition(get_edge_neighbour(1, 2) == VORONOI_BOX_TOP);

    assert_condition(_edges[2].size() == 3);
    assert_condition(get_edge_endpoint(2, 0) == 3);
    assert_condition(get_edge_endpoint(2, 1) == 4);
    assert_condition(get_edge_endpoint(2, 2) == 0);
    assert_condition(get_edge_endpoint_index(2, 0) == 0);
    assert_condition(get_edge_endpoint_index(2, 1) == 0);
    assert_condition(get_edge_endpoint_index(2, 2) == 1);
    assert_condition(get_edge_neighbour(2, 0) == VORONOI_BOX_LEFT);
    assert_condition(get_edge_neighbour(2, 1) == VORONOI_BOX_BACK);
    assert_condition(get_edge_neighbour(2, 2) == VORONOI_BOX_BOTTOM);

    assert_condition(_edges[3].size() == 3);
    assert_condition(get_edge_endpoint(3, 0) == 2);
    assert_condition(get_edge_endpoint(3, 1) == 1);
    assert_condition(get_edge_endpoint(3, 2) == 5);
    assert_condition(get_edge_endpoint_index(3, 0) == 0);
    assert_condition(get_edge_endpoint_index(3, 1) == 2);
    assert_condition(get_edge_endpoint_index(3, 2) == 0);
    assert_condition(get_edge_neighbour(3, 0) == VORONOI_BOX_BACK);
    assert_condition(get_edge_neighbour(3, 1) == VORONOI_BOX_LEFT);
    assert_condition(get_edge_neighbour(3, 2) == VORONOI_BOX_TOP);

    assert_condition(_edges[4].size() == 3);
    assert_condition(get_edge_endpoint(4, 0) == 2);
    assert_condition(get_edge_endpoint(4, 1) == 5);
    assert_condition(get_edge_endpoint(4, 2) == 7);
    assert_condition(get_edge_endpoint_index(4, 0) == 1);
    assert_condition(get_edge_endpoint_index(4, 1) == 2);
    assert_condition(get_edge_endpoint_index(4, 2) == 1);
    assert_condition(get_edge_neighbour(4, 0) == VORONOI_BOX_BOTTOM);
    assert_condition(get_edge_neighbour(4, 1) == VORONOI_BOX_BACK);
    assert_condition(get_edge_neighbour(4, 2) == VORONOI_BOX_RIGHT);

    assert_condition(_edges[5].size() == 3);
    assert_condition(get_edge_endpoint(5, 0) == 3);
    assert_condition(get_edge_endpoint(5, 1) == 8);
    assert_condition(get_edge_endpoint(5, 2) == 4);
    assert_condition(get_edge_endpoint_index(5, 0) == 2);
    assert_condition(get_edge_endpoint_index(5, 1) == 1);
    assert_condition(get_edge_endpoint_index(5, 2) == 1);
    assert_condition(get_edge_neighbour(5, 0) == VORONOI_BOX_BACK);
    assert_condition(get_edge_neighbour(5, 1) == VORONOI_BOX_TOP);
    assert_condition(get_edge_neighbour(5, 2) == VORONOI_BOX_RIGHT);

    assert_condition(_edges[6].size() == 3);
    assert_condition(get_edge_endpoint(6, 0) == 9);
    assert_condition(get_edge_endpoint(6, 1) == 0);
    assert_condition(get_edge_endpoint(6, 2) == 7);
    assert_condition(get_edge_endpoint_index(6, 0) == 2);
    assert_condition(get_edge_endpoint_index(6, 1) == 2);
    assert_condition(get_edge_endpoint_index(6, 2) == 0);
    assert_condition(get_edge_neighbour(6, 0) == 0);
    assert_condition(get_edge_neighbour(6, 1) == VORONOI_BOX_FRONT);
    assert_condition(get_edge_neighbour(6, 2) == VORONOI_BOX_BOTTOM);

    assert_condition(_edges[7].size() == 3);
    assert_condition(get_edge_endpoint(7, 0) == 6);
    assert_condition(get_edge_endpoint(7, 1) == 4);
    assert_condition(get_edge_endpoint(7, 2) == 8);
    assert_condition(get_edge_endpoint_index(7, 0) == 2);
    assert_condition(get_edge_endpoint_index(7, 1) == 2);
    assert_condition(get_edge_endpoint_index(7, 2) == 0);
    assert_condition(get_edge_neighbour(7, 0) == 0);
    assert_condition(get_edge_neighbour(7, 1) == VORONOI_BOX_BOTTOM);
    assert_condition(get_edge_neighbour(7, 2) == VORONOI_BOX_RIGHT);

    assert_condition(_edges[8].size() == 3);
    assert_condition(get_edge_endpoint(8, 0) == 7);
    assert_condition(get_edge_endpoint(8, 1) == 5);
    assert_condition(get_edge_endpoint(8, 2) == 9);
    assert_condition(get_edge_endpoint_index(8, 0) == 2);
    assert_condition(get_edge_endpoint_index(8, 1) == 1);
    assert_condition(get_edge_endpoint_index(8, 2) == 0);
    assert_condition(get_edge_neighbour(8, 0) == 0);
    assert_condition(get_edge_neighbour(8, 1) == VORONOI_BOX_RIGHT);
    assert_condition(get_edge_neighbour(8, 2) == VORONOI_BOX_TOP);

    assert_condition(_edges[9].size() == 3);
    assert_condition(get_edge_endpoint(9, 0) == 8);
    assert_condition(get_edge_endpoint(9, 1) == 1);
    assert_condition(get_edge_endpoint(9, 2) == 6);
    assert_condition(get_edge_endpoint_index(9, 0) == 2);
    assert_condition(get_edge_endpoint_index(9, 1) == 1);
    assert_condition(get_edge_endpoint_index(9, 2) == 0);
    assert_condition(get_edge_neighbour(9, 0) == 0);
    assert_condition(get_edge_neighbour(9, 1) == VORONOI_BOX_TOP);
    assert_condition(get_edge_neighbour(9, 2) == VORONOI_BOX_FRONT);

    break;
  }
  }
}

/**
 * @brief Unit test for the VoronoiGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// test static geometric functions in VoronoiCell
  {
    assert_condition(VoronoiCell::volume_tetrahedron(
                         CoordinateVector<>(0., 0., 0.),
                         CoordinateVector<>(0., 0., 1.),
                         CoordinateVector<>(0., 1., 0.),
                         CoordinateVector<>(1., 0., 0.)) == 1. / 6.);

    assert_condition(
        VoronoiCell::centroid_tetrahedron(
            CoordinateVector<>(0., 0., 0.), CoordinateVector<>(0., 0., 1.),
            CoordinateVector<>(0., 1., 0.), CoordinateVector<>(1., 0., 0.)) ==
        CoordinateVector<>(0.25, 0.25, 0.25));

    assert_condition(VoronoiCell::surface_area_triangle(
                         CoordinateVector<>(0., 0., 0.),
                         CoordinateVector<>(1., 0., 0.),
                         CoordinateVector<>(0., 1., 0.)) == 0.5);

    assert_condition(
        VoronoiCell::midpoint_triangle(CoordinateVector<>(0., 0., 0.),
                                       CoordinateVector<>(1., 0., 0.),
                                       CoordinateVector<>(0., 1., 0.)) ==
        CoordinateVector<>(1. / 3., 1. / 3., 0.));

    // vertex below plane
    auto vertex_test = VoronoiCell::test_vertex(
        CoordinateVector<>(0.5, 1., 0.), CoordinateVector<>(1., 0., 0.), 1.);
    assert_condition(vertex_test.first == -1);
    assert_condition(vertex_test.second == -0.5);

    // vertex above plane
    vertex_test = VoronoiCell::test_vertex(CoordinateVector<>(1.5, 1., 0.),
                                           CoordinateVector<>(1., 0., 0.), 1.);
    assert_condition(vertex_test.first == 1);
    assert_condition(vertex_test.second == 0.5);

    // vertex on plane
    vertex_test = VoronoiCell::test_vertex(CoordinateVector<>(1., 1., 0.),
                                           CoordinateVector<>(1., 0., 0.), 1.);
    assert_condition(vertex_test.first == 0);
    assert_condition(vertex_test.second == 0.);
  }

  /// test Voronoi cell volume, centroid and face calculation
  {
    Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    VoronoiCell cell(CoordinateVector<>(0.5), box);
    assert_condition(cell.get_max_radius_squared() == 0.75);
    cell.finalize();
    assert_values_equal_rel(cell.get_volume(), 1., 1.e-16);
    assert_condition(cell.get_centroid() == CoordinateVector<>(0.5));

    auto faces = cell.get_faces();
    assert_condition(VoronoiCell::get_face_neighbour(faces[0]) ==
                     VORONOI_BOX_FRONT);
    assert_condition(VoronoiCell::get_face_surface_area(faces[0]) == 1.);
    assert_condition(VoronoiCell::get_face_midpoint(faces[0]) ==
                     CoordinateVector<>(0.5, 0., 0.5));

    assert_condition(VoronoiCell::get_face_neighbour(faces[1]) ==
                     VORONOI_BOX_LEFT);
    assert_condition(VoronoiCell::get_face_surface_area(faces[1]) == 1.);
    assert_condition(VoronoiCell::get_face_midpoint(faces[1]) ==
                     CoordinateVector<>(0., 0.5, 0.5));

    assert_condition(VoronoiCell::get_face_neighbour(faces[2]) ==
                     VORONOI_BOX_BOTTOM);
    assert_condition(VoronoiCell::get_face_surface_area(faces[2]) == 1.);
    assert_condition(VoronoiCell::get_face_midpoint(faces[2]) ==
                     CoordinateVector<>(0.5, 0.5, 0.));

    assert_condition(VoronoiCell::get_face_neighbour(faces[3]) ==
                     VORONOI_BOX_TOP);
    assert_condition(VoronoiCell::get_face_surface_area(faces[3]) == 1.);
    assert_condition(VoronoiCell::get_face_midpoint(faces[3]) ==
                     CoordinateVector<>(0.5, 0.5, 1.));

    assert_condition(VoronoiCell::get_face_neighbour(faces[4]) ==
                     VORONOI_BOX_BACK);
    assert_condition(VoronoiCell::get_face_surface_area(faces[4]) == 1.);
    assert_condition(VoronoiCell::get_face_midpoint(faces[4]) ==
                     CoordinateVector<>(0.5, 1., 0.5));

    assert_condition(VoronoiCell::get_face_neighbour(faces[5]) ==
                     VORONOI_BOX_RIGHT);
    assert_condition(VoronoiCell::get_face_surface_area(faces[5]) == 1.);
    assert_condition(VoronoiCell::get_face_midpoint(faces[5]) ==
                     CoordinateVector<>(1., 0.5, 0.5));
  }

  /// test Voronoi cell intersection algorithm

  /// First part of the intersection algorithm: see anonymous enum above for
  /// a detailed description of each test.

  /// VORONOITEST_PATH_1_0
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_0);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.0: %i (%i %i %i %i)", status, varcheck[0], varcheck[1],
                 varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 0);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 1);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_1
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_1);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.1: %i (%i %i %i %i)", status, varcheck[0], varcheck[1],
                 varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 0);
    assert_condition(varcheck[1] == 1);
    assert_condition(varcheck[2] == 2);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_2
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_2);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.2: %i (%i %i %i %i)", status, varcheck[0], varcheck[1],
                 varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 0);
    assert_condition(varcheck[1] == 1);
    assert_condition(varcheck[2] == 2);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_3
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_3);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.3: %i", status);
    assert_condition(status == -1);
  }

  /// VORONOITEST_PATH_1_4_0
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_4_0);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.4.0: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 1);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 2);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_4_1
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_4_1);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.4.1: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 1);
    assert_condition(varcheck[1] == 1);
    assert_condition(varcheck[2] == 3);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_4_2
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_4_2);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.4.2: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 1);
    assert_condition(varcheck[1] == 1);
    assert_condition(varcheck[2] == 3);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_4_3
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_4_3);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.4.3: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 1);
    assert_condition(varcheck[1] == 1);
    assert_condition(varcheck[2] == 2);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_4_4
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_4_4);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.4.4: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 1);
    assert_condition(varcheck[1] == 2);
    assert_condition(varcheck[2] == 3);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_4_5
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_4_5);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.4.5: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 1);
    assert_condition(varcheck[1] == 2);
    assert_condition(varcheck[2] == 3);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_1_4_6
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_4_6);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.4.6: %i", status);
    assert_condition(status == -1);
  }

  /// VORONOITEST_PATH_1_5
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_1_5);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 1.5: %i (%i)", status, varcheck[0]);
    assert_condition(status == 2);
    assert_condition(varcheck[0] == 1);
  }

  /// VORONOITEST_PATH_2_0
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_0);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.0: %i (%i %i %i %i)", status, varcheck[0], varcheck[1],
                 varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 1);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 0);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_2_1
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_1);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.1: %i (%i %i %i %i)", status, varcheck[0], varcheck[1],
                 varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 2);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 0);
    assert_condition(varcheck[3] == 1);
  }

  /// VORONOITEST_PATH_2_2
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_2);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.2: %i (%i %i %i %i)", status, varcheck[0], varcheck[1],
                 varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 2);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 0);
    assert_condition(varcheck[3] == 1);
  }

  /// VORONOITEST_PATH_2_3
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_3);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.3: %i", status);
    assert_condition(status == 0);
  }

  /// VORONOITEST_PATH_2_4_0
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_4_0);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.4.0: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 2);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 1);
    assert_condition(varcheck[3] == 0);
  }

  /// VORONOITEST_PATH_2_4_1
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_4_1);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.4.1: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 3);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 1);
    assert_condition(varcheck[3] == 1);
  }

  /// VORONOITEST_PATH_2_4_2
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_4_2);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.4.2: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 3);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 1);
    assert_condition(varcheck[3] == 1);
  }

  /// VORONOITEST_PATH_2_4_3
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_4_3);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.4.3: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 2);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 1);
    assert_condition(varcheck[3] == 1);
  }

  /// VORONOITEST_PATH_2_4_4
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_4_4);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.4.4: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 3);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 1);
    assert_condition(varcheck[3] == 2);
  }

  /// VORONOITEST_PATH_2_4_5
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_4_5);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.4.5: %i (%i %i %i %i)", status, varcheck[0],
                 varcheck[1], varcheck[2], varcheck[3]);
    assert_condition(status == 1);
    assert_condition(varcheck[0] == 3);
    assert_condition(varcheck[1] == 0);
    assert_condition(varcheck[2] == 1);
    assert_condition(varcheck[3] == 2);
  }

  /// VORONOITEST_PATH_2_4_6
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_4_6);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.4.6: %i", status);
    assert_condition(status == 0);
  }

  /// VORONOITEST_PATH_2_5
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_PATH_2_5);
    CoordinateVector<> dx(1., 0., 0.);
    int varcheck[4];
    int status = cell.intersect(dx, 0, varcheck);
    cmac_warning("Path 2.5: %i (%i)", status, varcheck[0]);
    assert_condition(status == 2);
    assert_condition(varcheck[0] == 1);
  }

  /// VORONOITEST_VERTEX_DELETION
  {
    VoronoiCell cell;
    cell.setup_variables_for_test(VORONOITEST_VERTEX_DELETION);
    std::vector< bool > delete_stack(6, false);
    delete_stack[0] = true;
    cell.delete_vertices(delete_stack);
    cell.check_variables_after_test(VORONOITEST_VERTEX_DELETION);
  }

  /// VORONOITEST_INTERSECTION
  {
    VoronoiCell cell(CoordinateVector<>(0.5, 0.5, 0.5),
                     Box(CoordinateVector<>(0.), CoordinateVector<>(1.)));
    cell.check_variables_after_test(VORONOITEST_CONSTRUCTOR);
    int status = cell.intersect(CoordinateVector<>(0.5, -0.5, 0.), 0);
    assert_condition(status == 1);
    cell.check_variables_after_test(VORONOITEST_INTERSECTION);
  }

  /// test Voronoi grid
  {
    unsigned int numcell = 1000;
    VoronoiGrid grid(Box(CoordinateVector<>(0.), CoordinateVector<>(1.)), false,
                     numcell);
    for (unsigned int i = 0; i < numcell; ++i) {
      unsigned int new_index = grid.add_cell(Utilities::random_position());
      assert_condition(new_index == i);
    }
    grid.compute_grid();

    // print the grid while we have the grid connections
    std::ofstream file("test_voronoi_grid.txt");
    grid.print_grid(file);
    file.close();

    grid.finalize();
    const double tolerance = 1.e-11;
    double total_volume = 0.;
    for (unsigned int i = 0; i < numcell; ++i) {
      double cell_volume = grid.get_volume(i);
      // no single cell should fill the entire volume
      assert_condition(cell_volume < 1.);
      total_volume += cell_volume;

      // check that all neighbours of this cell have this cell as neighbour as
      // well
      const auto faces = grid.get_faces(i);
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        const unsigned int ngb = VoronoiCell::get_face_neighbour(*it);
        // some faces have the walls of the box as neighbour, we ignore these
        if (ngb < VORONOI_MAX_INDEX) {
          const double area = VoronoiCell::get_face_surface_area(*it);
          const CoordinateVector<> midpoint =
              VoronoiCell::get_face_midpoint(*it);
          const auto ngbfaces = grid.get_faces(ngb);
          auto ngbit = ngbfaces.begin();
          while (ngbit != ngbfaces.end() &&
                 VoronoiCell::get_face_neighbour(*ngbit) != i) {
            ++ngbit;
          }
          assert_condition(ngbit != ngbfaces.end());
          const double ngbarea = VoronoiCell::get_face_surface_area(*ngbit);
          const CoordinateVector<> ngbmidpoint =
              VoronoiCell::get_face_midpoint(*ngbit);
          assert_values_equal_rel(ngbarea, area, tolerance);
          assert_values_equal_rel(ngbmidpoint.x(), midpoint.x(), tolerance);
          assert_values_equal_rel(ngbmidpoint.y(), midpoint.y(), tolerance);
          assert_values_equal_rel(ngbmidpoint.z(), midpoint.z(), tolerance);
        }
      }
    }
    assert_values_equal_rel(total_volume, 1., 1.e-15);
  }

  /// test regular (degenerate) Voronoi grid
  {
    unsigned int numcell = 0;
    VoronoiGrid grid(Box(CoordinateVector<>(0.), CoordinateVector<>(1.)), false,
                     1000);
    for (unsigned int ix = 0; ix < 10; ++ix) {
      for (unsigned int iy = 0; iy < 10; ++iy) {
        for (unsigned int iz = 0; iz < 10; ++iz) {
          unsigned int new_index = grid.add_cell(CoordinateVector<>(
              (ix + 0.5) * 0.1, (iy + 0.5) * 0.1, (iz + 0.5) * 0.1));
          assert_condition(new_index == numcell);
          ++numcell;
        }
      }
    }
    grid.compute_grid();

    // print the grid while we have the grid connections
    std::ofstream file("test_voronoi_grid_regular.txt");
    grid.print_grid(file);
    file.close();

    grid.finalize();
    const double tolerance = 1.e-11;
    double total_volume = 0.;
    for (unsigned int i = 0; i < numcell; ++i) {
      double cell_volume = grid.get_volume(i);
      // no single cell should fill the entire volume
      assert_condition(cell_volume < 1.);
      total_volume += cell_volume;

      // check that all neighbours of this cell have this cell as neighbour as
      // well
      const auto faces = grid.get_faces(i);
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        const unsigned int ngb = VoronoiCell::get_face_neighbour(*it);
        // some faces have the walls of the box as neighbour, we ignore these
        if (ngb < VORONOI_MAX_INDEX) {
          const double area = VoronoiCell::get_face_surface_area(*it);
          const CoordinateVector<> midpoint =
              VoronoiCell::get_face_midpoint(*it);
          const auto ngbfaces = grid.get_faces(ngb);
          auto ngbit = ngbfaces.begin();
          while (ngbit != ngbfaces.end() &&
                 VoronoiCell::get_face_neighbour(*ngbit) != i) {
            ++ngbit;
          }
          assert_condition(ngbit != ngbfaces.end());
          const double ngbarea = VoronoiCell::get_face_surface_area(*ngbit);
          const CoordinateVector<> ngbmidpoint =
              VoronoiCell::get_face_midpoint(*ngbit);
          assert_values_equal_rel(ngbarea, area, tolerance);
          assert_values_equal_rel(ngbmidpoint.x(), midpoint.x(), tolerance);
          assert_values_equal_rel(ngbmidpoint.y(), midpoint.y(), tolerance);
          assert_values_equal_rel(ngbmidpoint.z(), midpoint.z(), tolerance);
        }
      }
    }
    assert_values_equal_rel(total_volume, 1., 1.e-15);
  }

  return 0;
}
