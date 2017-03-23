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
#include "VoronoiCell.hpp"
#include "VoronoiGrid.hpp"

/**
 * @brief Aliases for testcase setups.
 */
enum { VORONOITEST_PATH_1_0 = 0, VORONOITEST_PATH_1_1, VORONOITEST_PATH_1_2 };

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
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][0]) = 1;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[0][0]) = 0;
    _edges[1].resize(3);
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[1][0]) = 0;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[1][0]) = 0;
    break;
  }
  case VORONOITEST_PATH_1_1: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(2., 0., 0.);
    _vertices[2] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(3);
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][0]) = 1;
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][1]) = 2;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[0][1]) = 0;
    _edges[1].resize(3);
    _edges[2].resize(3);
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[2][0]) = 0;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[2][0]) = 1;
    break;
  }
  case VORONOITEST_PATH_1_2: {
    _vertices.resize(3);
    _vertices[0] = CoordinateVector<>(1., 0., 0.);
    _vertices[1] = CoordinateVector<>(2., 0., 0.);
    _vertices[2] = CoordinateVector<>(-1., 0., 0.);
    _edges.resize(3);
    _edges[0].resize(2);
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][0]) = 1;
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][1]) = 2;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[0][1]) = 0;
    _edges[1].resize(3);
    _edges[2].resize(3);
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[2][1]) = 0;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[2][0]) = 1;
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
    cell.finalize();
    assert_values_equal_rel(cell.get_volume(), 1., 1.e-16);
    assert_condition(cell.get_centroid() == CoordinateVector<>(0.5));

    auto faces = cell.get_faces();
    assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(
                         faces[0]) == VORONOI_BOX_FRONT);
    assert_condition(
        std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[0]) == 1.);
    assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[0]) ==
                     CoordinateVector<>(0.5, 0., 0.5));

    assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(
                         faces[1]) == VORONOI_BOX_LEFT);
    assert_condition(
        std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[1]) == 1.);
    assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[1]) ==
                     CoordinateVector<>(0., 0.5, 0.5));

    assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(
                         faces[2]) == VORONOI_BOX_BOTTOM);
    assert_condition(
        std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[2]) == 1.);
    assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[2]) ==
                     CoordinateVector<>(0.5, 0.5, 0.));

    assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(
                         faces[3]) == VORONOI_BOX_TOP);
    assert_condition(
        std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[3]) == 1.);
    assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[3]) ==
                     CoordinateVector<>(0.5, 0.5, 1.));

    assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(
                         faces[4]) == VORONOI_BOX_BACK);
    assert_condition(
        std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[4]) == 1.);
    assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[4]) ==
                     CoordinateVector<>(0.5, 1., 0.5));

    assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(
                         faces[5]) == VORONOI_BOX_RIGHT);
    assert_condition(
        std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[5]) == 1.);
    assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[5]) ==
                     CoordinateVector<>(1., 0.5, 0.5));
  }

  /// test Voronoi cell intersection algorithm
  {
    /// VORONOITEST_PATH_1_0
    {
      VoronoiCell cell;
      cell.setup_variables_for_test(VORONOITEST_PATH_1_0);
      CoordinateVector<> dx(1., 0., 0.);
      int status = cell.intersect(dx, 0, true);
      cmac_warning("Path 1.0: %i", status);
      assert_condition(status == 1);
    }

    /// VORONOITEST_PATH_1_1
    {
      VoronoiCell cell;
      cell.setup_variables_for_test(VORONOITEST_PATH_1_1);
      CoordinateVector<> dx(1., 0., 0.);
      int status = cell.intersect(dx, 0, true);
      cmac_warning("Path 1.1: %i", status);
      assert_condition(status == 1);
    }

    /// VORONOITEST_PATH_1_2
    {
      VoronoiCell cell;
      cell.setup_variables_for_test(VORONOITEST_PATH_1_2);
      CoordinateVector<> dx(1., 0., 0.);
      int status = cell.intersect(dx, 0, true);
      cmac_warning("Path 1.2: %i", status);
      assert_condition(status == 1);
    }
  }

  return 0;
}
