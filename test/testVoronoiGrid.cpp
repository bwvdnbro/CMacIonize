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
 * @brief Unit test for the VoronoiGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
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

  Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  VoronoiCell cell(CoordinateVector<>(0.5), box);
  cell.finalize();
  assert_values_equal_rel(cell.get_volume(), 1., 1.e-16);
  assert_condition(cell.get_centroid() == CoordinateVector<>(0.5));

  auto faces = cell.get_faces();
  assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(faces[0]) ==
                   VORONOI_BOX_FRONT);
  assert_condition(
      std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[0]) == 1.);
  assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[0]) ==
                   CoordinateVector<>(0.5, 0., 0.5));

  assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(faces[1]) ==
                   VORONOI_BOX_LEFT);
  assert_condition(
      std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[1]) == 1.);
  assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[1]) ==
                   CoordinateVector<>(0., 0.5, 0.5));

  assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(faces[2]) ==
                   VORONOI_BOX_BOTTOM);
  assert_condition(
      std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[2]) == 1.);
  assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[2]) ==
                   CoordinateVector<>(0.5, 0.5, 0.));

  assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(faces[3]) ==
                   VORONOI_BOX_TOP);
  assert_condition(
      std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[3]) == 1.);
  assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[3]) ==
                   CoordinateVector<>(0.5, 0.5, 1.));

  assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(faces[4]) ==
                   VORONOI_BOX_BACK);
  assert_condition(
      std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[4]) == 1.);
  assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[4]) ==
                   CoordinateVector<>(0.5, 1., 0.5));

  assert_condition(std::get< VoronoiCell::VORONOI_FACE_NEIGHBOUR >(faces[5]) ==
                   VORONOI_BOX_RIGHT);
  assert_condition(
      std::get< VoronoiCell::VORONOI_FACE_SURFACE_AREA >(faces[5]) == 1.);
  assert_condition(std::get< VoronoiCell::VORONOI_FACE_MIDPOINT >(faces[5]) ==
                   CoordinateVector<>(1., 0.5, 0.5));

  return 0;
}
