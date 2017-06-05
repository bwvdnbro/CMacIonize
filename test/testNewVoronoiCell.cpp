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
    NewVoronoiCell cell(0);
    cell.finalize(box, positions);

    double volume = cell.get_volume();
    assert_condition(volume == 1.);
    cmac_status("Cell volume computation works!");

    CoordinateVector<> centroid = cell.get_centroid();
    const double tolerance = 1.e-16;
    assert_values_equal_rel(centroid.x(), 0.5, tolerance);
    assert_values_equal_rel(centroid.y(), 0.5, tolerance);
    assert_values_equal_rel(centroid.z(), 0.5, tolerance);
    cmac_status("Cell centroid computation works!");

    std::vector< VoronoiFace > faces = cell.get_faces();

    assert_condition(faces[0].get_surface_area() == 1.);
    CoordinateVector<> midpoint = faces[0].get_midpoint();
    assert_condition(midpoint.x() == 0.);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[0].get_neighbour() == NEWVORONOICELL_BOX_LEFT);

    assert_condition(faces[1].get_surface_area() == 1.);
    midpoint = faces[1].get_midpoint();
    assert_condition(midpoint.x() == 1.);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[1].get_neighbour() == NEWVORONOICELL_BOX_RIGHT);

    assert_condition(faces[2].get_surface_area() == 1.);
    midpoint = faces[2].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_condition(midpoint.y() == 0.);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[2].get_neighbour() == NEWVORONOICELL_BOX_FRONT);

    assert_condition(faces[3].get_surface_area() == 1.);
    midpoint = faces[3].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_condition(midpoint.y() == 1.);
    assert_values_equal_rel(midpoint.z(), 0.5, tolerance);
    assert_condition(faces[3].get_neighbour() == NEWVORONOICELL_BOX_BACK);

    assert_condition(faces[4].get_surface_area() == 1.);
    midpoint = faces[4].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_condition(midpoint.z() == 0.);
    assert_condition(faces[4].get_neighbour() == NEWVORONOICELL_BOX_BOTTOM);

    assert_condition(faces[5].get_surface_area() == 1.);
    midpoint = faces[5].get_midpoint();
    assert_values_equal_rel(midpoint.x(), 0.5, tolerance);
    assert_values_equal_rel(midpoint.y(), 0.5, tolerance);
    assert_condition(midpoint.z() == 1.);
    assert_condition(faces[5].get_neighbour() == NEWVORONOICELL_BOX_TOP);

    cmac_status("Cell face computation works!");
  }

  /// test for NewVoronoiCell::find_tetrahedron
  {
    NewVoronoiCell cell(0);
    CoordinateVector< unsigned long > box_anchor(1000);
    CoordinateVector< unsigned long > box_sides(1000);
    std::vector< CoordinateVector< unsigned long > > positions(4);
    positions[0] = CoordinateVector< unsigned long >(1500);
    // general point
    positions[1] = CoordinateVector< unsigned long >(1250);
    // point on face
    positions[2] = CoordinateVector< unsigned long >(1500, 1250, 1250);
    // point on axis
    positions[3] = CoordinateVector< unsigned long >(1500, 1500, 1250);
    VoronoiBox< unsigned long > box(positions[0], box_anchor, box_sides);

    unsigned int tetrahedra[4];
    assert_condition(cell.find_tetrahedron(1, box, positions, tetrahedra) == 1);
    assert_condition(tetrahedra[0] == 7);

    assert_condition(cell.find_tetrahedron(2, box, positions, tetrahedra) == 2);
    assert_condition(tetrahedra[0] == 7);
    assert_condition(tetrahedra[1] == 6);

    assert_condition(cell.find_tetrahedron(3, box, positions, tetrahedra) == 4);
    assert_condition(tetrahedra[0] == 7);
    assert_condition(tetrahedra[1] == 4);
    assert_condition(tetrahedra[2] == 6);
    assert_condition(tetrahedra[3] == 5);

    cmac_status("Find tetrahedron works!");
  }

  return 0;
}
