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
