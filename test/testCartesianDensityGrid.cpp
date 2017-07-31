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
 * @file testCartesianDensityGrid.cpp
 *
 * @brief Unit test for the CartesianDensityGrid class
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "DensityFunction.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "Photon.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the CartesianDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit value: 0 on success.
 */
int main(int argc, char **argv) {
  HomogeneousDensityFunction testfunction(1., 2000.);
  CoordinateVector<> anchor;
  CoordinateVector<> sides(1., 1., 1.);
  Box<> box(anchor, sides);
  // this actually works! (note: it should not work, since 64 is not a
  // CoordinateVector<unsigned char>)
  // the reason it works is that we have defined a converting constructor
  // CoordinateVector<unsigned char>(unsigned char), which converts a single
  // unsigned char into a CoordinateVector<unsigned char>. The compiler is
  // smart enough to notice this, and automatically converts 64 to the required
  // CoordinateVector<unsigned char> argument.
  CartesianDensityGrid grid(box, 64, testfunction);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  assert_values_equal(1., grid.get_total_hydrogen_number());
  assert_values_equal(2000., grid.get_average_temperature());

  CoordinateVector<> photon_origin(0.51, 0.51, 0.51);
  unsigned long index = grid.get_cell_index(photon_origin);

  assert_condition(index == 32 * 64 * 64 + 32 * 64 + 32);

  Box<> cell = grid.get_cell(index);

  assert_condition(cell.get_anchor().x() == 0.5);
  assert_condition(cell.get_anchor().y() == 0.5);
  assert_condition(cell.get_anchor().z() == 0.5);
  assert_condition(cell.get_sides().x() == 0.015625);
  assert_condition(cell.get_sides().y() == 0.015625);
  assert_condition(cell.get_sides().z() == 0.015625);

  CoordinateVector<> cell_top_anchor = cell.get_top_anchor();

  assert_condition(cell_top_anchor.x() == 0.515625);
  assert_condition(cell_top_anchor.y() == 0.515625);
  assert_condition(cell_top_anchor.z() == 0.515625);

  auto ngbs = grid.get_neighbours(index);
  assert_condition(ngbs.size() == 6);
  unsigned long ngbindexexp =
      DensityGrid::iterator(
          grid.get_cell_index(CoordinateVector<>(0.49, 0.51, 0.51)), grid)
          .get_index();
  unsigned long ngbindex = std::get< 0 >(ngbs[0]).get_index();
  assert_condition(ngbindex == ngbindexexp);
  CoordinateVector<> midpoint = std::get< 1 >(ngbs[0]);
  assert_condition(midpoint.x() == 0.5);
  assert_condition(midpoint.y() == 0.5078125);
  assert_condition(midpoint.z() == 0.5078125);
  CoordinateVector<> normal = std::get< 2 >(ngbs[0]);
  assert_condition(normal.x() == -1.);
  assert_condition(normal.y() == 0.);
  assert_condition(normal.z() == 0.);
  double surface_area = std::get< 3 >(ngbs[0]);
  assert_condition(surface_area == 0.000244140625);

  ngbindexexp =
      DensityGrid::iterator(
          grid.get_cell_index(CoordinateVector<>(0.52, 0.51, 0.51)), grid)
          .get_index();
  ngbindex = std::get< 0 >(ngbs[1]).get_index();
  assert_condition(ngbindex == ngbindexexp);
  midpoint = std::get< 1 >(ngbs[1]);
  assert_condition(midpoint.x() == 0.515625);
  assert_condition(midpoint.y() == 0.5078125);
  assert_condition(midpoint.z() == 0.5078125);
  normal = std::get< 2 >(ngbs[1]);
  assert_condition(normal.x() == 1.);
  assert_condition(normal.y() == 0.);
  assert_condition(normal.z() == 0.);
  surface_area = std::get< 3 >(ngbs[1]);
  assert_condition(surface_area == 0.000244140625);

  ngbindexexp =
      DensityGrid::iterator(
          grid.get_cell_index(CoordinateVector<>(0.51, 0.49, 0.51)), grid)
          .get_index();
  ngbindex = std::get< 0 >(ngbs[2]).get_index();
  assert_condition(ngbindex == ngbindexexp);
  midpoint = std::get< 1 >(ngbs[2]);
  assert_condition(midpoint.x() == 0.5078125);
  assert_condition(midpoint.y() == 0.5);
  assert_condition(midpoint.z() == 0.5078125);
  normal = std::get< 2 >(ngbs[2]);
  assert_condition(normal.x() == 0.);
  assert_condition(normal.y() == -1.);
  assert_condition(normal.z() == 0.);
  surface_area = std::get< 3 >(ngbs[2]);
  assert_condition(surface_area == 0.000244140625);

  ngbindexexp =
      DensityGrid::iterator(
          grid.get_cell_index(CoordinateVector<>(0.51, 0.52, 0.51)), grid)
          .get_index();
  ngbindex = std::get< 0 >(ngbs[3]).get_index();
  assert_condition(ngbindex == ngbindexexp);
  midpoint = std::get< 1 >(ngbs[3]);
  assert_condition(midpoint.x() == 0.5078125);
  assert_condition(midpoint.y() == 0.515625);
  assert_condition(midpoint.z() == 0.5078125);
  normal = std::get< 2 >(ngbs[3]);
  assert_condition(normal.x() == 0.);
  assert_condition(normal.y() == 1.);
  assert_condition(normal.z() == 0.);
  surface_area = std::get< 3 >(ngbs[3]);
  assert_condition(surface_area == 0.000244140625);

  ngbindexexp =
      DensityGrid::iterator(
          grid.get_cell_index(CoordinateVector<>(0.51, 0.51, 0.49)), grid)
          .get_index();
  ngbindex = std::get< 0 >(ngbs[4]).get_index();
  assert_condition(ngbindex == ngbindexexp);
  midpoint = std::get< 1 >(ngbs[4]);
  assert_condition(midpoint.x() == 0.5078125);
  assert_condition(midpoint.y() == 0.5078125);
  assert_condition(midpoint.z() == 0.5);
  normal = std::get< 2 >(ngbs[4]);
  assert_condition(normal.x() == 0.);
  assert_condition(normal.y() == 0.);
  assert_condition(normal.z() == -1.);
  surface_area = std::get< 3 >(ngbs[4]);
  assert_condition(surface_area == 0.000244140625);

  ngbindexexp =
      DensityGrid::iterator(
          grid.get_cell_index(CoordinateVector<>(0.51, 0.51, 0.52)), grid)
          .get_index();
  ngbindex = std::get< 0 >(ngbs[5]).get_index();
  assert_condition(ngbindex == ngbindexexp);
  midpoint = std::get< 1 >(ngbs[5]);
  assert_condition(midpoint.x() == 0.5078125);
  assert_condition(midpoint.y() == 0.5078125);
  assert_condition(midpoint.z() == 0.515625);
  normal = std::get< 2 >(ngbs[5]);
  assert_condition(normal.x() == 0.);
  assert_condition(normal.y() == 0.);
  assert_condition(normal.z() == 1.);
  surface_area = std::get< 3 >(ngbs[5]);
  assert_condition(surface_area == 0.000244140625);

  // check faces
  {
    std::vector< Face > faces = grid.get_faces(index);

    assert_condition(faces[0].get_midpoint() ==
                     CoordinateVector<>(0.5, 0.5078125, 0.5078125));
    {
      auto vertices = faces[0].first_vertex();
      assert_condition(vertices != faces[0].last_vertex());
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.5, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.515625, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.515625, 0.515625));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.5, 0.515625));
      ++vertices;
      assert_condition(vertices == faces[0].last_vertex());
    }

    assert_condition(faces[1].get_midpoint() ==
                     CoordinateVector<>(0.515625, 0.5078125, 0.5078125));
    {
      auto vertices = faces[1].first_vertex();
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.5, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.515625, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.515625, 0.515625));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.5, 0.515625));
      ++vertices;
    }

    assert_condition(faces[2].get_midpoint() ==
                     CoordinateVector<>(0.5078125, 0.5, 0.5078125));
    {
      auto vertices = faces[2].first_vertex();
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.5, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.5, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.5, 0.515625));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.5, 0.515625));
      ++vertices;
    }

    assert_condition(faces[3].get_midpoint() ==
                     CoordinateVector<>(0.5078125, 0.515625, 0.5078125));
    {
      auto vertices = faces[3].first_vertex();
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.515625, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.515625, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.515625, 0.515625));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.515625, 0.515625));
      ++vertices;
    }

    assert_condition(faces[4].get_midpoint() ==
                     CoordinateVector<>(0.5078125, 0.5078125, 0.5));
    {
      auto vertices = faces[4].first_vertex();
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.5, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.5, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.515625, 0.5));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.515625, 0.5));
      ++vertices;
    }

    assert_condition(faces[5].get_midpoint() ==
                     CoordinateVector<>(0.5078125, 0.5078125, 0.515625));
    {
      auto vertices = faces[5].first_vertex();
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.5, 0.515625));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.5, 0.515625));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.515625, 0.515625, 0.515625));
      ++vertices;
      assert_condition(vertices.get_position() ==
                       CoordinateVector<>(0.5, 0.515625, 0.515625));
      ++vertices;
    }
  }

  // check different scenarios for the wall intersection algorithm
  CoordinateVector< char > next_index;
  double ds;
  {
    // positive x direction
    CoordinateVector<> photon_direction(1., 0., 0.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == 1);
    assert_condition(next_index.y() == 0);
    assert_condition(next_index.z() == 0);
    assert_condition(intersection.x() == 0.515625);
    assert_condition(intersection.y() == 0.51);
    assert_condition(intersection.z() == 0.51);
    assert_values_equal(ds, 0.005625);
  }
  {
    // negative x direction
    CoordinateVector<> photon_direction(-1., 0., 0.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == -1);
    assert_condition(next_index.y() == 0);
    assert_condition(next_index.z() == 0);
    assert_condition(intersection.x() == 0.5);
    assert_condition(intersection.y() == 0.51);
    assert_condition(intersection.z() == 0.51);
    assert_values_equal(ds, 0.01);
  }
  {
    // positive y direction
    CoordinateVector<> photon_direction(0., 1., 0.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == 0);
    assert_condition(next_index.y() == 1);
    assert_condition(next_index.z() == 0);
    assert_condition(intersection.x() == 0.51);
    assert_condition(intersection.y() == 0.515625);
    assert_condition(intersection.z() == 0.51);
    assert_values_equal(ds, 0.005625);
  }
  {
    // negative y direction
    CoordinateVector<> photon_direction(0., -1., 0.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == 0);
    assert_condition(next_index.y() == -1);
    assert_condition(next_index.z() == 0);
    assert_condition(intersection.x() == 0.51);
    assert_condition(intersection.y() == 0.5);
    assert_condition(intersection.z() == 0.51);
    assert_values_equal(ds, 0.01);
  }
  {
    // positive z direction
    CoordinateVector<> photon_direction(0., 0., 1.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == 0);
    assert_condition(next_index.y() == 0);
    assert_condition(next_index.z() == 1);
    assert_condition(intersection.x() == 0.51);
    assert_condition(intersection.y() == 0.51);
    assert_condition(intersection.z() == 0.515625);
    assert_values_equal(ds, 0.005625);
  }
  {
    // negative z direction
    CoordinateVector<> photon_direction(0., 0., -1.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == 0);
    assert_condition(next_index.y() == 0);
    assert_condition(next_index.z() == -1);
    assert_condition(intersection.x() == 0.51);
    assert_condition(intersection.y() == 0.51);
    assert_condition(intersection.z() == 0.5);
    assert_values_equal(ds, 0.01);
  }
  {
    // general direction
    CoordinateVector<> photon_direction(1., 2., -3.);
    // funny thing: the direction vector does not need to be normalized:
    // the norm cancels out in the equation for the intersection point
    // and all l parameters will always be in units of the norm, so values
    // can be safely compared
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    // the solution:
    //  l_x = 0.005625
    //  l_y = 0.005625/2 = 0.0028125
    //  l_z = 0.01/3 = 0.033...
    // ==> y closest
    assert_condition(next_index.x() == 0);
    assert_condition(next_index.y() == 1);
    assert_condition(next_index.z() == 0);
    // \vec{i} = \vec{o} + l_y \vec{d}
    // i_x = 0.51 + 0.0028125
    // i_y = 0.51 + 0.005625
    // i_z = 0.51 - 0.0084375
    assert_condition(intersection.x() == 0.5128125);
    assert_condition(intersection.y() == 0.515625);
    assert_condition(intersection.z() == 0.5015625);
    // ds = \sqrt{0.0028125^2 + 0.005625^2 + 0.0084375^2}
    assert_values_equal(ds, 0.0105234114);
  }
  {
    // intersection point on edge of two walls
    CoordinateVector<> photon_direction(0., 1., 1.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == 0);
    assert_condition(next_index.y() == 1);
    assert_condition(next_index.z() == 1);
    assert_condition(intersection.x() == 0.51);
    assert_condition(intersection.y() == 0.515625);
    assert_condition(intersection.z() == 0.515625);
    assert_values_equal(ds, 0.007954951288);
  }
  {
    // intersection point on corner of cell
    CoordinateVector<> photon_direction(1., 1., 1.);
    CoordinateVector<> intersection = grid.get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);
    assert_condition(next_index.x() == 1);
    assert_condition(next_index.y() == 1);
    assert_condition(next_index.z() == 1);
    assert_condition(intersection.x() == 0.515625);
    assert_condition(intersection.y() == 0.515625);
    assert_condition(intersection.z() == 0.515625);
    assert_values_equal(ds, 0.009742785793);
  }

  CoordinateVector<> photon_direction(1., 0., 0.);
  Photon photon(photon_origin, photon_direction, 1.);
  photon.set_cross_section(ION_H_n, 1.);
  photon.set_cross_section(ION_He_n, 1.);
  DensityGrid::iterator inside = grid.interact(photon, 0.125);

  assert_condition(inside == grid.end());

  return 0;
}
