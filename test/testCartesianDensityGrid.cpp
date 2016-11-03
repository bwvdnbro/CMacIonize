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
using namespace std;

/**
 * @brief Unit test for the CartesianDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit value: 0 on success.
 */
int main(int argc, char **argv) {
  HomogeneousDensityFunction testfunction;
  CoordinateVector<> anchor;
  CoordinateVector<> sides(1., 1., 1.);
  Box box(anchor, sides);
  // this actually works! (note: it should not work, since 64 is not a
  // CoordinateVector<unsigned char>)
  // the reason it works is that we have defined a converting constructor
  // CoordinateVector<unsigned char>(unsigned char), which converts a single
  // unsigned char into a CoordinateVector<unsigned char>. The compiler is
  // smart enough to notice this, and automatically converts 64 to the required
  // CoordinateVector<unsigned char> argument.
  CartesianDensityGrid grid(box, 64, 0.1, 8000., testfunction);

  assert_values_equal(1., grid.get_total_hydrogen_number());

  CoordinateVector<> photon_origin(0.51, 0.51, 0.51);
  unsigned long index = grid.get_cell_index(photon_origin);

  assert_condition(index == 32 * 64 * 64 + 32 * 64 + 32);

  Box cell = grid.get_cell(index);

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
  Photon photon(photon_origin, photon_direction, 1., 1., 1.);
  bool inside = grid.interact(photon, 0.125);

  assert_condition(inside == false);

  return 0;
}
