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
 * @file testCoordinateVector.cpp
 *
 * @brief Unit test for the CoordinateVector class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CoordinateVector.hpp"
#include <cmath>

/**
 * @brief Unit test for the CoordinateVector class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // test constructors
  {
    // empty constructor
    CoordinateVector<> a;
    assert_condition(a.x() == 0.);
    assert_condition(a.y() == 0.);
    assert_condition(a.z() == 0.);

    // single value constructor
    CoordinateVector<> b(1.);
    assert_condition(b.x() == 1.);
    assert_condition(b.y() == 1.);
    assert_condition(b.z() == 1.);

    // normal constructor
    CoordinateVector<> c(1., 2., 3.);
    assert_condition(c.x() == 1.);
    assert_condition(c.y() == 2.);
    assert_condition(c.z() == 3.);
  }

  // test subtraction
  {
    CoordinateVector<> a(2., 3., 4.);
    CoordinateVector<> b(1., 2., 3.);
    a -= b;
    assert_condition(a.x() == 1.);
    assert_condition(a.y() == 1.);
    assert_condition(a.z() == 1.);

    CoordinateVector<> c(2., 3., 4.);
    CoordinateVector<> d = c - b;
    assert_condition(d.x() == 1.);
    assert_condition(d.y() == 1.);
    assert_condition(d.z() == 1.);
  }

  // test addition
  {
    CoordinateVector<> a(1., 2., 3.);
    CoordinateVector<> b(2., 3., 4.);
    a += b;
    assert_condition(a.x() == 3.);
    assert_condition(a.y() == 5.);
    assert_condition(a.z() == 7.);

    CoordinateVector<> c = a + b;
    assert_condition(c.x() == 5.);
    assert_condition(c.y() == 8.);
    assert_condition(c.z() == 11.);
  }

  // test multiplication
  {
    CoordinateVector<> a(1., 2., 3.);
    a *= 2.;
    assert_condition(a.x() == 2.);
    assert_condition(a.y() == 4.);
    assert_condition(a.z() == 6.);

    CoordinateVector<> b = 3. * a;
    assert_condition(b.x() == 6.);
    assert_condition(b.y() == 12.);
    assert_condition(b.z() == 18.);

    CoordinateVector<> c = a * 2.;
    assert_condition(c.x() == 4.);
    assert_condition(c.y() == 8.);
    assert_condition(c.z() == 12.);
  }

  // test division
  {
    CoordinateVector<> a(1., 2., 3.);
    a /= 2.;
    assert_condition(a.x() == 0.5);
    assert_condition(a.y() == 1.);
    assert_condition(a.z() == 1.5);

    CoordinateVector<> b = a / 2.;
    assert_condition(b.x() == 0.25);
    assert_condition(b.y() == 0.5);
    assert_condition(b.z() == 0.75);
  }

  // test norm
  {
    CoordinateVector<> a(1.);
    assert_condition(a.norm2() == 3.);
    assert_condition(a.norm() == sqrt(3.));
  }

  return 0;
}
