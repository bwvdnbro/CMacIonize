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
#include "CoordinateVector.hpp"
#include <cassert>
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
    CoordinateVector a;
    assert(a.x() == 0.);
    assert(a.y() == 0.);
    assert(a.z() == 0.);

    // single value constructor
    CoordinateVector b(1.);
    assert(b.x() == 1.);
    assert(b.y() == 1.);
    assert(b.z() == 1.);

    // normal constructor
    CoordinateVector c(1., 2., 3.);
    assert(c.x() == 1.);
    assert(c.y() == 2.);
    assert(c.z() == 3.);
  }

  // test subtraction
  {
    CoordinateVector a(2., 3., 4.);
    CoordinateVector b(1., 2., 3.);
    a -= b;
    assert(a.x() == 1.);
    assert(a.y() == 1.);
    assert(a.z() == 1.);

    CoordinateVector c(2., 3., 4.);
    CoordinateVector d = c - b;
    assert(d.x() == 1.);
    assert(d.y() == 1.);
    assert(d.z() == 1.);
  }

  // test norm
  {
    CoordinateVector a(1.);
    assert(a.norm2() == 3.);
    assert(a.norm() == sqrt(3.));
  }

  return 0;
}
