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
 * @file testParallelCartesianDensityGrid.cpp
 *
 * @brief Unit test for the ParallelCartesianDensityGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Box.hpp"
#include "ParallelCartesianDensityGrid.hpp"

/**
 * @brief Unit test for the ParallelCartesianDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // this is shorthand for box(CoordinateVector<>(0.), CoordinateVector<>(1.))
  // and only works because we have defined the single value constructor for
  // CoordinateVector<>
  Box box(0., 1.);
  CoordinateVector< int > numcell(32, 32, 32);
  unsigned int numdomain = 64;
  std::pair< int, int > domain = std::make_pair(0, 64);

  ParallelCartesianDensityGrid grid(box, numcell, numdomain, domain);

  return 0;
}
