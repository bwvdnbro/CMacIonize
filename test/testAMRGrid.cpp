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
 * @file testAMRGrid.cpp
 *
 * @brief Unit test for the AMRGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AMRGrid.hpp"
#include "Box.hpp"
#include "CoordinateVector.hpp"
#include <fstream>

/**
 * @brief Unit test for the AMRGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  AMRGrid grid(Box(CoordinateVector<>(), CoordinateVector<>(1.)),
               CoordinateVector< int >(1));

  grid.create_cell(2, CoordinateVector<>(0.1));

  std::ofstream gfile("amrgrid.dump");
  grid.print(gfile);

  return 0;
}
