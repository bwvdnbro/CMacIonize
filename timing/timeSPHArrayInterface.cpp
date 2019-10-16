/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file timeSPHArrayInterface.cpp
 *
 * @brief Timing test for the SPHArrayInterface.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "TimingTools.hpp"

/**
 * @brief Timing test for the SPHArrayInterface.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  timingtools_init("timeSPHArrayInterface", argc, argv);

  // create text file with cell generator positions
  // create parameter file: generator position file, whether or not Lloyd,
  // type of mapping
  // arrays for library

  // call init
  // call library function

  return 0;
}
