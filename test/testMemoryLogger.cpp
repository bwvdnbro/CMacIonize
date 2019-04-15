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
 * @file testMemoryLogger.cpp
 *
 * @brief Unit test for the MemoryLogger class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "MemoryLogger.hpp"

#include <iostream>

/**
 * @brief Unit test for the MemoryLogger class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  MemoryLogger memory_logger;

  memory_logger.add_entry("start");

  memory_logger.add_entry("block");
  double *data = new double[100000];
  memory_logger.finalize_entry();

  memory_logger.add_entry("peak");

  delete[] data;

  memory_logger.add_entry("stop");

  std::cout << "Snapshots:\n";
  memory_logger.print(std::cout, true);
  std::cout << "Allocations:\n";
  memory_logger.print(std::cout, false);

  return 0;
}
