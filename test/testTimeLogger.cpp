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
 * @file testTimeLogger.cpp
 *
 * @brief Unit test for the TimeLogger class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "TimeLogger.hpp"

#include <cmath>
#include <iostream>

/**
 * @brief Unit test for the TimeLogger class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  TimeLogger time_logger;

  double value = 0.;

  time_logger.start("first loop");
  time_logger.start("first sub loop");
  for (uint_fast32_t i = 0; i < 1e6; ++i) {
    const double x = (0.5 + i);
    value += std::cos(2. * M_PI * x);
  }
  time_logger.end("first sub loop");
  time_logger.start("second sub loop");
  for (uint_fast32_t i = 0; i < 1e6; ++i) {
    const double x = (0.5 + i) * 0.1;
    value += std::sin(2. * M_PI * x);
  }
  time_logger.end("second sub loop");
  time_logger.end("first loop");

  time_logger.start("second loop");
  for (uint_fast32_t i = 0; i < 1e6; ++i) {
    const double x = (0.5 + i) * 0.1;
    value += std::sin(2. * M_PI * x);
  }
  time_logger.end("second loop");

  std::cout << "Result: " << value << std::endl;

  time_logger.output("test_timelogger.txt");

  return 0;
}
