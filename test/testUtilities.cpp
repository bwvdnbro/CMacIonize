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
 * @file testUtilities.cpp
 *
 * @brief Unit test for the Utilities namespace.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Utilities.hpp"

/**
 * @brief Unit test for the Utilities namespace.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  double time = 232354353.1233;
  std::string timestr = Utilities::human_readable_time(time);

  assert_condition(timestr == "7y 134d 6h 52m 33.1233s");

  unsigned long bytes = 125;
  assert_condition(Utilities::human_readable_bytes(bytes) == "125 bytes");
  bytes = 1253626623;
  assert_condition(Utilities::human_readable_bytes(bytes) == "1.17 GB");

  return 0;
}
