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

  std::string haystack;
  std::string needle = "dirty";
  haystack = "v0.1-dirty";
  assert_condition(Utilities::string_ends_with(haystack, needle) == true);
  haystack = "v0.1-eldkef";
  assert_condition(Utilities::string_ends_with(haystack, needle) == false);
  haystack = "v0.6";
  assert_condition(Utilities::string_ends_with(haystack, needle) == false);
  haystack = "v0.6.";
  assert_condition(Utilities::string_ends_with(haystack, needle) == false);
  haystack = "dirty";
  assert_condition(Utilities::string_ends_with(haystack, needle) == true);
  haystack = "dirty-adef22f";
  assert_condition(Utilities::string_ends_with(haystack, needle) == false);

  int number = 32768;
  std::vector< int > components = Utilities::decompose(number);
  int check = 1;
  for (unsigned int i = 0; i < components.size(); ++i) {
    check *= components[i];
  }
  assert_condition(check == number);

  number = 32769;
  components = Utilities::decompose(number);
  check = 1;
  for (unsigned int i = 0; i < components.size(); ++i) {
    check *= components[i];
  }
  assert_condition(check == number);

  number = 37;
  components = Utilities::decompose(number);
  check = 1;
  for (unsigned int i = 0; i < components.size(); ++i) {
    check *= components[i];
  }
  assert_condition(check == number);

  CoordinateVector< int > ncell(32);
  unsigned int numblock = 64;
  CoordinateVector< int > block = Utilities::subdivide(ncell, numblock);
  assert_condition(block.x() == 8);
  assert_condition(block.y() == 8);
  assert_condition(block.z() == 8);

  ncell = CoordinateVector< int >(33, 50, 100);
  numblock = 100;
  block = Utilities::subdivide(ncell, numblock);
  assert_condition(block.x() == 3);
  assert_condition(block.y() == 10);
  assert_condition(block.z() == 20);

  unsigned long test_value = 0xf300;
  assert_condition(Utilities::as_binary_sequence(test_value) ==
                   "0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 "
                   "0000 1111 0011 0000 0000")

      std::vector< unsigned int >
          test_vector(5);
  test_vector[0] = 3;
  test_vector[1] = 2;
  test_vector[2] = 4;
  test_vector[3] = 0;
  test_vector[4] = 1;
  std::vector< unsigned int > idx_test_vector = Utilities::argsort(test_vector);
  assert_condition(idx_test_vector[0] == 3);
  assert_condition(idx_test_vector[1] == 4);
  assert_condition(idx_test_vector[2] == 1);
  assert_condition(idx_test_vector[3] == 0);
  assert_condition(idx_test_vector[4] == 2);

  return 0;
}
