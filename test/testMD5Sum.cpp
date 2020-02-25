/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testMD5Sum.cpp
 *
 * @brief Unit test for the MD5Sum class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "MD5Sum.hpp"

/**
 * @brief Check that the given input string hashes into the given expected MD5
 * checksum.
 *
 * @param input Input string.
 * @param expected_output Expected output MD5 checksum.
 */
#define check_string(input, expected_output)                                   \
  {                                                                            \
    std::string real_output = MD5Sum::get_checksum(input);                     \
    cmac_status("\"%s\" -> %s", input, real_output.c_str());                   \
    assert_condition(real_output == expected_output);                          \
  }

/**
 * @brief Check that the given input file hashes into the given expected MD5
 * checksum.
 *
 * @param input Input file name.
 * @param expected_output Expected output MD5 checksum.
 */
#define check_file(input, expected_output)                                     \
  {                                                                            \
    std::string real_output = MD5Sum::get_file_checksum(input);                \
    cmac_status("\"%s\" -> %s", input, real_output.c_str());                   \
    assert_condition(real_output == expected_output);                          \
  }

/**
 * @brief Unit test for the MD5Sum class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// string checks

  // the empty message
  check_string("", "d41d8cd98f00b204e9800998ecf8427e");
  // the second example message on Wikipedia
  check_string("The quick brown fox jumps over the lazy dog",
               "9e107d9d372bb6826bd81d3542a419d6");
  // a message that is exactly 64 bytes long
  check_string(
      "0123456789abcdefghijklmnopqrstuvwxyz0123456789abcdefghijklmnopq\n",
      "ba31d6f258843f6d42e0cf4055feac7e");
  // a message that is short enough for the extra 1 bit, but too long to add
  // additional length data
  check_string(
      "0123456789abcdefghijklmnopqrstuvwxyz0123456789abcdefghijklmno\n",
      "05c48a78500e550faf80b077446e30f2");

  /// file checks

  // text file
  check_file("h0_testdata.txt", "2d2e4b19e63b7f0821ea1e11c2d7d684");
  // binary file
  check_file("FLASHtest.hdf5", "0341a923c893996e64508c69db21730b");

  return 0;
}
