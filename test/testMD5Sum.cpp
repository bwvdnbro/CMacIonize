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
 * @brief Unit test for the MD5Sum class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  //  assert_condition(MD5Sum::get_checksum("h0_testdata.txt") ==
  //  "2d2e4b19e63b7f0821ea1e11c2d7d684");

  cmac_status("Checksum: %s", MD5Sum::get_checksum("").c_str());

  assert_condition(MD5Sum::get_checksum("") ==
                   "d41d8cd98f00b204e9800998ecf8427e");

  return 0;
}
