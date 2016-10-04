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
 * @file testParameterFile.cpp
 *
 * @brief Unit test for the ParameterFile class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "ParameterFile.hpp"
#include <iostream>
using namespace std;

/**
 * @brief Unit test for the ParameterFile class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  ParameterFile params("test.param");

  params.print_contents(cout);

  assert_condition(params.get_value< int >("test_integer1") == 42);
  assert_condition(params.get_value< int >("test_integer2") == 42);
  assert_condition(params.get_value< int >("test_integer3") == 42);
  assert_condition(params.get_value< double >("test_float") == 3.14);
  assert_condition(params.get_value< bool >("test_bool1") == true);
  assert_condition(params.get_value< bool >("test_bool2") == true);
  assert_condition(params.get_value< bool >("test_bool3") == true);
  assert_condition(params.get_value< bool >("test_bool4") == true);
  assert_condition(params.get_value< bool >("test_bool5") == false);
  assert_condition(params.get_value< bool >("test_bool6") == false);
  assert_condition(params.get_value< bool >("test_bool7") == false);
  assert_condition(params.get_value< bool >("test_bool8") == false);
  assert_condition(params.get_value< string >("test_string") ==
                   "This is a test string.");
  assert_condition(params.get_value< int >("test_group.test_group_member") ==
                   42);
  assert_condition(
      params.get_value< string >("test_comments_group.test_comments_value") ==
      "test comments string");

  return 0;
}
