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
#include "CoordinateVector.hpp"
#include "ParameterFile.hpp"
#include "YAMLDictionary.hpp"
#include <iostream>

/**
 * @brief Unit test for the ParameterFile class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// YAMLDictionary test
  {
    std::string test_string = "test_value_1: 42.\n"
                              "test_value_2: 1. m";
    std::istringstream test_stream(test_string);
    YAMLDictionary dictionary(test_stream);

    dictionary.print_contents(std::cout);
  }

  ParameterFile params("test.param");

  assert_condition(params.get_value< int_fast32_t >("test_integer1") == 42);
  assert_condition(params.get_value< int_fast32_t >("test_integer2") == 42);
  assert_condition(params.get_value< int_fast32_t >("test_integer3") == 42);
  assert_condition(params.get_value< int_fast32_t >("test_integer4") == 1e6);
  assert_condition(params.get_value< int_fast32_t >("test_integer5") == 1e6);
  assert_condition(params.get_value< double >("test_float") == 3.14);
  assert_condition(params.get_physical_value< QUANTITY_LENGTH >("test_unit") ==
                   3.086e16);
  assert_condition(params.get_value< bool >("test_bool1") == true);
  assert_condition(params.get_value< bool >("test_bool2") == true);
  assert_condition(params.get_value< bool >("test_bool3") == true);
  assert_condition(params.get_value< bool >("test_bool4") == true);
  assert_condition(params.get_value< bool >("test_bool5") == false);
  assert_condition(params.get_value< bool >("test_bool6") == false);
  assert_condition(params.get_value< bool >("test_bool7") == false);
  assert_condition(params.get_value< bool >("test_bool8") == false);
  assert_condition(params.get_value< std::string >("test_string") ==
                   "This is a test string.");
  assert_condition(
      params.get_value< int_fast32_t >("test_group:test_group_member") == 42);
  assert_condition(params.get_value< std::string >(
                       "test_comments_group:test_comments_value") ==
                   "test comments string");

  CoordinateVector<> cvtest =
      params.get_value< CoordinateVector<> >("test_coordinatevector_double");
  assert_condition(cvtest.x() == 0.1);
  assert_condition(cvtest.y() == 0.2);
  assert_condition(cvtest.z() == 0.3);

  CoordinateVector< int_fast32_t > cvtest2 =
      params.get_value< CoordinateVector< int_fast32_t > >(
          "test_coordinatevector_int");
  assert_condition(cvtest2.x() == 42);
  assert_condition(cvtest2.y() == 42);
  assert_condition(cvtest2.z() == 40);

  CoordinateVector< bool > cvtest3 =
      params.get_value< CoordinateVector< bool > >(
          "test_coordinatevector_bool");
  assert_condition(cvtest3.x() == false);
  assert_condition(cvtest3.y() == true);
  assert_condition(cvtest3.z() == true);

  CoordinateVector<> cvtest_unit =
      params.get_physical_vector< QUANTITY_LENGTH >(
          "test_coordinatevector_unit");
  assert_condition(cvtest_unit.x() == 3.086e16);
  assert_condition(cvtest_unit.y() == 2.);
  assert_condition(cvtest_unit.z() == 2.4e19);

  assert_condition(
      params.get_value< int_fast32_t >(
          "test_group2:test_group_group:test_group_group_member") == 42);

  // default values
  assert_condition(params.get_value< int_fast32_t >("not_in_file1", 42) == 42);
  assert_condition(params.get_value< double >("not_in_file2", 3.14) == 3.14);
  assert_condition(params.get_physical_value< QUANTITY_LENGTH >(
                       "unit_not_in_file", "1. pc") == 3.086e16);
  assert_condition(params.get_value< std::string >("not_in", "file?") ==
                   "file?");
  assert_condition(params.get_value< bool >("not_in_file3", true) == true);

  cvtest = params.get_value< CoordinateVector<> >(
      "not_in_file4", CoordinateVector<>(0.1, 0.2, 0.3));
  assert_condition(cvtest.x() == 0.1);
  assert_condition(cvtest.y() == 0.2);
  assert_condition(cvtest.z() == 0.3);

  cvtest2 = params.get_value< CoordinateVector< int_fast32_t > >(
      "not_in_file5", CoordinateVector< int_fast32_t >(42, 42, 42));
  assert_condition(cvtest2.x() == 42);
  assert_condition(cvtest2.y() == 42);
  assert_condition(cvtest2.z() == 42);

  cvtest3 = params.get_value< CoordinateVector< bool > >(
      "not_in_file6", CoordinateVector< bool >(true, false, true));
  assert_condition(cvtest3.x() == true);
  assert_condition(cvtest3.y() == false);
  assert_condition(cvtest3.z() == true);

  cvtest = params.get_physical_vector< QUANTITY_LENGTH >(
      "coordinatevector_unit_not_in_file", "[1. pc,2.m,2.4e19m]");
  assert_condition(cvtest_unit.x() == 3.086e16);
  assert_condition(cvtest_unit.y() == 2.);
  assert_condition(cvtest_unit.z() == 2.4e19);

  assert_condition(params.get_value< std::string >(
                       "test_group2:test_group_group:test_group_str_member",
                       "hello!") == "hello!");

  params.print_contents(std::cout);

  return 0;
}
