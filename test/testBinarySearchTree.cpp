/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testBinarySearchTree.cpp
 *
 * @brief Unit test for the BinarySearchTree class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "BinarySearchTree.hpp"
#include "Utilities.hpp"
#include <iostream>
#include <vector>

/**
 * @brief Unit test for the BinarySearchTree class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  BinarySearchTree< unsigned int, unsigned int > btree;

  const unsigned int numpoint = 20;
  std::cout << "Array:\n";
  std::vector< unsigned int > array(numpoint);
  for (unsigned int i = 0; i < numpoint; ++i) {
    array[i] = Utilities::random_int(0, 200);
    std::cout << i << ": " << array[i] << "\n";
    btree.add_point(array[i], i);
  }

  std::cout << "Tree:\n";
  btree.print(std::cout);
  std::cout << std::endl;

  assert_condition(btree.get_max_depth() == 5);

  std::vector< unsigned int > range = btree.get_range(102, 143);
  assert_condition(range.size() == 6);
  assert_condition(range[0] == 13);
  assert_condition(range[1] == 9);
  assert_condition(range[2] == 19);
  assert_condition(range[3] == 11);
  assert_condition(range[4] == 16);
  assert_condition(range[5] == 17);

  return 0;
}
