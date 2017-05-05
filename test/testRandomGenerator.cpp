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
 * @file testRandomGenerator.cpp
 *
 * @brief Unit test for the RandomGenerator.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Unit test for the RandomGenerator.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  /// Basic test
  {
    RandomGenerator generator(42);

    double mean_random = 0.;
    unsigned int num = 1000000;
    double weight = 1. / num;
    for (unsigned int i = 0; i < num; ++i) {
      mean_random += weight * generator.get_uniform_random_double();
    }
    assert_values_equal_tol(mean_random, 0.5, 1.e-3);
  }

  /// Seed test: check that different seeds effectively lead to different
  /// results
  {
    RandomGenerator generator_A(42);
    RandomGenerator generator_B(512);

    assert_condition(generator_A.get_uniform_random_double() !=
                     generator_B.get_uniform_random_double());
  }

  return 0;
}
