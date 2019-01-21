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
 * @file testThreadSafeVector.cpp
 *
 * @brief Unit test for the ThreadSafeVector class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "ThreadSafeVector.hpp"

#include <omp.h>

/**
 * @brief Unit test for the ThreadSafeVector class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // make sure we use way too many threads, to force collisions
  omp_set_num_threads(512);

  // repeat the exercise a couple of times to increase the change of collisions
  // even more
  for (uint_fast32_t iloop = 0; iloop < 100; ++iloop) {

    // first fill the vector in a parallel environment
    // every thread requests one element and sets it to its thread number
    ThreadSafeVector< int_fast32_t > vector(512);
#pragma omp parallel default(shared)
    {
      const int_fast32_t this_thread = omp_get_thread_num();
      const size_t index = vector.get_free_element();
      vector[index] = this_thread;
    }

    // now check that all values are present
    // we make a flag for each individual element that changes from false to
    // true if that element is present
    bool flags[512];
    for (uint_fast32_t i = 0; i < 512; ++i) {
      flags[i] = false;
    }
    for (size_t i = 0; i < 512; ++i) {
      flags[vector[i]] = true;
    }
    for (uint_fast32_t i = 0; i < 512; ++i) {
      assert_condition(flags[i]);
    }

    // now check that the safe element function returns the size of the vector,
    // meaning it is full
    assert_condition(vector.get_free_element_safe() == vector.max_size());
  }

  return 0;
}
