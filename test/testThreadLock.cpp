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
 * @file testThreadLock.cpp
 *
 * @brief Unit test for the ThreadLock class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "ThreadLock.hpp"

#include <omp.h>

/**
 * @brief Unit test for the ThreadLock class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // make sure we use way too many threads, to force collisions
  omp_set_num_threads(512);

  // repeat the exercise a couple of times to increase the chance of collisions
  // even more
  for (uint_fast32_t iloop = 0; iloop < 100; ++iloop) {

    int_fast32_t sum = 0;
    ThreadLock lock;
#pragma omp parallel default(shared)
    {
      lock.lock();
      const int_fast32_t this_thread = omp_get_thread_num();
      sum += this_thread;
      lock.unlock();
    }

    int_fast32_t reference = 0;
    for (uint_fast32_t i = 0; i < 512; ++i) {
      reference += i;
    }
    assert_condition(sum == reference);
  }

  return 0;
}
