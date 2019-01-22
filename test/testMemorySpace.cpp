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
 * @file testMemorySpace.cpp
 *
 * @brief Unit test for the MemorySpace class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "MemorySpace.hpp"

#include <omp.h>

/**
 * @brief Unit test for the MemorySpace class.
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

    // create the memory space
    MemorySpace space(1000);

    // make all threads request a free buffer and then release it again
#pragma omp parallel default(shared)
    {
      size_t buffer = space.get_free_buffer();
      assert_condition(space[buffer].size() == 0);
      space.free_buffer(buffer);
    }

    // now check that all buffers were correctly released by serially trying
    // to request all of them again
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      space.get_free_buffer();
    }
  }

  return 0;
}
