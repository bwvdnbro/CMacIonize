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
 * @file testAtomicValue.cpp
 *
 * @brief Unit test for the AtomicValue class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "AtomicValue.hpp"
#include <cinttypes>
#include <omp.h>

/**
 * @brief Unit test for the AtomicValue class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // this test is a bit tricky, since it will only fail if two threads access
  // the sum variable at the same time, which does not always happen (it will
  // of course also only fail if Atomic is badly implemented)
  // to increase our chances of having a collision (two threads accessing sum at
  // the same time) we run the test a couple of times and we use a ridiculously
  // high number of threads
  omp_set_num_threads(512);
  for (uint_fast8_t i_loop = 0; i_loop < 100; ++i_loop) {
    AtomicValue< int_fast32_t > sum(0);
    AtomicValue< int_fast32_t > max(0);
    int_fast32_t nthread = 0;
#pragma omp parallel shared(sum, nthread)
    {
      const int_fast32_t this_thread = omp_get_thread_num();
      if (this_thread % 2 == 0) {
        sum.post_add(this_thread);
        sum.pre_increment();
      } else {
        sum.pre_add(this_thread);
        sum.post_increment();
      }
      max.max(this_thread);
#pragma omp single
      { nthread = omp_get_num_threads(); }
    }

    int_fast32_t reference = 0;
    for (int_fast32_t i = 0; i < nthread; ++i) {
      reference += i + 1;
    }

    assert_condition(sum.value() == reference);
    assert_condition(max.value() == nthread - 1);
  }

  return 0;
}
