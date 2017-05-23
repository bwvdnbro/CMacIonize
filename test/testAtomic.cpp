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
 * @file testAtomic.cpp
 *
 * @brief Unit test for the Atomic class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Atomic.hpp"
#include <omp.h>

/**
 * @brief Unit test for the Atomic class.
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
  // the same time) we run the test a couple of times
  for (unsigned int i = 0; i < 100; ++i) {
    int sum = 0;
    int nthread = 0;
#pragma omp parallel shared(sum, nthread)
    {
      Atomic::add(sum, omp_get_thread_num());
#pragma omp single
      { nthread = omp_get_num_threads(); }
    }

    if (nthread > 1) {
      int reference = 0;
      for (int i = 0; i < nthread; ++i) {
        reference += i;
      }

      assert_condition(sum == reference);
    } else {
      cmac_status("This test only works if OMP_NUM_THREADS > 1.");
    }
  }

#ifdef HAVE_ATOMIC
  for (unsigned int i = 0; i < 100; ++i) {
    double sum = 0.;
    int nthread = 0.;
#pragma omp parallel shared(sum, nthread)
    {
      double this_thread = omp_get_thread_num();
      Atomic::add(sum, this_thread);
#pragma omp single
      { nthread = omp_get_num_threads(); }
    }

    if (nthread > 1) {
      double reference = 0.;
      for (int i = 0; i < nthread; ++i) {
        reference += i;
      }

      assert_condition(sum == reference);
    } else {
      cmac_status("This test only works if OMP_NUM_THREADS > 1.");
    }
  }
#endif

  return 0;
}
