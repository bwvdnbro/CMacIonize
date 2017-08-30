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
 * @file WorkEnvironment.hpp
 *
 * @brief Class that is responsible for managing the global number of threads.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef WORKENVIRONMENT_HPP
#define WORKENVIRONMENT_HPP

#include "Configuration.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#ifdef HAVE_OUTPUT_CYCLES
#include "Error.hpp"
#include "Timer.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#endif

/**
 * @brief Class that is responsible for managing the global number of threads.
 */
class WorkEnvironment {
public:
  /**
   * @brief Set the system maximum number of threads.
   *
   * This number is used by default by all WorkDistributor instances. A
   * WorkDistributor instance can use a lower number of threads.
   *
   * If set to a negative or zero value, the number of threads will be set to
   * the maximum number of threads available on the system.
   *
   * @param max_num_threads Maximum number of threads to use.
   * @return Number of threads that will be available for use.
   */
  inline static int set_max_num_threads(int max_num_threads) {
    int num_threads = 1;
#ifdef HAVE_OPENMP
    if (max_num_threads > 0) {
      omp_set_num_threads(max_num_threads);
    }

#pragma omp parallel
    {
#pragma omp single
      { num_threads = omp_get_num_threads(); }
    }
#endif

#ifdef HAVE_OUTPUT_CYCLES
    // make sure the jobtimes_*.txt files are empty
    for (int i = 0; i < num_threads; ++i) {
      std::stringstream filename;
      filename << "jobtimes_" << i << ".txt";
      std::ofstream file(filename.str(), std::ofstream::trunc);
    }

    // calibrate the cycle count
    Timer timer;
    unsigned int lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    unsigned long first_count = ((unsigned long)hi << 32) | lo;
    timer.start();

    // do some useless work that keeps the cpu busy for a while
    double result = 0.;
    for (unsigned int i = 0; i < 1000000; ++i) {
      double x = 0.01 * i + 0.2;
      result += std::sqrt(x);
    }

    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    unsigned long last_count = ((unsigned long)hi << 32) | lo;
    timer.stop();

    cmac_warning("Time/CPU cycle: %g",
                 timer.value() / (last_count - first_count));
#endif

    return num_threads;
  }
};

#endif // WORKENVIRONMENT_HPP
