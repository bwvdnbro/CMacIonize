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
   * @param max_num_threads Maximum number of threads to use.
   */
  inline static void set_max_num_threads(int max_num_threads) {
#ifdef HAVE_OPENMP
    omp_set_num_threads(max_num_threads);
#endif

#ifdef HAVE_OUTPUT_CYCLES
    // make sure the jobtimes_*.txt files are empty
    for (int i = 0; i < max_num_threads; ++i) {
      std::stringstream filename;
      filename << "jobtimes_" << i << ".txt";
      std::ofstream file(filename.str(), std::ofstream::trunc);
    }
#endif
  }
};

#endif // WORKENVIRONMENT_HPP
