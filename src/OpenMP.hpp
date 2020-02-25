/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file OpenMP.hpp
 *
 * @brief Wrapper around OpenMP functionality that also works on systems that
 * don't have OpenMP.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef OPENMP_HPP
#define OPENMP_HPP

#include "Configuration.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>

/**
 * @brief Macro to set the number of threads.
 *
 * @param number_of_threads Number of shared memory parallel threads.
 */
#define set_number_of_threads(number_of_threads)                               \
  omp_set_num_threads(number_of_threads)

/**
 * @brief Get the index of the thread that is executing this piece of code.
 *
 * @return Index of the calling thread.
 */
#define get_thread_index() omp_get_thread_num()
#else

/**
 * @brief Macro to set the number of threads.
 *
 * @param number_of_threads Number of shared memory parallel threads.
 */
#define set_number_of_threads(number_of_threads)

/**
 * @brief Get the index of the thread that is executing this piece of code.
 *
 * @return Index of the calling thread.
 */
#define get_thread_index() 0
#endif

#endif // OPENMP_HPP
