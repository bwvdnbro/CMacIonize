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
 * @file Lock.hpp
 *
 * @brief Lock used to ensure secure access to variables in a shared memory
 * parallel context.
 *
 * This class and the WorkDistributor class are the only classes that should
 * explicitly contain OpenMP statements.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LOCK_HPP
#define LOCK_HPP

#include "Configuration.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/**
 * @brief Lock used to ensure secure access to variables in a shared memory
 * parallel context.
 *
 * If OpenMP is not available, this class does nothing at all.
 */
class Lock {
private:
#ifdef HAVE_OPENMP
  /*! @brief Underlying OpenMP lock. */
  omp_lock_t _lock;
#endif

public:
  /**
   * @brief Constructor.
   */
  inline Lock() {
#ifdef HAVE_OPENMP
    omp_init_lock(&_lock);
#endif
  }

  inline ~Lock() {
#ifdef HAVE_OPENMP
    omp_destroy_lock(&_lock);
#endif
  }

  /**
   * @brief Obtain exclusive access of the Lock.
   *
   * This routine tries to get exclusive access of the lock. If the lock is
   * already in use by another thread, the calling thread will hang until the
   * lock is released.
   */
  inline void lock() {
#ifdef HAVE_OPENMP
    omp_set_lock(&_lock);
#endif
  }

  /**
   * @brief Release exclusive access of the Lock.
   *
   * This routine will release the lock. If another thread has been waiting for
   * exclusive access to the lock, that thread will now be able to get access.
   */
  inline void unlock() {
#ifdef HAVE_OPENMP
    omp_unset_lock(&_lock);
#endif
  }
};

#endif // LOCK_HPP
