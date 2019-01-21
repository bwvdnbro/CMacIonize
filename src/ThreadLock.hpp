/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017,2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ThreadLock.hpp
 *
 * @brief Implementation independent lock.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef THREADLOCK_HPP
#define THREADLOCK_HPP

/*! @brief Use OpenMP locks. */
//#define LOCK_OPENMP

/*! @brief Use Atomic locks. */
#define LOCK_ATOMIC

/*! @brief Use PThread spinlocks. */
//#define LOCK_PTHREAD_SPIN

#if defined(LOCK_OPENMP)
#include <omp.h>
#elif defined(LOCK_ATOMIC)
#include "AtomicValue.hpp"
#elif defined(LOCK_PTHREAD_SPIN)
#include <pthread.h>
#else
#error "No locking implementation chosen!"
#endif

/**
 * @brief Implementation independent lock.
 */
class ThreadLock {
private:
/*! @brief Lock variable. */
#if defined(LOCK_OPENMP)
  omp_lock_t _lock;
#elif defined(LOCK_ATOMIC)
  AtomicValue< bool > _lock;
#elif defined(LOCK_PTHREAD_SPIN)
  pthread_spinlock_t _lock;
#endif

public:
  /**
   * @brief Constructor.
   */
  inline ThreadLock()
#if defined(LOCK_ATOMIC)
      : _lock(false)
#endif
  {
#if defined(LOCK_OPENMP)
    omp_init_lock(&_lock);
#elif defined(LOCK_PTHREAD_SPIN)
    pthread_spin_init(&_lock, PTHREAD_PROCESS_PRIVATE);
#endif
  }

  /**
   * @brief Destructor.
   */
  inline ~ThreadLock() {
#if defined(LOCK_OPENMP)
    omp_destroy_lock(&_lock);
#elif defined(LOCK_PTHREAD_SPIN)
    pthread_spin_destroy(&_lock);
#endif
  }

  /**
   * @brief Try to lock the lock.
   *
   * @return True if the lock succeeded.
   */
  inline bool try_lock() {
#if defined(LOCK_OPENMP)
    return omp_test_lock(&_lock);
#elif defined(LOCK_ATOMIC)
    return _lock.lock();
#elif defined(LOCK_PTHREAD_SPIN)
    return pthread_spin_trylock(&_lock);
#endif
  }

  /**
   * @brief Lock the lock.
   *
   * Hangs until the lock succeeds.
   */
  inline void lock() {
#if defined(LOCK_OPENMP)
    omp_set_lock(&_lock);
#elif defined(LOCK_ATOMIC)
    while (!_lock.lock()) {
    }
#elif defined(LOCK_PTHREAD_SPIN)
    pthread_spin_lock(&_lock);
#endif
  }

  /**
   * @brief Unlock the lock.
   */
  inline void unlock() {
#if defined(LOCK_OPENMP)
    omp_unset_lock(&_lock);
#elif defined(LOCK_ATOMIC)
    _lock.unlock();
#elif defined(LOCK_PTHREAD_SPIN)
    pthread_spin_unlock(&_lock);
#endif
  }
};

#endif // THREADLOCK_HPP
