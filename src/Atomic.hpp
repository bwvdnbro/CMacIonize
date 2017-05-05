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
 * @file Atomic.hpp
 *
 * @brief Atomic (lock-free) operations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ATOMIC_HPP
#define ATOMIC_HPP

#include <atomic>

/**
 * @brief Atomic (lock-free) operations.
 */
class Atomic {
public:
  /**
   * @brief Add the second argument to the first argument in a thread safe way.
   *
   * This code was based on the LockFree::add routine in SKIRT.
   *
   * @param a First argument.
   * @param b Second argument.
   */
  template < typename _datatype_ >
  static inline void add(_datatype_ &a, _datatype_ b) {
    std::atomic< _datatype_ > *atom = new (&a) std::atomic< _datatype_ >;
    _datatype_ old = *atom;
    while (!atom->compare_exchange_weak(old, old + b)) {
    }
  }
};

#endif // ATOMIC_HPP
