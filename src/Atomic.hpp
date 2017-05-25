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

#include "Configuration.hpp"

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

#ifndef HAVE_ATOMIC
/**
 * @brief Add the second value to the first value, assuming the first value is
 * a binary cast of a double precision floating point value.
 *
 * @param a Binary cast of double precision floating point value.
 * @param b Double precision floating point value to add.
 * @return Double precision floating point sum of a and b, binary cast back to
 * an unsigned long integer type.
 */
static unsigned long add_floating_point(unsigned long a, double b) {
  union {
    unsigned long l;
    double d;
  } vala, valb;
  vala.l = a;
  valb.d = vala.d + b;
  return valb.l;
}

/**
 * @brief Atomic::add specialization for double precision floating point values
 * for systems that do not support atomic floating point operations.
 *
 * @param a Reference to the double precision floating point value that should
 * be atomically incremented.
 * @param b Double precision floating point value that should be added to a.
 */
template <> inline void Atomic::add< double >(double &a, double b) {
  std::atomic< unsigned long > *atom = new (&a) std::atomic< unsigned long >;
  unsigned long old = *atom;
  while (!atom->compare_exchange_weak(old, add_floating_point(old, b))) {
  }
}
#endif

#endif // ATOMIC_HPP
