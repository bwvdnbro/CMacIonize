/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AtomicValue.hpp
 *
 * @brief Value that can only be changed by thread-safe atomic operations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ATOMICVALUE_HPP
#define ATOMICVALUE_HPP

/*! @brief Activate this to use standard C++11 atomic operations. */
#define CPP_ATOMIC

/*! @brief Activate this to use GNU GCC specific atomic directives. */
//#define GCC_ATOMIC

#if defined(CPP_ATOMIC)
#include <atomic>
#elif defined(GCC_ATOMIC)
// nothing to include
#else
#error "No atomic implementation chosen!"
#endif

#include <algorithm>

/**
 * @brief Value that can only be changed by thread-safe atomic operations.
 */
template < typename _type_ > class AtomicValue {
private:
/*! @brief Underlying value. */
#if defined(CPP_ATOMIC)
  std::atomic< _type_ > _value;
#elif defined(GCC_ATOMIC)
  volatile _type_ _value;
#endif

public:
  /**
   * @brief Empty constructor.
   */
  inline AtomicValue() : _value(0) {}

  /**
   * @brief Constructor.
   *
   * @param value Value.
   */
  inline AtomicValue(_type_ value) : _value(value) {}

  /**
   * @brief Read the value of the atomic variable.
   *
   * @return Current value of the variable.
   */
  inline const _type_ value() const {
#if defined(CPP_ATOMIC)
    return _value.load();
#elif defined(GCC_ATOMIC)
    return _value;
#endif
  }

  /**
   * @brief Set the value of the atomic variable.
   *
   * @param value New value for the variable.
   */
  inline void set(const _type_ value) {
#if defined(CPP_ATOMIC)
    _value.store(value);
#elif defined(GCC_ATOMIC)
    _value = value;
#endif
  }

  /**
   * @brief Lock the value atomically, making sure only one thread is allowed to
   * set it.
   *
   * @return True if the flag was set, False if it could not be set, meaning the
   * flag has already been set by another thread.
   */
  inline bool lock() {
#if defined(CPP_ATOMIC)
    _type_ expected(false);
    return _value.compare_exchange_strong(expected, true);
#elif defined(GCC_ATOMIC)
    return __sync_bool_compare_and_swap(&_value, false, true);
#endif
  }

  /**
   * @brief Unlock the value atomically, making sure only one thread is allowed
   * to set it.
   */
  inline void unlock() {
#if defined(CPP_ATOMIC)
    _type_ expected(true);
    _value.compare_exchange_strong(expected, false);
#elif defined(GCC_ATOMIC)
    __sync_bool_compare_and_swap(&_value, true, false);
#endif
  }

  /**
   * @brief Atomically increment the variable by 1 and return the original
   * value.
   *
   * @return Original value of the variable.
   */
  inline _type_ post_increment() {
#if defined(CPP_ATOMIC)
    return _value++;
#elif defined(GCC_ATOMIC)
    return __sync_fetch_and_add(&_value, 1);
#endif
  }

  /**
   * @brief Atomically increment the variable by 1 and return the new value.
   *
   * @return New value of the variable.
   */
  inline _type_ pre_increment() {
#if defined(CPP_ATOMIC)
    return ++_value;
#elif defined(GCC_ATOMIC)
    return __sync_add_and_fetch(&_value, 1);
#endif
  }

  /**
   * @brief Atomically decrement the variable by 1 and return the new value.
   *
   * @return New value of the variable.
   */
  inline _type_ pre_decrement() {
#if defined(CPP_ATOMIC)
    return --_value;
#elif defined(GCC_ATOMIC)
    return __sync_sub_and_fetch(&_value, 1);
#endif
  }

  /**
   * @brief Atomically add the given value to the variable and return the old
   * value of the variable.
   *
   * @param increment Value to add to the variable.
   * @return Original value of the variable.
   */
  inline _type_ post_add(const _type_ increment) {
#if defined(CPP_ATOMIC)
    // we cannot use (_value += increment), as the Intel compiler seems to
    // replace this (incorrectly) with an add_fetch rather than a fetch_add...
    return _value.fetch_add(increment);
#elif defined(GCC_ATOMIC)
    return __sync_fetch_and_add(&_value, increment);
#endif
  }

  /**
   * @brief Atomically add the given value to the variable and return the new
   * value of the variable.
   *
   * @param increment Value to add to the variable.
   * @return New value of the variable.
   */
  inline _type_ pre_add(const _type_ increment) {
#if defined(CPP_ATOMIC)
    // std::atomic does not have direct support for add_fetch, so we have to
    // emulate it
    return _value.fetch_add(increment) + increment;
#elif defined(GCC_ATOMIC)
    return __sync_add_and_fetch(&_value, increment);
#endif
  }

  /**
   * @brief Atomically subtract the given value from the variable and return the
   * new value of the variable.
   *
   * @param decrement Value to subtract from the variable.
   * @return New value of the variable.
   */
  inline _type_ pre_subtract(const _type_ decrement) {
#if defined(CPP_ATOMIC)
    // std::atomic does not have direct support for add_fetch, so we have to
    // emulate it
    return _value.fetch_sub(decrement) - decrement;
#elif defined(GCC_ATOMIC)
    return __sync_sub_and_fetch(&_value, decrement);
#endif
  }

  /**
   * @brief Atomically update the variable with the maximum of its current value
   * and the given value.
   *
   * @param value Potential new maximum for the variable.
   */
  inline void max(const _type_ value) {
#if defined(CPP_ATOMIC)
    _type_ old_value = _value.load();
    _type_ new_value = std::max(value, old_value);
    while (!_value.compare_exchange_strong(old_value, new_value)) {
      old_value = _value.load();
      new_value = std::max(value, old_value);
    }
#elif defined(GCC_ATOMIC)
    _type_ old_value = _value;
    _type_ new_value = std::max(value, old_value);
    while (!__sync_bool_compare_and_swap(&_value, old_value, new_value)) {
      old_value = _value;
      new_value = std::max(value, old_value);
    }
#endif
  }
};

#endif // ATOMICVALUE_HPP
