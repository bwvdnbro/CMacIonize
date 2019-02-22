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
 * @file ThreadSafeVector.hpp
 *
 * @brief Thread safe fixed size vector.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef THREADSAFEVECTOR_HPP
#define THREADSAFEVECTOR_HPP

#include "AtomicValue.hpp"
#include "Error.hpp"

#include <string>

/*! @brief Total size of the variables whose size is known at compile time. */
#define THREADSAFEVECTOR_FIXED_SIZE sizeof(ThreadSafeVector< _datatype_ >)

/*! @brief Size per element of the variables whose size is unknown at compile
 *  time. */
#define THREADSAFEVECTOR_ELEMENT_SIZE                                          \
  (sizeof(_datatype_) + sizeof(AtomicValue< bool >))

/**
 * @brief Thread safe fixed size vector.
 */
template < typename _datatype_ > class ThreadSafeVector {
private:
  /*! @brief Current active index in the vector. */
  AtomicValue< size_t > _current_index;

  /*! @brief Size of the vector. */
  const size_t _size;

  /*! @brief Number of elements that have been taken. */
  AtomicValue< size_t > _number_taken;

  /*! @brief Vector itself. */
  _datatype_ *_vector;

  /*! @brief Atomic locks for the vector elements. */
  AtomicValue< bool > *_locks;

#ifdef THREADSAFEVECTOR_STATS
  /*! @brief Maximum number of elements taken simultaneously at any given point
   * in time. */
  AtomicValue< size_t > _max_number_taken;
#endif

  /*! @brief Label appended to error messages. */
  const std::string _label;

public:
  /**
   * @brief Constructor.
   *
   * @param size Size of the vector.
   * @param label Label used in error messages.
   */
  inline ThreadSafeVector(const size_t size, const std::string label = "")
      : _current_index(0), _size(size), _number_taken(0), _label(label) {
    _vector = new _datatype_[size];
    // Atomic values are automatically initialized to 0 (= false) by the
    // constructor
    _locks = new AtomicValue< bool >[size];
  }

  /**
   * @brief Destructor.
   */
  inline ~ThreadSafeVector() {
    delete[] _vector;
    delete[] _locks;
  }

  /**
   * @brief Clear the contents of the vector.
   *
   * This method is not meant to be thread safe.
   */
  inline void clear() {
    for (size_t i = 0; i < _size; ++i) {
      _locks[i].unlock();
    }
    _number_taken.set(0);
    _current_index.set(0);

    // clear the elements
    delete[] _vector;
    _vector = new _datatype_[_size];

#ifdef THREADSAFEVECTOR_STATS
    _max_number_taken.set(0);
#endif
  }

  /**
   * @brief Get a continuous block of elements with the given size at the start
   * of the vector.
   *
   * This method is not meant to be thread safe.
   *
   * @param size Size of the block.
   */
  inline void get_free_elements(const size_t size) {

    cmac_assert_message(size < _size,
                        "Not enough elements available (%zu < %zu)! (%s)", size,
                        _size, _label.c_str());

    for (size_t i = 0; i < size; ++i) {
      _locks[i].lock();
    }
    _current_index.set(size);
    _number_taken.set(size);

#ifdef THREADSAFEVECTOR_STATS
    _max_number_taken.max(_number_taken.value());
#endif
  }

  /**
   * @brief Clear the contents of the vector from the given offset index
   * onwards.
   *
   * This method is not meant to be thread safe. We assume all values before the
   * given offset are in use.
   *
   * @param offset First index that should be cleared.
   */
  inline void clear_after(const size_t offset) {
    for (size_t i = offset; i < _size; ++i) {
      _locks[i].unlock();
      _vector[i] = _datatype_();
    }
    _number_taken.set(offset);
    _current_index.set(offset);

#ifdef THREADSAFEVECTOR_STATS
    _max_number_taken.set(offset);
#endif
  }

  /**
   * @brief Access the element with the given index.
   *
   * @param index Index of an element.
   * @return Read/write reference to the element with that index.
   */
  inline _datatype_ &operator[](const size_t index) {
    cmac_assert_message(index < _size,
                        "Element out of range (index: %zu, size: %zu)! (%s)",
                        index, _size, _label.c_str());
    cmac_assert_message(_locks[index].value(),
                        "Element not in use (index: %zu)! (%s)", index,
                        _label.c_str());
    return _vector[index];
  }

  /**
   * @brief Read-only access to the element with the given index.
   *
   * @param index Index of an element.
   * @return Read-only reference to the element with that index.
   */
  inline const _datatype_ &operator[](const size_t index) const {
    cmac_assert_message(index < _size,
                        "Element out of range (index: %zu, size: %zu)! (%s)",
                        index, _size, _label.c_str());
    cmac_assert_message(_locks[index].value(),
                        "Element not in use (index: %zu)! (%s)", index,
                        _label.c_str());
    return _vector[index];
  }

  /**
   * @brief Get the index of a free element in the vector.
   *
   * This element will be locked and needs to be freed later by calling
   * free_element().
   *
   * @return Index of a free element.
   */
  inline size_t get_free_element() {
    cmac_assert_message(_number_taken.value() < _size,
                        "No more free elements in vector (%zu < %zu)! (%s)",
                        _number_taken.value(), _size, _label.c_str());
    size_t index = _current_index.post_increment() % _size;
    while (!_locks[index].lock()) {
      index = _current_index.post_increment() % _size;
    }
#ifdef THREADSAFEVECTOR_STATS
    const size_t number_taken = _number_taken.pre_increment();
    _max_number_taken.max(number_taken);
#else
    _number_taken.pre_increment();
#endif
    return index;
  }

  /**
   * @brief Get the index of a free element in the vector, if available.
   *
   * This element will be locked and needs to be freed later by calling
   * free_element(). If no free element is available, this function returns
   * max_size().
   *
   * @return Index of a free element.
   */
  inline size_t get_free_element_safe() {
    if (_number_taken.value() < _size) {
      size_t index = _current_index.post_increment() % _size;
      while (!_locks[index].lock()) {
        index = _current_index.post_increment() % _size;
      }
#ifdef THREADSAFEVECTOR_STATS
      const size_t number_taken = _number_taken.pre_increment();
      _max_number_taken.max(number_taken);
#else
      _number_taken.pre_increment();
#endif
      return index;
    } else {
      return _size;
    }
  }

  /**
   * @brief Free the element with the given index.
   *
   * The element can be overwritten after this method has been called.
   *
   * @param index Index of an element that was in use.
   */
  inline void free_element(const size_t index) {
    cmac_assert_message(_locks[index].value(),
                        "Element not in use (index: %zu)! (%s)", index,
                        _label.c_str());
    _locks[index].unlock();
    _number_taken.pre_decrement();
  }

  /**
   * @brief Get the size of the vector.
   *
   * Only makes sense if the vector is continuous, i.e. non of the elements has
   * ever been removed.
   *
   * @return Number of used elements in the vector.
   */
  inline size_t size() const {
    cmac_assert_message(_number_taken.value() == _current_index.value(),
                        "Non continuous vector (%zu =/= %zu)! (%s)",
                        _number_taken.value(), _current_index.value(),
                        _label.c_str());
    return _number_taken.value();
  }

  /**
   * @brief Maximum number of elements in the vector.
   *
   * @return Maximum number of elements that can be stored in the vector.
   */
  inline size_t max_size() const { return _size; }

  /**
   * @brief Get the size in memory of the vector.
   *
   * @return Size in memory of the vector (in bytes).
   */
  inline size_t get_memory_size() const {
    return THREADSAFEVECTOR_FIXED_SIZE + _size * THREADSAFEVECTOR_ELEMENT_SIZE;
  }

/**
 * @brief Get the maximum number of elements that was taken in this vector at
 * any given time.
 *
 * @return Maximum number of elements that was taken.
 */
#ifdef THREADSAFEVECTOR_STATS
  inline size_t get_max_number_taken() const {
    return _max_number_taken.value();
  }
#endif

/**
 * @brief Reset the counter for the maximum number of elements that was taken.
 */
#ifdef THREADSAFEVECTOR_STATS
  inline void reset_max_number_taken() { _max_number_taken.set(0); }
#endif
};

#endif // THREADSAFEVECTOR_HPP
