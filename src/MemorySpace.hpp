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
 * @file MemorySpace.hpp
 *
 * @brief Buffer with pre-allocated PhotonBuffers that can be used.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MEMORYSPACE_HPP
#define MEMORYSPACE_HPP

#include "AtomicValue.hpp"
#include "PhotonBuffer.hpp"
#include "ThreadLock.hpp"
#include "ThreadSafeVector.hpp"

/**
 * @brief Buffer with pre-allocated PhotonBuffers that can be used.
 */
class MemorySpace {
private:
  /*! @brief Memory space itself. */
  ThreadSafeVector< PhotonBuffer > _memory_space;

public:
  /**
   * @brief Constructor.
   *
   * @param size Size of the memory space.
   */
  inline MemorySpace(const size_t size) : _memory_space(size) {}

  /**
   * @brief Access the buffer with the given index.
   *
   * @param index Index of an element in the buffer.
   * @return Corresponding PhotonBuffer.
   */
  inline PhotonBuffer &operator[](const size_t index) {
    // note that the ThreadSafeVector already asserts if this element is in use
    return _memory_space[index];
  }

  /**
   * @brief Get the index of a free buffer in the memory space.
   *
   * @return Index of a free buffer.
   */
  inline size_t get_free_buffer() {
    const size_t index = _memory_space.get_free_element_safe();
    cmac_assert_message(index < _memory_space.max_size(),
                        "No more free elements in memory space!");
    return index;
  }

  /**
   * @brief Free the buffer with the given index.
   *
   * @param index Index of a buffer that was in use.
   */
  inline void free_buffer(const size_t index) {
    _memory_space[index].reset();
    _memory_space.free_element(index);
  }

  /**
   * @brief Copy photons from the given buffer into the buffer with the given
   * index.
   *
   * This method assumes the buffer was successfully locked before it was
   * called. It copies photons until the buffer is full, and starts a new buffer
   * if that happens.
   *
   * @param index Index of the buffer we want to add to.
   * @param buffer Buffer to copy over.
   * @return The index of the last buffer we copied into. If this index is not
   * the same as the original index, the original index should be queued.
   */
  inline size_t add_photons(const size_t index, PhotonBuffer &buffer) {

    PhotonBuffer &buffer_target = _memory_space[index];
    const uint_fast32_t size_in = buffer.size();
    uint_fast32_t counter_in = 0;
    while (buffer_target.size() < PHOTONBUFFER_SIZE && counter_in < size_in) {
      const uint_fast32_t target_index = buffer_target.get_next_free_photon();
      buffer_target[target_index] = buffer[counter_in];
      ++counter_in;
    }
    size_t index_out = index;
    if (buffer_target.size() == PHOTONBUFFER_SIZE) {
      index_out = get_free_buffer();
      // note that the hungry other threads might already be attacking this
      // buffer. We need to make sure they release it without deleting it if
      // it does not (yet) contain any photons.
      PhotonBuffer &buffer_target_new = _memory_space[index_out];
      // copy the old buffer properties
      buffer_target_new.set_subgrid_index(buffer_target.get_subgrid_index());
      buffer_target_new.set_direction(buffer_target.get_direction());
      // maybe assert that this value is indeed 0?
      while (counter_in < size_in) {
        const uint_fast32_t target_index =
            buffer_target_new.get_next_free_photon();
        buffer_target_new[target_index] = buffer[counter_in];
        ++counter_in;
      }
    }
    return index_out;
  }

  /**
   * @brief Get the size of the memory space in memory.
   *
   * @return Size of the memory space in memory (in bytes).
   */
  inline size_t get_memory_size() const {
    return _memory_space.get_memory_size();
  }

#ifdef THREADSAFEVECTOR_STATS
  /**
   * @brief Get the maximum number of elements in the memory space at any given
   * time.
   *
   * @return Maximum number of elements in the memory space.
   */
  inline size_t get_max_number_elements() const {
    return _memory_space.get_max_number_taken();
  }

  /**
   * @brief Reset the counter for the maximum number of elements in the memory
   * space.
   */
  inline void reset_max_number_elements() {
    _memory_space.reset_max_number_taken();
  }
#endif
};

#endif // MEMORYSPACE_HPP
