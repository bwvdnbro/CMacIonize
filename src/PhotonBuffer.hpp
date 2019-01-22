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
 * @file PhotonBuffer.hpp
 *
 * @brief Photon buffer.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#ifndef PHOTONBUFFER_HPP
#define PHOTONBUFFER_HPP

// project includes
#include "AtomicValue.hpp"
#include "Configuration.hpp"
#include "Error.hpp"
#include "PhotonPacket.hpp"
#include "ThreadLock.hpp"

// standard library includes
#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*! @brief Number of photons that can be stored in a single buffer. */
#define PHOTONBUFFER_SIZE 200u

/*! @brief Size of the MPI buffer necessary to store a PhotonBuffer. */
#define PHOTONBUFFER_MPI_SIZE                                                  \
  (sizeof(size_t) + sizeof(uint_least32_t) + sizeof(int_least32_t) +           \
   PHOTONBUFFER_SIZE * PHOTON_MPI_SIZE)

/**
 * @brief Photon buffer.
 */
class PhotonBuffer {
private:
  /*! @brief Subgrid with which this buffer is associated. */
  size_t _subgrid_index;

  /*! @brief TravelDirection of the photons in the buffer. */
  int_least32_t _direction;

  /*! @brief Number of photons in the buffer. */
  uint_least32_t _actual_size;

  /*! @brief Actual photon buffer. */
  PhotonPacket _photons[PHOTONBUFFER_SIZE];

public:
  /**
   * @brief Empty constructor.
   */
  PhotonBuffer() : _actual_size(0) {}

#ifdef HAVE_MPI
  /**
   * @brief Store the contents of the PhotonBuffer in the given MPI
   * communication buffer.
   *
   * @param buffer Buffer to use (should be preallocated and have at least size
   * PHOTONBUFFER_MPI_SIZE).
   */
  inline void pack(char buffer[PHOTONBUFFER_MPI_SIZE]) {
    int buffer_position = 0;
    MPI_Pack(&_subgrid_index, 1, MPI_UNSIGNED, buffer, PHOTONBUFFER_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_direction, 1, MPI_INT, buffer, PHOTONBUFFER_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_actual_size, 1, MPI_UNSIGNED, buffer, PHOTONBUFFER_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    for (uint_fast32_t i = 0; i < PHOTONBUFFER_SIZE; ++i) {
      _photons[i].pack(&buffer[buffer_position]);
      buffer_position += PHOTON_MPI_SIZE;
    }
  }

  /**
   * @brief Copy the contents of the given MPI communication buffer into this
   * PhotonBuffer.
   *
   * @param buffer MPI communication buffer (should have at least size
   * PHOTONBUFFER_MPI_SIZE).
   */
  inline void unpack(char buffer[PHOTONBUFFER_MPI_SIZE]) {
    int buffer_position = 0;
    MPI_Unpack(buffer, PHOTONBUFFER_MPI_SIZE, &buffer_position, &_subgrid_index,
               1, MPI_UNSIGNED, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTONBUFFER_MPI_SIZE, &buffer_position, &_direction, 1,
               MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTONBUFFER_MPI_SIZE, &buffer_position, &_actual_size,
               1, MPI_UNSIGNED, MPI_COMM_WORLD);
    for (uint_fast32_t i = 0; i < PHOTONBUFFER_SIZE; ++i) {
      _photons[i].unpack(&buffer[buffer_position]);
      buffer_position += PHOTON_MPI_SIZE;
    }
  }
#endif

  /**
   * @brief Check that the given PhotonBuffer is equal to this one.
   *
   * @param other Other PhotonBuffer.
   */
  inline void check_equal(const PhotonBuffer &other) {
    cmac_assert_message(_subgrid_index == other._subgrid_index,
                        "Subgrid indices do not match!");
    cmac_assert_message(_direction == other._direction,
                        "Directions do not match!");
    cmac_assert_message(_actual_size == other._actual_size,
                        "Sizes do not match!");
    for (uint_fast32_t i = 0; i < PHOTONBUFFER_SIZE; ++i) {
      _photons[i].check_equal(other._photons[i]);
    }
  }

  /**
   * @brief Get the size of the buffer.
   *
   * @return Number of photons in the buffer.
   */
  inline uint_fast32_t size() const { return _actual_size; }

  /**
   * @brief Reset the buffer by setting its size to 0.
   */
  inline void reset() { _actual_size = 0; }

  /**
   * @brief Get read-only access to the Photon with the given index.
   *
   * @param index Index of a photon in the buffer.
   * @return Read-only reference to the Photon that corresponds to that index.
   */
  inline const PhotonPacket &operator[](const uint_fast32_t index) const {
    return _photons[index];
  }

  /**
   * @brief Get read/write access to the Photon with the given index.
   *
   * @param index Index of a photon in the buffer.
   * @return Reference to the Photon that corresponds to that index.
   */
  inline PhotonPacket &operator[](const uint_fast32_t index) {
    return _photons[index];
  }

  /**
   * @brief Get the index of the next available Photon in the buffer.
   *
   * @return Index of the next available Photon in the buffer.
   */
  inline uint_fast32_t get_next_free_photon() {
    cmac_assert_message(_actual_size < PHOTONBUFFER_SIZE,
                        "No more free elements in PhotonBuffer!");
    return _actual_size++;
  }

  /**
   * @brief Grow the buffer to the given size.
   *
   * New photons are not initialized by this method.
   *
   * @param size New buffer size.
   */
  inline void grow(const uint_fast32_t size) { _actual_size = size; }

  /**
   * @brief Set the subgrid index for this buffer.
   *
   * @param subgrid_index Subgrid index for this buffer.
   */
  inline void set_subgrid_index(const size_t subgrid_index) {
    _subgrid_index = subgrid_index;
  }

  /**
   * @brief Get the subgrid index for this buffer.
   *
   * @return Subgrid index for this buffer.
   */
  inline size_t get_subgrid_index() const { return _subgrid_index; }

  /**
   * @brief Set the direction for this buffer.
   *
   * @param direction Direction for this buffer.
   */
  inline void set_direction(const int_least32_t direction) {
    _direction = direction;
  }

  /**
   * @brief Get the direction for this buffer.
   *
   * @return Direction for this buffer.
   */
  inline int_least32_t get_direction() const { return _direction; }
};

#endif // PHOTONBUFFER_HPP
