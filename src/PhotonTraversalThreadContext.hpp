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
 * @file PhotonTraversalThreadContext.hpp
 *
 * @brief Thread context responsible for storing thread-local data structures
 * used by the photon propagation task.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PHOTONTRAVERSALTHREADCONTEXT_HPP
#define PHOTONTRAVERSALTHREADCONTEXT_HPP

#include "DensitySubGrid.hpp"
#include "PhotonBuffer.hpp"
#include "TravelDirections.hpp"

/**
 * @brief Thread context responsible for storing thread-local data structures
 * used by the photon propagation task.
 */
class PhotonTraversalThreadContext {
private:
  /*! @brief Buffers used to store outgoing photon packets during grid
   *  traversal. */
  PhotonBuffer _local_buffers[TRAVELDIRECTION_NUMBER];

  /*! @brief Flags used to distinguish active and inactive directions. */
  bool _local_buffer_flags[TRAVELDIRECTION_NUMBER];

public:
  /**
   * @brief Constructor.
   */
  inline PhotonTraversalThreadContext() {
    for (int_fast8_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      _local_buffers[i].set_direction(
          TravelDirections::output_to_input_direction(i));
      _local_buffers[i].reset();
      _local_buffer_flags[i] = true;
    }
  }

  /**
   * @brief Get a pointer to the local buffer array.
   *
   * @return Pointer to the local buffer array.
   */
  inline PhotonBuffer *get_local_buffers() { return _local_buffers; }

  /**
   * @brief Get a pointer to the local buffer flag array.
   *
   * @return Pointer to the local buffer flag array.
   */
  inline bool *get_local_buffer_flags() { return _local_buffer_flags; }

  /**
   * @brief Initialize the context for traversal of the given subgrid.
   *
   * This function disables all outgoing buffers that correspond to directions
   * that are not inside the simulation box. If reemission is disabled, the
   * local outgoing buffer that stores absorbed photon packets is also disabled.
   *
   * This function also ensures that all active outgoing buffers are correctly
   * initialised.
   *
   * @param grid Subgrid to traverse.
   * @param do_reemission Whether or not diffuse reemission is active.
   */
  inline void initialize(const DensitySubGrid &grid, const bool do_reemission) {

    for (int_fast8_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      const uint_fast32_t ngb = grid.get_neighbour(i);
      if (ngb != NEIGHBOUR_OUTSIDE) {
        _local_buffer_flags[i] = true;
        _local_buffers[i].reset();
      } else {
        _local_buffer_flags[i] = false;
      }
    }
    // if reemission is disabled, disable output to the internal buffer
    if (!do_reemission) {
      _local_buffer_flags[TRAVELDIRECTION_INSIDE] = false;
    }
  }

  /**
   * @brief Store the given photon packet that leaves the subgrid in the given
   * direction.
   *
   * If the buffer for the given direction does not exist, the photon is not
   * stored.
   *
   * @param output_direction TravelDirection of the outgoing photon packet.
   * @param photon Photon packet.
   * @return True if the photon was stored, false otherwise.
   */
  inline bool store_photon(const int_fast32_t output_direction,
                           const PhotonPacket &photon) {

    if (_local_buffer_flags[output_direction]) {
      // get the correct output buffer
      PhotonBuffer &output_buffer = _local_buffers[output_direction];

      // add the photon
      const uint_fast32_t index = output_buffer.get_next_free_photon();
      output_buffer[index] = photon;
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief Check if the buffer for the given outgoing direction contains any
   * photon packets.
   *
   * @param output_direction TravelDirection for outgoing photons.
   * @return True if the corresponding buffer contains photons.
   */
  inline bool has_outgoing_photons(const int_fast32_t output_direction) const {
    return _local_buffer_flags[output_direction] &&
           _local_buffers[output_direction].size() > 0;
  }

  /**
   * @brief Get the buffer corresponding to the given outgoing direction.
   *
   * @param output_direction TravelDirection for outgoing photons.
   * @return Buffer for that direction.
   */
  inline const PhotonBuffer &
  get_outgoing_buffer(const int_fast32_t output_direction) const {
    return _local_buffers[output_direction];
  }
};

#endif // PHOTONTRAVERSALTHREADCONTEXT_HPP
