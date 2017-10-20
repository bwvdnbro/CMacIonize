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
 * @file PhotonBucket.hpp
 *
 * @brief List of PhotonBatch instances for a single small region of the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONBUCKET_HPP
#define PHOTONBUCKET_HPP

#include "Lock.hpp"
#include "PhotonBatch.hpp"

#include <vector>

/**
 * @brief List of PhotonBatch instances for a single small region of the grid.
 */
class PhotonBucket {
private:
  /*! @brief Internal PhotonBatch list. */
  std::vector< PhotonBatch > _batches;

  /*! @brief Maximum size of a single PhotonBatch. */
  uint_least32_t _max_size;

  /*! @brief Lock used to ensure thread safe access to the internal PhotonBatch
   *  list. */
  Lock _lock;

public:
  /**
   * @brief Constructor.
   *
   * @param max_size Maximum size of a single PhotonBatch.
   */
  PhotonBucket(uint_fast32_t max_size) : _max_size(max_size) {}

  /**
   * @brief Add the given Photon to the bucket.
   *
   * This method uses a locking mechanism to guarantee thread safety.
   *
   * @param photon Photon to add.
   */
  void add_photon(Photon photon) {
    // lock the bucket: everything below can only be done by one thread at a
    // time
    _lock.lock();

    if (_batches.size() == 0) {
      // create the first batch
      _batches.push_back(PhotonBatch(_max_size));
    }
    bool success = _batches.back().add_photon(photon);
    if (!success) {
      // create a new batch
      _batches.push_back(PhotonBatch(_max_size));
      success = _batches.back().add_photon(photon);
      if (!success) {
        cmac_error("Something went wrong while adding a photon to a bucket!");
      }
    }

    // ready changing the batches: unlock the bucket
    _lock.unlock();
  }

  /**
   * @brief Get a PhotonBatch to process.
   *
   * @return Pointer to a PhotonBatch that can be processed. Memory management
   * for this pointer should be done by the calling process.
   */
  PhotonBatch *get_batch() {
    PhotonBatch *first_batch = nullptr;

    // lock the bucket: we change the internal array
    _lock.lock();

    if (_batches.size() > 0) {
      first_batch = new PhotonBatch(_batches[0]);
      _batches.erase(_batches.begin());
    }

    // ready changing the internal array
    _lock.unlock();

    return first_batch;
  }
};

#endif // PHOTONBUCKET_HPP
