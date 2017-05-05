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
 * @file PhotonPool.hpp
 *
 * @brief System wide pool of photons across the entire grid. Contains separate
 * PhotonBucket instances for each sub region, which in turn contain small
 * PhotonBatch instances for each Worker.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONPOOL_HPP
#define PHOTONPOOL_HPP

#include "PhotonBucket.hpp"

#include <vector>

/**
 * @brief System wide pool of photons across the entire grid. Contains separate
 * PhotonBucket instances for each sub region, which in turn contain small
 * PhotonBatch instances for each Worker.
 */
class PhotonPool {
private:
  /*! @brief PhotonBuckets for each sub region of the grid. */
  std::vector< PhotonBucket > _buckets;

public:
  /**
   * @brief Constructor.
   *
   * @param num_bucket Number of buckets in the pool.
   * @param max_size Maximum size of a single PhotonBatch in a bucket.
   */
  PhotonPool(unsigned int num_bucket, unsigned int max_size)
      : _buckets(num_bucket, max_size) {}

  /**
   * @brief Add a photon to the bucket with the given index.
   *
   * @param bucket_index Index of a PhotonBucket.
   * @param photon Photon to add.
   */
  void add_photon(unsigned int bucket_index, Photon photon) {
    _buckets[bucket_index].add_photon(photon);
  }

  /**
   * @brief Get a PhotonBatch to process.
   *
   * @param bucket_index Index of a PhotonBucket.
   * @return Pointer to a PhotonBatch that can be processed. Memory management
   * for this pointer should be done by the calling process.
   */
  PhotonBatch *get_batch(unsigned int bucket_index) {
    return _buckets[bucket_index].get_batch();
  }
};

#endif // PHOTONPOOL_HPP
