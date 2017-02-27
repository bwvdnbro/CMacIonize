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
 * @file PhotonBatch.hpp
 *
 * @brief Batch of Photons that is propagated by a single thread on a single
 * process, on a small, local part of the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONBATCH_HPP
#define PHOTONBATCH_HPP

#include "Photon.hpp"

#include <vector>

/**
 * @brief Batch of Photons that is propagated by a single thread on a single
 * process, on a small, local part of the DensityGrid.
 */
class PhotonBatch {
private:
  /*! @brief List of Photons contained in this batch. */
  std::vector< Photon > _photons;

  /*! @brief Indices of the next sub regions in which the photons need to be
   *  propagated. */
  std::vector< int > _indices;

  /*! @brief Maximum size of the internal arrays. */
  unsigned int _size;

public:
  /**
   * @brief Constructor.
   *
   * @param size Maximum size of the internal arrays.
   */
  PhotonBatch(unsigned int size) : _size(size) {
    _photons.reserve(size);
    _indices.reserve(size);
  }

  /**
   * @brief Add the given Photon to the internal list.
   *
   * @param photon Photon to add.
   * @return True if the Photon was successfully added.
   */
  bool add_photon(Photon photon) {
    if (_photons.size() < _size) {
      _photons.push_back(photon);
      _indices.push_back(-1);
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief Iterator to loop over the Photons contained in the PhotonBatch.
   */
  class iterator {
  private:
    /*! @brief Reference to the PhotonBatch over which we loop. */
    PhotonBatch &_batch;

    /*! @brief Index of the element the iterator is currently pointing to. */
    unsigned int _index;

  public:
    /// Basic iterator functionality

    /**
     * @brief Constructor.
     *
     * @param batch Reference to the PhotonBatch over which we loop.
     * @param index Index of the element the iterator should point to.
     */
    iterator(PhotonBatch &batch, unsigned int index)
        : _batch(batch), _index(index) {}

    /**
     * @brief Increment operator.
     *
     * @return Reference to the incremented iterator.
     */
    iterator &operator++() {
      ++_index;
      return *this;
    }

    /**
     * @brief Comparison operator.
     *
     * @param it iterator to compare with.
     * @return True if both iterators point to the same element of the same
     * PhotonBatch.
     */
    bool operator==(iterator it) const {
      return (&_batch == &it._batch && _index == it._index);
    }

    /**
     * @brief Comparison operator.
     *
     * @param it iterator to compare with.
     * @return True if both iterators point to another element, or if they point
     * to another PhotonBatch.
     */
    bool operator!=(iterator it) const { return !(*this == it); }

    /// Extra iterator methods

    /**
     * @brief Get the Photon the iterator is currently pointing to.
     *
     * @return Photon iterator is currently pointing to.
     */
    Photon get_photon() { return _batch._photons[_index]; }

    /**
     * @brief Set the new index for the Photon the iterator is currently
     * pointing to.
     *
     * The new index is either the index of the next region of the grid through
     * which the Photon should be propagated, or -1 if the Photon leaves the
     * system.
     *
     * @param index New index for the Photon the iterator is pointing to.
     */
    void set_new_index(int index) { _batch._indices[_index] = index; }

    /**
     * @brief Get the new index for the Photon the iterator is currently
     * pointing to.
     *
     * The new index is either the index of the next region of the grid through
     * which the Photon should be propagated, or -1 if the Photon leaves the
     * system.
     *
     * @return New index for the Photon the iterator is pointing to.
     */
    int get_new_index() const { return _batch._indices[_index]; }
  };

  /**
   * @brief Get an iterator to the first Photon in the batch.
   *
   * @return iterator to the first Photon in the batch.
   */
  iterator begin() { return iterator(*this, 0); }

  /**
   * @brief Get an iterator to the beyond last Photon in the batch.
   *
   * @return iterator to the beyond last Photon in the batch.
   */
  iterator end() { return iterator(*this, _photons.size()); }
};

#endif // PHOTONBATCH_HPP
