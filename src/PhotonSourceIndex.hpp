/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file PhotonSourceIndex.hpp
 *
 * @brief Index pointing to a source in the PhotonSource.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCEINDEX_HPP
#define PHOTONSOURCEINDEX_HPP

/**
 * @brief Index pointing to a source in the PhotonSource.
 */
class PhotonSourceIndex {
private:
  /*! @brief Number of photons emitted by the currently active discrete source.
   */
  unsigned int _discrete_active_number_of_photons;

  /*! @brief Currently emitted discrete photon index. */
  unsigned int _discrete_active_photon_index;

  /*! @brief Currently active discrete photon source index. */
  unsigned int _discrete_active_source_index;

  /*! @brief Number of photons already emitted by the continuous source. */
  unsigned int _continuous_active_number_of_photons;

public:
  /**
   * @brief Constructor.
   *
   * @param discrete_active_number_of_photons Number of photons emitted by the
   * currently active discrete source.
   * @param discrete_active_photon_index Number of photons already emitted by
   * the currently active discrete source.
   * @param discrete_active_source_index Index of the currently active discrete
   * source.
   * @param continuous_active_number_of_photons Number of photons already
   * emitted by the continuous source.
   */
  inline PhotonSourceIndex(unsigned int discrete_active_number_of_photons,
                           unsigned int discrete_active_photon_index,
                           unsigned int discrete_active_source_index,
                           unsigned int continuous_active_number_of_photons)
      : _discrete_active_number_of_photons(discrete_active_number_of_photons),
        _discrete_active_photon_index(discrete_active_photon_index),
        _discrete_active_source_index(discrete_active_source_index),
        _continuous_active_number_of_photons(
            continuous_active_number_of_photons) {}

  /**
   * @brief Get the number of photons already emitted by the continuous source.
   *
   * @return Number of photons already emitted by the continuous source.
   */
  inline unsigned int get_continuous_active_number_of_photons() {
    return _continuous_active_number_of_photons;
  }

  /**
   * @brief Increment the number of photons already emitted by the continuous
   * source.
   */
  inline void increment_continuous_active_number_of_photons() {
    ++_continuous_active_number_of_photons;
  }

  /**
   * @brief Reset the number of photons already emitted by the continuous
   * source.
   */
  inline void reset_continuous_active_number_of_photons() {
    _continuous_active_number_of_photons = 0;
  }

  /**
   * @brief Get the index of the currently active discrete source.
   *
   * @return Index of the currently active discrete source.
   */
  inline unsigned int get_discrete_active_source_index() {
    return _discrete_active_source_index;
  }

  /**
   * @brief Increment the index of the currently active discrete source.
   */
  inline void increment_discrete_active_source_index() {
    ++_discrete_active_source_index;
  }

  /**
   * @brief Reset the index of the currently active discrete source.
   */
  inline void reset_discrete_active_source_index() {
    _discrete_active_source_index = 0;
  }

  /**
   * @brief Get the number of photons already emitted by the currently active
   * discrete source.
   *
   * @return Number of photons already emitted by the currently active discrete
   * source.
   */
  inline unsigned int get_discrete_active_photon_index() {
    return _discrete_active_photon_index;
  }

  /**
   * @brief Increment the number of photons already emitted by the currently
   * active discrete source.
   */
  inline void increment_discrete_active_photon_index() {
    ++_discrete_active_photon_index;
  }

  /**
   * @brief Reset the number of photons already emitted by the currently active
   * discrete source.
   */
  inline void reset_discrete_active_photon_index() {
    _discrete_active_photon_index = 0;
  }

  /**
   * @brief Get the total number of photons that should be emitted by the
   * currently active discrete source.
   *
   * @return Total number of photons to be emitted by the currently active
   * discrete source.
   */
  inline unsigned int get_discrete_active_number_of_photons() {
    return _discrete_active_number_of_photons;
  }

  /**
   * @brief Set the total number of photons that should be emitted by the
   * currently active discrete source.
   *
   * @param discrete_active_number_of_photons New number of photons to be
   * emitted by the currently active discrete source.
   */
  inline void set_discrete_active_number_of_photons(
      unsigned int discrete_active_number_of_photons) {
    _discrete_active_number_of_photons = discrete_active_number_of_photons;
  }
};

#endif // PHOTONSOURCEINDEX_HPP
