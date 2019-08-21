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
 * @file PhotonPacket.hpp
 *
 * @brief Photon packet.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONPACKET_HPP
#define PHOTONPACKET_HPP

/*! @brief Size of the MPI buffer necessary to store a single Photon. */
#define PHOTON_MPI_SIZE (9 * sizeof(double))

#include "Configuration.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/**
 * @brief Photon packet.
 */
class PhotonPacket {
private:
  /*! @brief Current position of the photon packet (in m). */
  double _position[3];

  /*! @brief Propagation direction of the photon packet. */
  double _direction[3];

  /*! @brief Target optical depth for the photon packet. */
  double _target_optical_depth;

  /*! @brief Photoionization cross section of the photons in the photon packet
   *  (in m^2). */
  double _photoionization_cross_section;

  /*! @brief Weight of the photon packet. */
  double _weight;

public:
#ifdef HAVE_MPI
  /**
   * @brief Store the contents of the PhotonPacket in the given MPI
   * communication buffer.
   *
   * @param buffer Buffer to use (should be preallocated and have at least size
   * PHOTON_MPI_SIZE).
   */
  inline void pack(char buffer[PHOTON_MPI_SIZE]) {
    int_least32_t buffer_position = 0;
    MPI_Pack(_position, 3, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(_direction, 3, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_target_optical_depth, 1, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_photoionization_cross_section, 1, MPI_DOUBLE, buffer,
             PHOTON_MPI_SIZE, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&_weight, 1, MPI_DOUBLE, buffer, PHOTON_MPI_SIZE, &buffer_position,
             MPI_COMM_WORLD);
  }

  /**
   * @brief Copy the contents of the given MPI communication buffer into this
   * PhotonPacket.
   *
   * @param buffer MPI communication buffer (should have at least size
   * PHOTON_MPI_SIZE).
   */
  inline void unpack(char buffer[PHOTON_MPI_SIZE]) {
    int_least32_t buffer_position = 0;
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position, _position, 3,
               MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position, _direction, 3,
               MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position,
               &_target_optical_depth, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position,
               &_photoionization_cross_section, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, PHOTON_MPI_SIZE, &buffer_position, &_weight, 1,
               MPI_DOUBLE, MPI_COMM_WORLD);
  }
#endif

  /**
   * @brief Check that the given PhotonPacket is equal to this one.
   *
   * @param other Other PhotonPacket.
   */
  inline void check_equal(const PhotonPacket &other) {

    cmac_assert_message(_position[0] == other._position[0],
                        "Positions do not match!");
    cmac_assert_message(_position[1] == other._position[1],
                        "Positions do not match!");
    cmac_assert_message(_position[2] == other._position[2],
                        "Positions do not match!");
    cmac_assert_message(_direction[0] == other._direction[0],
                        "Directions do not match!");
    cmac_assert_message(_direction[1] == other._direction[1],
                        "Directions do not match!");
    cmac_assert_message(_direction[2] == other._direction[2],
                        "Directions do not match!");
    cmac_assert_message(_target_optical_depth == other._target_optical_depth,
                        "Target optical depths do not match!");
    cmac_assert_message(_photoionization_cross_section ==
                            other._photoionization_cross_section,
                        "Cross sections do not match!");
    cmac_assert_message(_weight == other._weight, "Weights do not match!");
  }

  /**
   * @brief Get a constant pointer to the direction array.
   *
   * @return Constant pointer to the direction array.
   */
  inline const double *get_direction() const { return _direction; }

  /**
   * @brief Get a writable pointer to the direction array.
   *
   * @return Writable pointer to the direction array.
   */
  inline double *get_direction() { return _direction; }

  /**
   * @brief Set the direction vector directly.
   *
   * @param x x component.
   * @param y y component.
   * @param z z component.
   */
  inline void set_direction(const double x, const double y, const double z) {
    _direction[0] = x;
    _direction[1] = y;
    _direction[2] = z;
  }

  /**
   * @brief Get a constant pointer to the position array.
   *
   * @return Constant pointer to the position array.
   */
  inline const double *get_position() const { return _position; }

  /**
   * @brief Update the position of the photon packet.
   *
   * @param x New x coordinate (in m).
   * @param y New y coordinate (in m).
   * @param z New z coordinate (in m).
   */
  inline void set_position(const double x, const double y, const double z) {
    _position[0] = x;
    _position[1] = y;
    _position[2] = z;
  }

  /**
   * @brief Get the target optical depth for the photon packet.
   *
   * @return Target optical depth.
   */
  inline double get_target_optical_depth() const {
    return _target_optical_depth;
  }

  /**
   * @brief Set the target optical depth for the photon packet.
   *
   * @param target_optical_depth New target optical depth.
   */
  inline void set_target_optical_depth(const double target_optical_depth) {
    _target_optical_depth = target_optical_depth;
  }

  /**
   * @brief Get the photoionization cross section for the photon packet.
   *
   * @return Photoionization cross section for the photon packet (in m^2).
   */
  inline double get_photoionization_cross_section() const {
    return _photoionization_cross_section;
  }

  /**
   * @brief Set the photoionization cross section for the photon packet.
   *
   * @param photoionization_cross_section New photoionization cross section
   * (in m^2).
   */
  inline void set_photoionization_cross_section(
      const double photoionization_cross_section) {
    _photoionization_cross_section = photoionization_cross_section;
  }

  /**
   * @brief Get the weight for the photon packet.
   *
   * @return Weight for the photon packet.
   */
  inline double get_weight() const { return _weight; }

  /**
   * @brief Set the weight for the photon packet.
   *
   * @param weight New weight for the photon packet.
   */
  inline void set_weight(const double weight) { _weight = weight; }
};

#endif // PHOTONPACKET_HPP
