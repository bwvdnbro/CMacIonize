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
 * @file Photon.hpp
 *
 * @brief Photon package.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTON_HPP
#define PHOTON_HPP

#include "CoordinateVector.hpp"

/**
 * @brief Photon package.
 */
class Photon {
private:
  /*! @brief Current position of the photon. */
  CoordinateVector<> _position;

  /*! @brief Current direction the photon is moving in. */
  CoordinateVector<> _direction;

  /*! @brief Current energy contents of the photon. */
  double _energy;

  /*! @brief Hydrogen ionization cross section. */
  double _xsecH;

  /*! @brief Helium ionization cross section. */
  double _xsecHe;

public:
  /**
   * @brief Constructor.
   *
   * @param position Initial position of the photon.
   * @param direction Initial direction of the photon.
   * @param energy Initial energy of the photon.
   * @param xsecH Hydrogen photoionization cross section of the photon.
   * @param xsecHe Helium photoionization cross section of the photon.
   */
  inline Photon(CoordinateVector<> position, CoordinateVector<> direction,
                double energy, double xsecH, double xsecHe)
      : _position(position), _direction(direction), _energy(energy),
        _xsecH(xsecH), _xsecHe(xsecHe) {}

  /**
   * @brief Get the current position of the photon.
   *
   * @return Current position of the photon.
   */
  inline CoordinateVector<> get_position() { return _position; }

  /**
   * @brief Get the current direction the photon is moving in.
   *
   * @return Current movement direction of the photon.
   */
  inline CoordinateVector<> get_direction() { return _direction; }

  /**
   * @brief Get the current energy of the photon.
   *
   * @return Current energy of the photon.
   */
  inline double get_energy() { return _energy; }

  /**
   * @brief Get the ionization cross section for hydrogen.
   *
   * @return Hydrogen ionization cross section.
   */
  inline double get_hydrogen_cross_section() { return _xsecH; }

  /**
   * @brief Get the ionization cross section for helium.
   *
   * @return Helium ionization cross section.
   */
  inline double get_helium_cross_section() { return _xsecHe; }

  /**
   * @brief Set the position of the photon.
   *
   * @param position New position of the photon.
   */
  inline void set_position(CoordinateVector<> position) {
    _position = position;
  }

  /**
   * @brief Set the direction of the photon.
   *
   * @param direction New direction of the photon.
   */
  inline void set_direction(CoordinateVector<> direction) {
    _direction = direction;
  }

  /**
   * @brief Set the energy of the photon.
   *
   * @param energy Energy of the photon.
   */
  inline void set_energy(double energy) { _energy = energy; }

  /**
   * @brief Set the hydrogen ionization cross section.
   *
   * @param xsecH Hydrogen ionization cross section.
   */
  inline void set_hydrogen_cross_section(double xsecH) { _xsecH = xsecH; }

  /**
   * @brief Set the helium ionization cross section.
   *
   * @param xsecHe Helium ionization cross section.
   */
  inline void set_helium_cross_section(double xsecHe) { _xsecHe = xsecHe; }
};

#endif // PHOTON_HPP
