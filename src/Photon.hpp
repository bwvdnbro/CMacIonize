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
#include "ElementNames.hpp"

/**
 * @brief Photon types.
 *
 * All photons start as type 0 (PHOTONTYPE_PRIMARY), but during scattering
 * events their type changes.
 * By recording a type, we can check how an observer would see the photon.
 */
enum PhotonType {
  PHOTONTYPE_PRIMARY = 0,
  PHOTONTYPE_DIFFUSE_HI,
  PHOTONTYPE_DIFFUSE_HeI,
  PHOTONTYPE_ABSORBED,
  // THIS ELEMENT SHOULD ALWAYS BE LAST!
  // It is used to initialize arrays that have an entry for each PhotonType.
  // By putting it last, PHOTONTYPE_NUMBER will have an integer value equal to
  // the number of defined types above.
  PHOTONTYPE_NUMBER
};

/**
 * @brief Photon package.
 */
class Photon {
private:
  /*! @brief Current position of the photon (in m). */
  CoordinateVector<> _position;

  /*! @brief Current direction the photon is moving in. */
  CoordinateVector<> _direction;

  /*! @brief Current energy contents of the photon (in Hz). */
  double _energy;

  /*! @brief Ionization cross sections (in m^2). */
  double _cross_sections[NUMBER_OF_ELEMENTS];

  /*! @brief Type of the photon. All photons start off as PHOTONTYPE_PRIMARY,
   *  but their type can change during reemission events. */
  PhotonType _type;

public:
  /**
   * @brief Constructor.
   *
   * @param position Initial position of the photon (in m).
   * @param direction Initial direction of the photon.
   * @param energy Initial energy of the photon (in Hz).
   * @param xsecH Hydrogen photoionization cross section of the photon (in m^2).
   * @param xsecHe Helium photoionization cross section of the photon (in m^2).
   */
  inline Photon(CoordinateVector<> position, CoordinateVector<> direction,
                double energy, double xsecH, double xsecHe)
      : _position(position), _direction(direction), _energy(energy),
        _type(PHOTONTYPE_PRIMARY) {
    _cross_sections[ELEMENT_H] = xsecH;
    _cross_sections[ELEMENT_He] = xsecHe;
  }

  /**
   * @brief Get the current position of the photon.
   *
   * @return Current position of the photon (in m).
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
   * @return Current energy of the photon (in Hz).
   */
  inline double get_energy() { return _energy; }

  /**
   * @brief Get the ionization cross section for the given element.
   *
   * @param element ElementName of a valid element.
   * @return Ionization cross section (in m^2).
   */
  inline double get_cross_section(ElementName element) {
    return _cross_sections[element];
  }

  /**
   * @brief Get the type of the photon.
   *
   * @return PhotonType type identifier.
   */
  inline PhotonType get_type() { return _type; }

  /**
   * @brief Set the position of the photon.
   *
   * @param position New position of the photon (in m).
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
   * @param energy Energy of the photon (in Hz).
   */
  inline void set_energy(double energy) { _energy = energy; }

  /**
   * @brief Set the ionization cross section for the given element.
   *
   * @param element ElementName of a valid element.
   * @param cross_section Ionization cross section (in m^2).
   */
  inline void set_cross_section(ElementName element, double cross_section) {
    _cross_sections[element] = cross_section;
  }

  /**
   * @brief Set the photon type.
   *
   * @param type PhotonType type identifier.
   */
  inline void set_type(PhotonType type) { _type = type; }
};

#endif // PHOTON_HPP
