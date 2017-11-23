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

  /*! @brief Inverse of the movement direction. */
  CoordinateVector<> _inverse_direction;

  /*! @brief Current energy contents of the photon (in Hz). */
  double _energy;

  /*! @brief Ionization cross sections (in m^2). */
  double _cross_sections[NUMBER_OF_IONNAMES];

  /*! @brief Abundance corrected helium cross section (in m^2). */
  double _cross_section_He_corr;

  /*! @brief Type of the photon. All photons start off as PHOTONTYPE_PRIMARY,
   *  but their type can change during reemission events. */
  PhotonType _type;

  /*! @brief Weight of the photon. */
  double _weight;

  /*! @brief Stokes parameters for scattering. */
  double _stokes_parameters[4];

  /*! @brief Direction parameters: @f$\sin(\theta{})@f$, @f$\cos(\theta{})@f$,
   *  @f$\phi{}@f$ (in radians), @f$\sin(\phi{})@f$, and @f$\cos(\phi{})@f$. */
  double _direction_parameters[5];

public:
  /**
   * @brief Constructor.
   *
   * @param position Initial position of the photon (in m).
   * @param direction Initial direction of the photon.
   * @param energy Initial energy of the photon (in Hz).
   */
  inline Photon(CoordinateVector<> position, CoordinateVector<> direction,
                double energy)
      : _position(position), _direction(direction),
        _inverse_direction(1. / direction.x(), 1. / direction.y(),
                           1. / direction.z()),
        _energy(energy), _cross_section_He_corr(0.), _type(PHOTONTYPE_PRIMARY),
        _weight(1.) {

    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _cross_sections[i] = 0.;
    }
    _stokes_parameters[0] = 1.;
    _stokes_parameters[1] = 0.;
    _stokes_parameters[2] = 0.;
    _stokes_parameters[3] = 0.;
  }

  /**
   * @brief Get the current position of the photon.
   *
   * @return Current position of the photon (in m).
   */
  inline CoordinateVector<> get_position() const { return _position; }

  /**
   * @brief Get the current direction the photon is moving in.
   *
   * @return Current movement direction of the photon.
   */
  inline CoordinateVector<> get_direction() const { return _direction; }

  /**
   * @brief Get the inverse of the direction the photon is currently moving in.
   *
   * @return Inverse current movement direction of the photon.
   */
  inline CoordinateVector<> get_inverse_direction() const {
    return _inverse_direction;
  }

  /**
   * @brief Get the current energy of the photon.
   *
   * @return Current energy of the photon (in Hz).
   */
  inline double get_energy() const { return _energy; }

  /**
   * @brief Get the ionization cross section for the given ion.
   *
   * @param ion IonName of a valid ion.
   * @return Ionization cross section (in m^2).
   */
  inline double get_cross_section(IonName ion) const {
    return _cross_sections[ion];
  }

  /**
   * @brief Get the abundance corrected helium cross section.
   *
   * @return Abundance corrected helium cross section (in m^2).
   */
  inline double get_cross_section_He_corr() const {
    return _cross_section_He_corr;
  }

  /**
   * @brief Get the type of the photon.
   *
   * @return PhotonType type identifier.
   */
  inline PhotonType get_type() const { return _type; }

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
    _inverse_direction[0] = 1. / direction.x();
    _inverse_direction[1] = 1. / direction.y();
    _inverse_direction[2] = 1. / direction.z();
  }

  /**
   * @brief Set the energy of the photon.
   *
   * @param energy Energy of the photon (in Hz).
   */
  inline void set_energy(double energy) { _energy = energy; }

  /**
   * @brief Set the ionization cross section for the given ion.
   *
   * @param ion IonName of a valid ion.
   * @param cross_section Ionization cross section (in m^2).
   */
  inline void set_cross_section(IonName ion, double cross_section) {
    _cross_sections[ion] = cross_section;
  }

  /**
   * @brief Set the abundance corrected helium cross section.
   *
   * @param cross_section_He_corr Abundance corrected helium cross section (in
   * m^2).
   */
  inline void set_cross_section_He_corr(double cross_section_He_corr) {
    _cross_section_He_corr = cross_section_He_corr;
  }

  /**
   * @brief Set the photon type.
   *
   * @param type PhotonType type identifier.
   */
  inline void set_type(PhotonType type) { _type = type; }

  /**
   * @brief Set the weight of the Photon.
   *
   * @param weight New weight for the Photon.
   */
  inline void set_weight(double weight) { _weight = weight; }

  /**
   * @brief Get the weight of the Photon.
   *
   * @return Weight of the Photon.
   */
  inline double get_weight() const { return _weight; }

  /**
   * @brief Get the Stokes parameters for scattering and store them in the given
   * variables.
   *
   * @param fi First parameter.
   * @param fq Second parameter.
   * @param fu Third parameter.
   * @param fv Fourth parameter.
   */
  inline void get_stokes_parameters(double &fi, double &fq, double &fu,
                                    double &fv) const {
    fi = _stokes_parameters[0];
    fq = _stokes_parameters[1];
    fu = _stokes_parameters[2];
    fv = _stokes_parameters[3];
  }

  /**
   * @brief Set new values for the Stokes parameters.
   *
   * @param fi First parameter.
   * @param fq Second parameter.
   * @param fu Third parameter.
   * @param fv Fourth parameter.
   */
  inline void set_stokes_parameters(double fi, double fq, double fu,
                                    double fv) {
    _stokes_parameters[0] = fi;
    _stokes_parameters[1] = fq;
    _stokes_parameters[2] = fu;
    _stokes_parameters[3] = fv;
  }

  /**
   * @brief Get the direction components of the photon direction.
   *
   * @param sin_theta Variable to store @f$\sin(\theta{})@f$ in.
   * @param cos_theta Variable to store @f$\cos(\theta{})@f$ in.
   * @param phi Variable to store @f$\phi@f$ in (in radians).
   * @param sin_phi Variable to store @f$\sin(\phi{})@f$ in.
   * @param cos_phi Variable to store @f$\cos(\phi{})@f$ in.
   */
  inline void get_direction_parameters(double &sin_theta, double &cos_theta,
                                       double &phi, double &sin_phi,
                                       double &cos_phi) const {
    sin_theta = _direction_parameters[0];
    cos_theta = _direction_parameters[1];
    phi = _direction_parameters[2];
    sin_phi = _direction_parameters[3];
    cos_phi = _direction_parameters[4];
  }

  /**
   * @brief Set new values for the direction components.
   *
   * @param sin_theta New value for @f$\sin(\theta{})@f$.
   * @param cos_theta New value for @f$\cos(\theta{})@f$.
   * @param phi New value for @f$\phi{}@f$ (in radians).
   * @param sin_phi New value for @f$\sin(\phi{})@f$.
   * @param cos_phi New value for @f$\cos(\phi{})@f$.
   */
  inline void set_direction_parameters(double sin_theta, double cos_theta,
                                       double phi, double sin_phi,
                                       double cos_phi) {
    _direction_parameters[0] = sin_theta;
    _direction_parameters[1] = cos_theta;
    _direction_parameters[2] = phi;
    _direction_parameters[3] = sin_phi;
    _direction_parameters[4] = cos_phi;
  }
};

#endif // PHOTON_HPP
