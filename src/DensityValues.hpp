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
 * @file DensityValues.hpp
 *
 * @brief Density values at a point in space.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYVALUES_HPP
#define DENSITYVALUES_HPP

#include "CoordinateVector.hpp"
#include "ElementNames.hpp"

/**
 * @brief Density values at a point in space.
 */
class DensityValues {
private:
  /*! @brief Number density of hydrogen (in m^-3). */
  double _number_density;

  /*! @brief Ionic fractions. For hydrogen and helium, these are the neutral
   *  fractions. For other elements, they are the fraction of the end product
   *  of ionization (e.g. _ionic_fraction[ION_C_p1] is the fraction of C that
   *  is in the form of C++). */
  double _ionic_fraction[NUMBER_OF_IONNAMES];

  /*! @brief Temperature (in K). */
  double _temperature;

  /*! @brief Fluid velocity (in m s^-1). */
  CoordinateVector<> _velocity;

  /*! @brief Magnetic field strength (in kg A^-1 s^-2). */
  CoordinateVector<> _magnetic_field;

  /*! @brief Cosmic ray energy (in m^2 s^-2). */
  double _cosmic_ray_energy;

  /*! @brief Cosmic ray heating factor (in kg m A^-1 s^-4). */
  double _cosmic_ray_factor;

public:
  /**
   * @brief Empty constructor.
   */
  inline DensityValues()
      : _number_density(0.), _temperature(0.), _cosmic_ray_energy(0.),
        _cosmic_ray_factor(-1.) {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _ionic_fraction[i] = 0.;
    }
  }

  /**
   * @brief Set the number density of hydrogen.
   *
   * @param number_density Value for the number density (in m^-3).
   */
  inline void set_number_density(double number_density) {
    _number_density = number_density;
  }

  /**
   * @brief Set the ionic fraction of the given ion.
   *
   * @param ion IonName of a valid ion.
   * @param ionic_fraction New value for the ionic fraction.
   */
  inline void set_ionic_fraction(IonName ion, double ionic_fraction) {
    _ionic_fraction[ion] = ionic_fraction;
  }

  /**
   * @brief Set the temperature.
   *
   * @param temperature Temperature value (in K).
   */
  inline void set_temperature(double temperature) {
    _temperature = temperature;
  }

  /**
   * @brief Set the fluid velocity.
   *
   * @param velocity Fluid velocity (in m s^-1).
   */
  inline void set_velocity(CoordinateVector<> velocity) {
    _velocity = velocity;
  }

  /**
   * @brief Set the magnetic field.
   *
   * @param magnetic_field Magnetic field (in kg A^-1 s^-2).
   */
  inline void set_magnetic_field(CoordinateVector<> magnetic_field) {
    _magnetic_field = magnetic_field;
  }

  /**
   * @brief Set the cosmic ray energy.
   *
   * @param cosmic_ray_energy Cosmic ray energy (in m^2 s^-2).
   */
  inline void set_cosmic_ray_energy(double cosmic_ray_energy) {
    _cosmic_ray_energy = cosmic_ray_energy;
  }

  /**
   * @brief Set the cosmic ray factor.
   *
   * @param cosmic_ray_factor Cosmic ray factor (in kg m A^-1 s^-4).
   */
  inline void set_cosmic_ray_factor(double cosmic_ray_factor) {
    _cosmic_ray_factor = cosmic_ray_factor;
  }

  /**
   * @brief Get the number density of hydrogen.
   *
   * @return Number density (in m^-3).
   */
  inline double get_number_density() const { return _number_density; }

  /**
   * @brief Get the ionic fraction of the given ion.
   *
   * @param ion IonName of a valid ion.
   * @return Ionic fraction.
   */
  inline double get_ionic_fraction(IonName ion) const {
    return _ionic_fraction[ion];
  }

  /**
   * @brief Get the temperature.
   *
   * @return Temperature (in K).
   */
  inline double get_temperature() const { return _temperature; }

  /**
   * @brief Get the fluid velocity.
   *
   * @return Fluid velocity (in m s^-1).
   */
  inline CoordinateVector<> get_velocity() const { return _velocity; }

  /**
   * @brief Get the magnetic field.
   *
   * @return Magnetic field (in kg A^-1 s^-2).
   */
  inline CoordinateVector<> get_magnetic_field() const {
    return _magnetic_field;
  }

  /**
   * @brief Get the cosmic ray energy.
   *
   * @return Cosmic ray energy (in m^2 s^-2).
   */
  inline double get_cosmic_ray_energy() const { return _cosmic_ray_energy; }

  /**
   * @brief Get the cosmic ray factor.
   *
   * @return Cosmic ray factor (in kg m A^-1 s^-4).
   */
  inline double get_cosmic_ray_factor() const { return _cosmic_ray_factor; }
};

#endif // DENSITYVALUES_HPP
