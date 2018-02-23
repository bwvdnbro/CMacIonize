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
 * @file PhysicalConstants.hpp
 *
 * @brief Values for physical constants in SI units.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHYSICALCONSTANTS_HPP
#define PHYSICALCONSTANTS_HPP

#include "Error.hpp"

/**
 * @brief Names for physical constants.
 */
enum PhysicalConstantName {
  /*! @brief Planck constant \f$h\f$ (in J s). */
  PHYSICALCONSTANT_PLANCK = 0,
  /*! @brief Boltzmann constant \f$k\f$ (in J s^-1). */
  PHYSICALCONSTANT_BOLTZMANN,
  /*! @brief Speed of light \f$c\f$ (in m s^-1). */
  PHYSICALCONSTANT_LIGHTSPEED,
  /*! @brief Electron volt eV (in J). */
  PHYSICALCONSTANT_ELECTRONVOLT,
  /*! @brief Mass of a proton (in kg). */
  PHYSICALCONSTANT_PROTON_MASS,
  /*! @brief Mass of an electron (in kg). */
  PHYSICALCONSTANT_ELECTRON_MASS,
  /*! @brief Rydberg energy (in J). */
  PHYSICALCONSTANT_RYDBERG_ENERGY,
  /*! @brief Newton gravitational constant (in m^3 kg^-1 s^-2). */
  PHYSICALCONSTANT_NEWTON_CONSTANT,
  /*! @brief Solar mass (in kg). */
  PHYSICALCONSTANT_SOLAR_MASS,
  /*! @brief Astronomical unit (in m). */
  PHYSICALCONSTANT_ASTRONOMICAL_UNIT,
  /*! @brief Counter (should always be last element!). */
  NUMBER_OF_PHYSICALCONSTANTS
};

/**
 * @brief Values for physical constants in SI units.
 */
class PhysicalConstants {
public:
  /**
   * @brief Get the value of the given physical constant.
   *
   * @param name PhysicalConstantName.
   * @return Value of that constant (in SI units).
   */
  static inline double get_physical_constant(PhysicalConstantName name) {

    switch (name) {

    case PHYSICALCONSTANT_PLANCK:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?h
      // in J s
      return 6.626070040e-34;

    case PHYSICALCONSTANT_BOLTZMANN:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?k
      // in J s^-1
      return 1.38064852e-23;

    case PHYSICALCONSTANT_LIGHTSPEED:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?c
      // in m s^-1
      return 299792458.;

    case PHYSICALCONSTANT_ELECTRONVOLT:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?tevj
      // in J
      return 1.6021766208e-19;

    case PHYSICALCONSTANT_PROTON_MASS:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?mp
      // in kg
      return 1.672621898e-27;

    case PHYSICALCONSTANT_ELECTRON_MASS:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?me
      // in kg
      return 9.10938356e-31;

    case PHYSICALCONSTANT_RYDBERG_ENERGY:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?rydhcj
      // in J
      return 2.179872325e-18;

    case PHYSICALCONSTANT_NEWTON_CONSTANT:
      // NIST 2014 CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?bg
      // in m^3 kg^-1 s^-2
      return 6.67408e-11;

    case PHYSICALCONSTANT_SOLAR_MASS:
      // http://asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
      // in kg
      return 1.9884e30;

    case PHYSICALCONSTANT_ASTRONOMICAL_UNIT:
      // http://asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
      // in m
      return 149597870700.;

    default:
      cmac_error("Unknown physical constant: %i!", name);
      return 0.;
    }
  }
};

#endif // PHYSICALCONSTANTS_HPP
