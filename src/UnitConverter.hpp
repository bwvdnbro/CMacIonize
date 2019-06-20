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
 * @file UnitConverter.hpp
 *
 * @brief Class that can be used to convert values from one unit system to
 * another.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UNITCONVERTER_HPP
#define UNITCONVERTER_HPP

#include "Error.hpp"
#include "PhysicalConstants.hpp"
#include "Unit.hpp"

#include <cmath>
#include <string>
#include <vector>

/// To add a new unit (such as "cm" or "kK"...):
///   - add an entry in get_single_unit() below
/// To add a new quantity (such as QUANTITY_TEMPERATURE_CHANGE...):
///   - add an entry in the Quantity enum, AND
///   - add an entry in get_SI_unit() that gives the SI unit of the new quantity
///     (e.g. K s^-1)
/// To add a new quantity conversion (such as QUANTITY_LENGTH to
/// QUANTITY_FREQUENCY for photons):
///   - add an entry in the designated part of try_conversion()
///
/// Also try adding checks to testUnitConverter. Especially when adding new
/// units or new quantity conversions.
///
/// You should not normally change anything else, unless you really know what
/// you're doing.

/**
 * @brief List of supported quantities.
 *
 * These are not only the basic SI quantities, but also all derived quantities
 * that are used in the program.
 *
 * For convenience, they are ordered alphabetically.
 */
enum Quantity {
  QUANTITY_ACCELERATION,
  QUANTITY_ANGLE,
  QUANTITY_DENSITY,
  QUANTITY_ENERGY,
  QUANTITY_ENERGY_CHANGE_RATE,
  QUANTITY_ENERGY_RATE,
  QUANTITY_FLUX,
  QUANTITY_FREQUENCY,
  QUANTITY_FREQUENCY_PER_MASS,
  QUANTITY_LENGTH,
  QUANTITY_MASS,
  QUANTITY_MASS_RATE,
  QUANTITY_NUMBER_DENSITY,
  QUANTITY_OPACITY,
  QUANTITY_REACTION_RATE,
  QUANTITY_SURFACE_AREA,
  QUANTITY_TEMPERATURE,
  QUANTITY_TIME,
  QUANTITY_VELOCITY
};

/**
 * @brief Class that can be used to convert values from one unit system to
 * another.
 */
class UnitConverter {
public:
  /**
   * @brief Get the Unit corresponding to a single, non composed unit string.
   *
   * @param name std::string containing the unique name of a unit.
   * @return Unit corresponding to that name.
   */
  static inline Unit get_single_unit(std::string name) {
    /// length units
    if (name == "m") {
      return Unit(1., 1, 0, 0, 0, 0, 0);
    } else if (name == "cm") {
      return Unit(0.01, 1, 0, 0, 0, 0, 0);
    } else if (name == "pc") {
      return Unit(3.086e16, 1, 0, 0, 0, 0, 0);
    } else if (name == "kpc") {
      return Unit(3.086e19, 1, 0, 0, 0, 0, 0);
    } else if (name == "angstrom") {
      return Unit(1.e-10, 1, 0, 0, 0, 0, 0);
    } else if (name == "km") {
      return Unit(1000., 1, 0, 0, 0, 0, 0);
      /// time units
    } else if (name == "s") {
      return Unit(1., 0, 1, 0, 0, 0, 0);
    } else if (name == "Gyr") {
      return Unit(3.154e16, 0, 1, 0, 0, 0, 0);
    } else if (name == "Myr") {
      return Unit(3.154e13, 0, 1, 0, 0, 0, 0);
    } else if (name == "yr") {
      return Unit(3.154e7, 0, 1, 0, 0, 0, 0);
      /// mass units
    } else if (name == "kg") {
      return Unit(1., 0, 0, 1, 0, 0, 0);
    } else if (name == "g") {
      return Unit(0.001, 0, 0, 1, 0, 0, 0);
    } else if (name == "Msol") {
      return Unit(1.98855e30, 0, 0, 1, 0, 0, 0);
      /// temperature units
    } else if (name == "K") {
      return Unit(1., 0, 0, 0, 1, 0, 0);
      /// angle units
    } else if (name == "radians") {
      return Unit(1., 0, 0, 0, 0, 0, 1);
    } else if (name == "degrees") {
      return Unit(M_PI / 180., 0, 0, 0, 0, 0, 1);
      /// alias units
      /// frequency units
    } else if (name == "Hz") {
      return Unit(1., 0, -1, 0, 0, 0, 0);
      /// multi-quantity units
      /// energy units
    } else if (name == "J") {
      return Unit(1., 2, -2, 1, 0, 0, 0);
    } else if (name == "erg") {
      return Unit(1.e-7, 2, -2, 1, 0, 0, 0);
    } else if (name == "eV") {
      return Unit(PhysicalConstants::get_physical_constant(
                      PHYSICALCONSTANT_ELECTRONVOLT),
                  2, -2, 1, 0, 0, 0);
    } else {
      /// error handler
      cmac_error("Unknown unit: \"%s\"!", name.c_str());
      return Unit(0., 0, 0, 0, 0, 0, 0);
    }
  }

  /**
   * @brief Get the name of the SI unit of the given quantity.
   *
   * @param quantity Quantity.
   * @return Name of the SI unit corresponding to the given quantity.
   */
  static inline std::string get_SI_unit_name(Quantity quantity) {
    switch (quantity) {
    case QUANTITY_ACCELERATION:
      return "m s^-2";
    case QUANTITY_ANGLE:
      return "radians";
    case QUANTITY_DENSITY:
      return "kg m^-3";
    case QUANTITY_ENERGY:
      return "J";
    case QUANTITY_ENERGY_CHANGE_RATE:
      return "J m^-3 s^-1";
    case QUANTITY_ENERGY_RATE:
      return "J s^-1";
    case QUANTITY_FLUX:
      return "m^-2 s^-1";
    case QUANTITY_FREQUENCY:
      return "Hz";
    case QUANTITY_FREQUENCY_PER_MASS:
      return "Hz kg^-1";
    case QUANTITY_LENGTH:
      return "m";
    case QUANTITY_MASS:
      return "kg";
    case QUANTITY_MASS_RATE:
      return "kg s^-1";
    case QUANTITY_NUMBER_DENSITY:
      return "m^-3";
    case QUANTITY_OPACITY:
      return "m^-1";
    case QUANTITY_REACTION_RATE:
      return "m^3 s^-1";
    case QUANTITY_SURFACE_AREA:
      return "m^2";
    case QUANTITY_TEMPERATURE:
      return "K";
    case QUANTITY_TIME:
      return "s";
    case QUANTITY_VELOCITY:
      return "m s^-1";
    default:
      cmac_error("Unknown quantity: %i!", quantity);
      return "";
    }
  }

  /**
   * @brief Get the SI unit of the given quantity.
   *
   * @param quantity Quantity.
   * @return SI Unit corresponding to the given quantity.
   */
  static inline Unit get_SI_unit(Quantity quantity) {
    return get_unit(get_SI_unit_name(quantity));
  }

  /**
   * @brief Try to convert a value in the given units to other given units,
   * assuming the basic quantities of both units are not the same.
   *
   * This only works for some quantities, such as e.g. energies and frequencies.
   *
   * @param value Value.
   * @param unit_from Unit of the value.
   * @param unit_to Unit to convert to.
   * @return Value in the new Unit.
   */
  static inline double try_conversion(double value, Unit unit_from,
                                      Unit unit_to) {
    std::vector< Quantity > Aunits;
    std::vector< Quantity > Bunits;
    std::vector< double > A_in_B_fac;
    std::vector< int_fast32_t > A_in_B_pow;

    /// add new quantity conversions below
    // energy to frequency conversion for photons
    Aunits.push_back(QUANTITY_ENERGY);
    Bunits.push_back(QUANTITY_FREQUENCY);
    A_in_B_fac.push_back(
        1. / PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK));
    A_in_B_pow.push_back(1);
    // wavelength to frequency conversion for photons
    Aunits.push_back(QUANTITY_LENGTH);
    Bunits.push_back(QUANTITY_FREQUENCY);
    A_in_B_fac.push_back(
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED));
    A_in_B_pow.push_back(-1);

    /// don't change the part below unless you know what you're doing
    // just try every unit combination in the lists
    for (size_t i = 0; i < Aunits.size(); ++i) {
      const Unit Aunit = get_SI_unit(Aunits[i]);
      const Unit Bunit = get_SI_unit(Bunits[i]);
      if (unit_from.is_same_quantity(Aunit) &&
          unit_to.is_same_quantity(Bunit)) {
        const double fval = 1. * unit_from;
        const double tval = 1. * unit_to;
        double Sval = value * fval;
        if (A_in_B_pow[i] > 0) {
          int_fast32_t j = 1;
          const double Sfac = Sval;
          while (j < A_in_B_pow[i]) {
            ++j;
            Sval *= Sfac;
          }
        } else {
          const double Sfac = Sval;
          Sval = 1.;
          int j = 0;
          while (j < std::abs(A_in_B_pow[i])) {
            ++j;
            Sval /= Sfac;
          }
        }
        return Sval * A_in_B_fac[i] / tval;
      }
      if (unit_from.is_same_quantity(Bunit) &&
          unit_to.is_same_quantity(Aunit)) {
        const double fval = 1. * unit_from;
        const double tval = 1. * unit_to;
        double Sval = value * fval / A_in_B_fac[i];
        const double A_in_B_pow_inv = 1. / A_in_B_pow[i];
        if (A_in_B_pow_inv > 0) {
          if (A_in_B_pow_inv == 1.) {
            int_fast32_t j = 1;
            const double Sfac = Sval;
            while (j < A_in_B_pow_inv) {
              ++j;
              Sval *= Sfac;
            }
          } else {
            Sval = std::pow(Sval, A_in_B_pow_inv);
          }
        } else {
          if (A_in_B_pow_inv == -1.) {
            const double Sfac = Sval;
            Sval = 1.;
            int_fast32_t j = 0;
            while (j < std::abs(A_in_B_pow_inv)) {
              ++j;
              Sval /= Sfac;
            }
          } else {
            Sval = std::pow(Sval, A_in_B_pow_inv);
          }
        }
        return Sval / tval;
      }
    }

    // if we did not find any matches, we throw an error
    cmac_error("No known conversion from \"%s\" to \"%s\"!",
               unit_from.to_string().c_str(), unit_to.to_string().c_str());
    return 0.;
  }

  /**
   * @brief Get the Unit corresponding to the given string.
   *
   * The string should consist of valid unit names (present in
   * get_single_unit()), and powers of these. Units can be in arbitrary order.
   * Units without a power need to be separated from other units by at least 1
   * space. Examples of valid unit names (these are tested in the unit test):
   *  - " m"
   *  - "K kg^3 s^-1m "
   *
   * Example of an invalid unit name:
   *  - "kgs": there should be at least one space in between "kg" and "s"
   *
   * @param name std::string representation of a Unit.
   * @return Unit corresponding to that name.
   */
  static inline Unit get_unit(std::string name) {
    // find the first unit name
    size_t pos1 = 0;
    while (pos1 < name.size() && !isalpha(name[pos1])) {
      ++pos1;
    }
    if (pos1 == name.size()) {
      cmac_error("Empty unit provided!");
    }
    // find either a space or a '^' or the end of the string
    size_t pos2 = pos1 + 1;
    while (pos2 < name.size() && (name[pos2] != ' ' && name[pos2] != '^')) {
      ++pos2;
    }
    Unit unit = get_single_unit(name.substr(pos1, pos2 - pos1));
    if (pos2 == name.size()) {
      return unit;
    }
    if (name[pos2] == '^') {
      ++pos2;
      pos1 = pos2;
      ++pos2;
      while (pos2 < name.size() &&
             (isdigit(name[pos2]) || name[pos2] == '+' || name[pos2] == '-')) {
        ++pos2;
      }
      const std::string powstr = name.substr(pos1, pos2 - pos1);
      const int_fast32_t power = strtod(powstr.c_str(), nullptr);
      unit ^= power;
      if (pos2 == name.size()) {
        return unit;
      }
    }
    while (pos2 < name.size()) {
      pos1 = pos2;
      while (pos1 < name.size() && !isalpha(name[pos1])) {
        ++pos1;
      }
      if (pos1 == name.size()) {
        // it could happen that there are spaces at the end of the string
        return unit;
      }
      // find either a space or a '^' or the end of the string
      pos2 = pos1 + 1;
      while (pos2 < name.size() && (name[pos2] != ' ' && name[pos2] != '^')) {
        ++pos2;
      }
      Unit unit2 = get_single_unit(name.substr(pos1, pos2 - pos1));
      if (pos2 == name.size()) {
        unit *= unit2;
        return unit;
      }
      if (name[pos2] == '^') {
        ++pos2;
        pos1 = pos2;
        ++pos2;
        while (pos2 < name.size() && (isdigit(name[pos2]) ||
                                      name[pos2] == '+' || name[pos2] == '-')) {
          ++pos2;
        }
        const std::string powstr = name.substr(pos1, pos2 - pos1);
        const int_fast32_t power = strtod(powstr.c_str(), nullptr);
        unit2 ^= power;
      }
      unit *= unit2;
    }
    return unit;
  }

  /**
   * @brief Convert the given value with the given unit to the equivalent value
   * in SI units.
   *
   * This function needs to be specialized for all quantities.
   *
   * @param value Value in obscure non SI unit.
   * @param unit Name of the obscure non SI unit.
   * @return Value in generally accepted SI units.
   */
  template < Quantity _quantity_ >
  static inline double to_SI(double value, std::string unit) {
    const Unit SI_unit = get_SI_unit(_quantity_);
    Unit strange_unit = get_unit(unit);
    if (!SI_unit.is_same_quantity(strange_unit)) {
      return try_conversion(value, strange_unit, SI_unit);
    }
    return strange_unit * value;
  }

  /**
   * @brief Convert the given value in SI units to the equivalent value in the
   * given non SI unit.
   *
   * This function needs to be specialized for all quantities.
   *
   * @param value Value in SI units.
   * @param unit Name of an obscure non SI unit.
   * @return Value in that obscure non SI unit.
   */
  template < Quantity _quantity_ >
  static inline double to_unit(double value, std::string unit) {
    const Unit SI_unit = get_SI_unit(_quantity_);
    Unit strange_unit = get_unit(unit);
    if (!SI_unit.is_same_quantity(strange_unit)) {
      return try_conversion(value, SI_unit, strange_unit);
    }
    return value / strange_unit;
  }

  /**
   * @brief Convert the given value from a given unit into another given unit.
   *
   * @param value Value to convert.
   * @param unit_from Unit the value currently has.
   * @param unit_to Unit the value should have.
   * @return Value in new unit.
   */
  static inline double convert(double value, std::string unit_from,
                               std::string unit_to) {
    Unit strange_unit_from = get_unit(unit_from);
    const Unit strange_unit_to = get_unit(unit_to);
    if (!strange_unit_to.is_same_quantity(strange_unit_from)) {
      return try_conversion(value, strange_unit_from, strange_unit_to);
      cmac_error("Units are not compatible: \"%s\" and \"%s\"!",
                 strange_unit_from.to_string().c_str(),
                 strange_unit_to.to_string().c_str());
    }
    strange_unit_from /= strange_unit_to;
    return value * strange_unit_from;
  }
};

#endif // UNITCONVERTER_HPP
