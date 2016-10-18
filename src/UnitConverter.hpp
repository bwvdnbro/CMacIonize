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

#include <string>

/**
 * @brief List of supported quantities.
 *
 * These are not only the basic SI quantities, but also all derived quantities
 * that are used in the program.
 */
enum Quantity {
  QUANTITY_LENGTH = 0,
  QUANTITY_MASS,
  QUANTITY_TIME,
  QUANTITY_TEMPERATURE,
  QUANTITY_NUMBER_DENSITY
};

/**
 * @brief Class that can be used to convert values from one unit system to
 * another.
 */
template < Quantity quantity > class UnitConverter {
public:
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
  static inline double to_SI(double value, std::string unit);

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
  static inline double to_unit(double value, std::string unit);

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
    return to_unit(to_SI(value, unit_from), unit_to);
  }
};

/// QUANTITY_LENGTH

/**
 * @brief UnitConverter::to_SI() specialization for a length.
 *
 * @param value Length value in strange units.
 * @param unit Strange units.
 * @return Length value in metres.
 */
template <>
inline double UnitConverter< QUANTITY_LENGTH >::to_SI(double value,
                                                      std::string unit) {
  if (unit == "m") {
    // quantity is already in SI units
    return value;
  } else if (unit == "pc") {
    return 3.086e16 * value;
  } else {
    error("Unknown length unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/**
 * @brief UnitConverter::to_unit() specialization for a length.
 *
 * @param value Length value in metres.
 * @param unit Strange units.
 * @return Length value in strange units.
 */
template <>
inline double UnitConverter< QUANTITY_LENGTH >::to_unit(double value,
                                                        std::string unit) {
  if (unit == "m") {
    // quantity is already in requested units
    return value;
  } else if (unit == "pc") {
    return value / 3.086e16;
  } else {
    error("Unknown length unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/// QUANTITY_TEMPERATURE

/**
 * @brief UnitConverter::to_SI() specialization for a temperature.
 *
 * @param value Temperature value in strange units.
 * @param unit Strange units.
 * @return Temperature value in Kelvin.
 */
template <>
inline double UnitConverter< QUANTITY_TEMPERATURE >::to_SI(double value,
                                                           std::string unit) {
  if (unit == "K") {
    // quantity is already in SI units
    return value;
  } else {
    error("Unknown temperature unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/**
 * @brief UnitConverter::to_unit() specialization for a temperature.
 *
 * @param value Temperature value in metres.
 * @param unit Strange units.
 * @return Temperature value in strange units.
 */
template <>
inline double UnitConverter< QUANTITY_TEMPERATURE >::to_unit(double value,
                                                             std::string unit) {
  if (unit == "K") {
    // quantity is already in requested units
    return value;
  } else {
    error("Unknown temperature unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/// QUANTITY_NUMBER_DENSITY

/**
 * @brief UnitConverter::to_SI() specialization for a number density.
 *
 * @param value Number density value in strange units.
 * @param unit Strange units.
 * @return Number density value in Kelvin.
 */
template <>
inline double
UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(double value,
                                                std::string unit) {
  if (unit == "m^-3") {
    // quantity is already in SI units
    return value;
  } else if (unit == "cm^-3") {
    return value * 1.e-6;
  } else {
    error("Unknown number density unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/**
 * @brief UnitConverter::to_unit() specialization for a number density.
 *
 * @param value Number density value in metres.
 * @param unit Strange units.
 * @return Number density value in strange units.
 */
template <>
inline double
UnitConverter< QUANTITY_NUMBER_DENSITY >::to_unit(double value,
                                                  std::string unit) {
  if (unit == "m^-3") {
    // quantity is already in requested units
    return value;
  } else if (unit == "cm^-3") {
    return value * 1.e6;
  } else {
    error("Unknown number density unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

#endif // UNITCONVERTER_HPP
