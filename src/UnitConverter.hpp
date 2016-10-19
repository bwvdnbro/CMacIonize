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
 *
 * For convenience, they are ordered alphabetically.
 */
enum Quantity {
  QUANTITY_FREQUENCY,
  QUANTITY_LENGTH,
  QUANTITY_MASS,
  QUANTITY_NUMBER_DENSITY,
  QUANTITY_REACTION_RATE,
  QUANTITY_SURFACE_AREA,
  QUANTITY_TEMPERATURE,
  QUANTITY_TIME
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

/// QUANTITY_FREQUENCY

/**
 * @brief UnitConverter::to_SI() specialization for a frequency.
 *
 * @param value Frequency value in strange units.
 * @param unit Strange units.
 * @return Frequency value in per second (Hertz).
 */
template <>
inline double UnitConverter< QUANTITY_FREQUENCY >::to_SI(double value,
                                                         std::string unit) {
  if (unit == "s^-1") {
    // quantity is already in SI units
    return value;
  } else if (unit == "Hz") {
    // just an alias for the SI unit name
    return value;
  } else if (unit == "eV") {
    return value / 4.135668e-15;
  } else {
    error("Unknown frequency unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/**
 * @brief UnitConverter::to_unit() specialization for a frequency.
 *
 * @param value Frequency value in per second (Hertz).
 * @param unit Strange units.
 * @return Frequency value in strange units.
 */
template <>
inline double UnitConverter< QUANTITY_FREQUENCY >::to_unit(double value,
                                                           std::string unit) {
  if (unit == "s^-1") {
    // quantity is already in requested units
    return value;
  } else if (unit == "Hz") {
    return value;
  } else if (unit == "eV") {
    return value * 4.135668e-15;
  } else {
    error("Unknown frequency unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

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
  } else if (unit == "cm") {
    return 0.01 * value;
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
  } else if (unit == "cm") {
    return 100. * value;
  } else if (unit == "pc") {
    return value / 3.086e16;
  } else {
    error("Unknown length unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/// QUANTITY_MASS

/**
 * @brief UnitConverter::to_SI() specialization for a mass.
 *
 * @param value Mass value in strange units.
 * @param unit Strange units.
 * @return Mass value in kilograms.
 */
template <>
inline double UnitConverter< QUANTITY_MASS >::to_SI(double value,
                                                    std::string unit) {
  if (unit == "kg") {
    // quantity is already in SI units
    return value;
  } else if (unit == "g") {
    return 0.001 * value;
  } else {
    error("Unknown mass unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/**
 * @brief UnitConverter::to_unit() specialization for a mass.
 *
 * @param value Mass value in kilograms.
 * @param unit Strange units.
 * @return Mass value in strange units.
 */
template <>
inline double UnitConverter< QUANTITY_MASS >::to_unit(double value,
                                                      std::string unit) {
  if (unit == "kg") {
    // quantity is already in requested units
    return value;
  } else if (unit == "g") {
    return value * 1000.;
  } else {
    error("Unknown mass unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/// QUANTITY_NUMBER_DENSITY

/**
 * @brief UnitConverter::to_SI() specialization for a number density.
 *
 * @param value Number density value in strange units.
 * @param unit Strange units.
 * @return Number density value in per cubic metres.
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
 * @param value Number density value in per cubic metres.
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

/// QUANTITY_REACTION_RATE

/**
 * @brief UnitConverter::to_SI() specialization for a reaction rate.
 *
 * @param value Reaction rate value in strange units.
 * @param unit Strange units.
 * @return Reaction rate value in cubic metres per second.
 */
template <>
inline double UnitConverter< QUANTITY_REACTION_RATE >::to_SI(double value,
                                                             std::string unit) {
  if (unit == "m^3s^-1") {
    // quantity is already in SI units
    return value;
  } else if (unit == "cm^3s^-1") {
    return 1.e-6 * value;
  } else {
    error("Unknown reaction rate unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/**
 * @brief UnitConverter::to_unit() specialization for a reaction rate.
 *
 * @param value Reaction rate value in cubic metres per second.
 * @param unit Strange units.
 * @return Reaction rate value in strange units.
 */
template <>
inline double
UnitConverter< QUANTITY_REACTION_RATE >::to_unit(double value,
                                                 std::string unit) {
  if (unit == "m^3s^-1") {
    // quantity is already in requested units
    return value;
  } else if (unit == "cm^3s^-1") {
    return 1.e6 * value;
  } else {
    error("Unknown reaction rate unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/// QUANTITY_SURFACE_AREA

/**
 * @brief UnitConverter::to_SI() specialization for a surface area.
 *
 * @param value Surface area value in strange units.
 * @param unit Strange units.
 * @return Surface area value in squared metres.
 */
template <>
inline double UnitConverter< QUANTITY_SURFACE_AREA >::to_SI(double value,
                                                            std::string unit) {
  if (unit == "m^2") {
    // quantity is already in SI units
    return value;
  } else if (unit == "cm^2") {
    return 1.e-4 * value;
  } else {
    error("Unknown surface area unit: \"%s\".", unit.c_str());
    return 0.;
  }
}

/**
 * @brief UnitConverter::to_unit() specialization for a surface area.
 *
 * @param value Surface area value in squared metres.
 * @param unit Strange units.
 * @return Surface area value in strange units.
 */
template <>
inline double
UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(double value,
                                                std::string unit) {
  if (unit == "m^2") {
    // quantity is already in requested units
    return value;
  } else if (unit == "cm^2") {
    return 1.e4 * value;
  } else {
    error("Unknown surface area unit: \"%s\".", unit.c_str());
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
 * @param value Temperature value in Kelvin.
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

#endif // UNITCONVERTER_HPP
