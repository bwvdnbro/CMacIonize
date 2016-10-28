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
 * @file testUnitConverter.cpp
 *
 * @brief Unit test for the UnitConverter class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "UnitConverter.hpp"

/**
 * @brief Unit test for the UnitConverter class.
 *
 * The UnitConverter class supports three methods for each quantity:
 *   - a method that converts from a strange unit to SI units
 *   - a method that converts from SI units to a strange unit
 *   - a method that converts from one strange unit into another strange unit
 * The last method uses the other two and is independent of the quantity or
 * unit (and we know it works properly). To check whether the methods are
 * correctly implemented for a given unit, we can hence simply convert a value
 * from that unit into itself, and check if the resulting value is still the
 * same. This will fail if the two other methods are inconsistent, indicating
 * that at least one of them is wrong.
 *
 * To check for both methods being wrong, we implement a simple sanity check
 * test: we get the value of the unit in SI units, and then check if this value
 * is larger or smaller than the SI unit.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  /// DENSITY
  assert_condition(UnitConverter< QUANTITY_DENSITY >::convert(1., "kg m^-3",
                                                              "kg m^-3") == 1.);

  assert_condition(UnitConverter< QUANTITY_DENSITY >::convert(1., "g cm^-3",
                                                              "g cm^-3") == 1.);
  // 1 g in a cubic centimetre is more than 1 kg in a cubic metre
  assert_condition(UnitConverter< QUANTITY_DENSITY >::to_SI(1, "g cm^-3") > 1.);

  /// FREQUENCY
  assert_condition(
      UnitConverter< QUANTITY_FREQUENCY >::convert(1., "s^-1", "s^-1") == 1.);

  assert_condition(
      UnitConverter< QUANTITY_FREQUENCY >::convert(1., "Hz", "Hz") == 1.);
  // Hz is just another name for s^-1
  assert_condition(UnitConverter< QUANTITY_FREQUENCY >::to_SI(1., "Hz") == 1.);

  assert_condition(
      UnitConverter< QUANTITY_FREQUENCY >::convert(1., "eV", "eV") == 1.);
  // 1 eV corresponds to a very high frequency in Hz
  assert_condition(UnitConverter< QUANTITY_FREQUENCY >::to_SI(1., "eV") > 1.);

  /// LENGTH
  assert_condition(UnitConverter< QUANTITY_LENGTH >::convert(1., "m", "m") ==
                   1.);

  assert_condition(UnitConverter< QUANTITY_LENGTH >::convert(1., "cm", "cm") ==
                   1.);
  // cm is smaller than m
  assert_condition(UnitConverter< QUANTITY_LENGTH >::to_SI(1., "cm") < 1.);

  assert_condition(UnitConverter< QUANTITY_LENGTH >::convert(1., "pc", "pc") ==
                   1.);
  // a parsec is a *bit* larger than a metre
  assert_condition(UnitConverter< QUANTITY_LENGTH >::to_SI(1., "pc") > 1.);

  assert_condition(
      UnitConverter< QUANTITY_LENGTH >::convert(1., "kpc", "kpc") == 1.);
  // a kpc is even larger
  assert_condition(UnitConverter< QUANTITY_LENGTH >::to_SI(1., "kpc") > 1.);

  /// MASS
  assert_condition(UnitConverter< QUANTITY_MASS >::convert(1., "kg", "kg") ==
                   1.);

  assert_condition(UnitConverter< QUANTITY_MASS >::convert(1., "g", "g") == 1.);
  // a gramme is less than a kilogram
  assert_condition(UnitConverter< QUANTITY_MASS >::to_SI(1., "g") < 1.);

  /// NUMBER DENSITY
  assert_condition(
      UnitConverter< QUANTITY_NUMBER_DENSITY >::convert(1., "m^-3", "m^-3"));

  assert_condition(
      UnitConverter< QUANTITY_NUMBER_DENSITY >::convert(1., "cm^-3", "cm^-3"));
  // 1 particle in a cubic centimetre is much more than 1 particle in a cubic
  // metre
  assert_condition(
      UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(1., "cm^-3") > 1.);

  /// REACTION RATE
  assert_condition(UnitConverter< QUANTITY_REACTION_RATE >::convert(
                       1., "m^3s^-1", "m^3s^-1") == 1.);

  assert_condition(UnitConverter< QUANTITY_REACTION_RATE >::convert(
                       1., "cm^3s^-1", "cm^3s^-1") == 1.);
  // a cubic centimetre per second is less than a cubic metre per second
  assert_condition(
      UnitConverter< QUANTITY_REACTION_RATE >::to_SI(1., "cm^3s^-1") < 1.);

  /// SURFACE AREA
  assert_condition(
      UnitConverter< QUANTITY_SURFACE_AREA >::convert(1., "m^2", "m^2") == 1.);

  assert_condition(UnitConverter< QUANTITY_SURFACE_AREA >::convert(
                       1., "cm^2", "cm^2") == 1.);
  // a square centimetre is less than a square metre
  assert_condition(UnitConverter< QUANTITY_SURFACE_AREA >::to_SI(1., "cm^2") <
                   1.);

  /// TEMPERATURE
  assert_condition(
      UnitConverter< QUANTITY_TEMPERATURE >::convert(1., "K", "K") == 1.);

  return 0;
}
