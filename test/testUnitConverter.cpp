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
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  assert_condition(
      UnitConverter< QUANTITY_FREQUENCY >::convert(1., "s^-1", "s^-1") == 1.);
  assert_condition(
      UnitConverter< QUANTITY_FREQUENCY >::convert(1., "Hz", "Hz") == 1.);
  assert_condition(
      UnitConverter< QUANTITY_FREQUENCY >::convert(1., "eV", "eV") == 1.);

  assert_condition(UnitConverter< QUANTITY_LENGTH >::convert(1., "m", "m") ==
                   1.);
  assert_condition(UnitConverter< QUANTITY_LENGTH >::convert(1., "cm", "cm") ==
                   1.);
  assert_condition(UnitConverter< QUANTITY_LENGTH >::convert(1., "pc", "pc") ==
                   1.);
  assert_condition(
      UnitConverter< QUANTITY_LENGTH >::convert(1., "kpc", "kpc") == 1.);

  assert_condition(UnitConverter< QUANTITY_MASS >::convert(1., "kg", "kg") ==
                   1.);
  assert_condition(UnitConverter< QUANTITY_MASS >::convert(1., "g", "g") == 1.);

  assert_condition(
      UnitConverter< QUANTITY_NUMBER_DENSITY >::convert(1., "m^-3", "m^-3"));
  assert_condition(
      UnitConverter< QUANTITY_NUMBER_DENSITY >::convert(1., "cm^-3", "cm^-3"));

  assert_condition(UnitConverter< QUANTITY_REACTION_RATE >::convert(
                       1., "m^3s^-1", "m^3s^-1") == 1.);
  assert_condition(UnitConverter< QUANTITY_REACTION_RATE >::convert(
                       1., "cm^3s^-1", "cm^3s^-1") == 1.);

  assert_condition(
      UnitConverter< QUANTITY_TEMPERATURE >::convert(1., "K", "K") == 1.);

  return 0;
}
