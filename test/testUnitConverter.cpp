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
#include <string>
#include <vector>

/**
 * @brief Perform unit tests for the given template quantity.
 *
 * For every unit name in the given list, we perform 3 checks:
 *  - we check if the conversion of the unit name into itself yields 1.
 *  - we do the same by explicitly converting from the unit to SI units and then
 *    back to the unit
 *  - we check if the SI value of the unit matches the expected logic: it is
 *    either larger, smaller, or equal to the default SI unit for the quantity
 *
 * @param unitnames Names of the units to check.
 * @param unitlogic Flags indicating if the SI value of the corresponding unit
 * is equal, smaller or larger than 1.
 */
template < Quantity _quantity_ >
void check_quantity(std::vector< std::string > unitnames,
                    std::vector< int > unitlogic) {
  for (unsigned int i = 0; i < unitnames.size(); ++i) {
    assert_condition(UnitConverter::convert(1., unitnames[i], unitnames[i]) ==
                     1.);
    assert_condition(UnitConverter::to_unit< _quantity_ >(
                         UnitConverter::to_SI< _quantity_ >(1., unitnames[i]),
                         unitnames[i]) == 1.);
    double SI_value = UnitConverter::to_SI< _quantity_ >(1., unitnames[i]);
    if (unitlogic[i] > 0) {
      assert_condition(SI_value > 1.);
    } else if (unitlogic[i] < 0) {
      assert_condition(SI_value < 1.);
    } else {
      assert_condition(SI_value == 1.);
    }
  }
}

/**
 * @brief Unit test for the UnitConverter class.
 *
 * The UnitConverter class supports three methods:
 *   - a method that converts from a strange unit to SI units
 *   - a method that converts from SI units to a strange unit
 *   - a method that converts from one strange unit into another strange unit
 * To check whether the methods are correctly implemented, we only need to check
 * one quantity. To check that the separate units are correctly implemented, we
 * just check all of them.
 *
 * To check if the conversion from one unit into SI units is correct, we apply
 * larger/smaller/equal logic: for each unit we check if the SI value of that
 * unit is larger, smaller, or equal to 1. We have to specify the correct logic
 * for each unit separately.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  {
    Unit unit = UnitConverter::get_unit(" m");
    Unit unit2(1., 1, 0, 0, 0, 0);
    assert_condition(unit == unit2);
  }

  {
    Unit unit = UnitConverter::get_unit("K kg^3 s^-1m ");
    Unit unit2(1., 1, -1, 3, 1, 0);
    assert_condition(unit == unit2);
  }
  {
    assert_values_equal_rel(UnitConverter::convert(1., "eV", "Hz"), 2.417989e14,
                            1.e-7);
    assert_values_equal_rel(UnitConverter::convert(1., "Hz", "eV"),
                            4.135668e-15, 1.e-7);

    // the conversion from angstrom to Hz is non-linear, so we test with a value
    // different than 1.
    assert_values_equal_rel(UnitConverter::convert(2., "angstrom", "Hz"),
                            1.499e18, 1.e-4);
    assert_values_equal_rel(UnitConverter::convert(2., "Hz", "angstrom"),
                            1.499e18, 1.e-4);
  }

  /// DENSITY
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("kg m^-3");
    unitlogic.push_back(0);
    unitnames.push_back("g cm^-3");
    // 1 g in a cubic centimetre is more than 1 kg in a cubic metre
    unitlogic.push_back(1);

    check_quantity< QUANTITY_DENSITY >(unitnames, unitlogic);
  }

  /// ENERGY
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("J");
    unitlogic.push_back(0);
    unitnames.push_back("eV");
    unitlogic.push_back(-1);
    unitnames.push_back("erg");
    unitlogic.push_back(-1);

    check_quantity< QUANTITY_ENERGY >(unitnames, unitlogic);
  }

  /// ENERGY CHANGE RATE
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("kg m^-1s^-3");
    unitlogic.push_back(0);
    unitnames.push_back("J m^-3s^-1");
    // J m^-3s^-1 is just another name for kg m^-1s^-3
    unitlogic.push_back(0);
    unitnames.push_back("erg cm^-3s^-1");
    // according to WolframAlpha, an erg per cubic centimetres per second is
    // less than a Joule per cubic metres per second
    unitlogic.push_back(-1);

    check_quantity< QUANTITY_ENERGY_CHANGE_RATE >(unitnames, unitlogic);
  }

  /// ENERGY_RATE
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("kg m^2s^-3");
    unitlogic.push_back(0);
    unitnames.push_back("J s^-1");
    // J s^-1 is just another name for kg m^2s^-3
    unitlogic.push_back(0);
    unitnames.push_back("erg s^-1");
    // an erg per second is a lot less than a Joule per second
    unitlogic.push_back(-1);

    check_quantity< QUANTITY_ENERGY_RATE >(unitnames, unitlogic);
  }

  /// FLUX
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("m^-2s^-1");
    unitlogic.push_back(0);
    unitnames.push_back("cm^-2s^-1");
    // 1 photon per square centimetres per second is more than 1 photon per
    // square metres per second
    unitlogic.push_back(1);

    check_quantity< QUANTITY_FLUX >(unitnames, unitlogic);
  }

  /// FREQUENCY
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("s^-1");
    unitlogic.push_back(0);
    unitnames.push_back("Hz");
    // Hz is just another name for s^-1
    unitlogic.push_back(0);
    unitnames.push_back("eV");
    // 1 eV corresponds to a very high frequency in Hz
    unitlogic.push_back(1);

    check_quantity< QUANTITY_FREQUENCY >(unitnames, unitlogic);
  }

  /// LENGTH
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("m");
    unitlogic.push_back(0);
    unitnames.push_back("cm");
    // cm is smaller than m
    unitlogic.push_back(-1);
    unitnames.push_back("pc");
    // a parsec is a *bit* larger than a metre
    unitlogic.push_back(1);
    unitnames.push_back("kpc");
    // a kpc is even larger
    unitlogic.push_back(1);
    unitnames.push_back("angstrom");
    // an angstrom is orders of magnitude smaller than a metre
    unitlogic.push_back(-1);

    check_quantity< QUANTITY_LENGTH >(unitnames, unitlogic);
  }

  /// MASS
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("kg");
    unitlogic.push_back(0);
    unitnames.push_back("g");
    // a gramme is less than a kilogram
    unitlogic.push_back(-1);
    unitnames.push_back("Msol");
    // a solar mass weighs much more than a kilogram
    unitlogic.push_back(1);

    check_quantity< QUANTITY_MASS >(unitnames, unitlogic);
  }

  /// NUMBER DENSITY
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("m^-3");
    unitlogic.push_back(0);
    unitnames.push_back("cm^-3");
    // 1 particle in a cubic centimetre is much more than 1 particle in a cubic
    // metre
    unitlogic.push_back(1);

    check_quantity< QUANTITY_NUMBER_DENSITY >(unitnames, unitlogic);
  }

  /// OPACITY
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("m^-1");
    unitlogic.push_back(0);

    check_quantity< QUANTITY_OPACITY >(unitnames, unitlogic);
  }

  /// REACTION RATE
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("m^3s^-1");
    unitlogic.push_back(0);
    unitnames.push_back("cm^3s^-1");
    // a cubic centimetre per second is less than a cubic metre per second
    unitlogic.push_back(-1);

    check_quantity< QUANTITY_REACTION_RATE >(unitnames, unitlogic);
  }

  /// SURFACE AREA
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("m^2");
    unitlogic.push_back(0);
    unitnames.push_back("cm^2");
    // a square centimetre is less than a square metre
    unitlogic.push_back(-1);

    check_quantity< QUANTITY_SURFACE_AREA >(unitnames, unitlogic);
  }

  /// TEMPERATURE
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("K");
    unitlogic.push_back(0);

    check_quantity< QUANTITY_TEMPERATURE >(unitnames, unitlogic);
  }

  /// TIME
  {
    std::vector< std::string > unitnames;
    std::vector< int > unitlogic;
    unitnames.push_back("s");
    unitlogic.push_back(0);
    unitnames.push_back("Gyr");
    // a Gyr takes longer than a second
    unitlogic.push_back(1);

    check_quantity< QUANTITY_TIME >(unitnames, unitlogic);
  }

  return 0;
}
