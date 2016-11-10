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
 * @file testIonizationStateCalculator.cpp
 *
 * @brief Unit test for the IonizationStateCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "IonizationStateCalculator.hpp"
#include "UnitConverter.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the IonizationStateCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // test find_H0
  std::ifstream file("h0_testdata.txt");
  std::string line;
  while (getline(file, line)) {
    std::stringstream linestream(line);
    double jH, jHe, h0f, he0f, h0, he0;
    linestream >> jH >> jHe >> h0f >> he0f;
    IonizationStateCalculator::find_H0(
        UnitConverter< QUANTITY_REACTION_RATE >::to_SI(3.12e-13, "cm^3s^-1"),
        UnitConverter< QUANTITY_REACTION_RATE >::to_SI(3.51e-13, "cm^3s^-1"),
        UnitConverter< QUANTITY_FREQUENCY >::to_SI(jH, "s^-1"),
        UnitConverter< QUANTITY_FREQUENCY >::to_SI(jHe, "s^-1"),
        UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(100., "cm^-3"), 0.1,
        1.5e4, h0, he0);
    assert_values_equal(h0, h0f);
    assert_values_equal(he0, he0f);

    double h0s;
    IonizationStateCalculator::find_H0(
        UnitConverter< QUANTITY_REACTION_RATE >::to_SI(3.12e-13, "cm^3s^-1"),
        0., UnitConverter< QUANTITY_FREQUENCY >::to_SI(jH, "s^-1"), 0.,
        UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(100., "cm^-3"), 0.,
        1.5e4, h0, he0);
    IonizationStateCalculator::find_H0_simple(
        UnitConverter< QUANTITY_REACTION_RATE >::to_SI(3.12e-13, "cm^3s^-1"),
        UnitConverter< QUANTITY_FREQUENCY >::to_SI(jH, "s^-1"),
        UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(100., "cm^-3"), 1.5e4,
        h0s);
    assert_values_equal_tol(h0, h0s, 1.e-7);
  }

  return 0;
}
