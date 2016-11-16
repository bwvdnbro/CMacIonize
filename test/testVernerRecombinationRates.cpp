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
 * @brief testVernerRecombinationRates.cpp
 *
 * @brief Unit test for the VernerRecombinationRates class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "UnitConverter.hpp"
#include "VernerRecombinationRates.hpp"
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

/**
 * @brief Unit test for the VernerRecombinationRates class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  VernerRecombinationRates recombination_rates;

  ifstream file("verner_rec_testdata.txt");
  string line;
  while (getline(file, line)) {
    if (line[0] != '#') {

      stringstream linestream(line);

      double T, alphaH, alphaHe, alphaCp1, alphaCp2, alphaN, alphaNp1, alphaNp2,
          alphaO, alphaOp1, alphaNe, alphaNep1, alphaSp1, alphaSp2, alphaSp3;

      linestream >> T >> alphaH >> alphaHe >> alphaCp1 >> alphaCp2 >> alphaN >>
          alphaNp1 >> alphaNp2 >> alphaO >> alphaOp1 >> alphaNe >> alphaNep1 >>
          alphaSp1 >> alphaSp2 >> alphaSp3;

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_H, T),
              "cm^3s^-1"),
          alphaH, 1.e-15);

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_He, T),
              "cm^3s^-1"),
          alphaHe, 1.e-15);

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Cp1, T),
              "cm^3s^-1"),
          alphaCp1, 1.e-15);

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Cp2, T),
              "cm^3s^-1"),
          alphaCp2, 1.e-15);

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_N, T),
              "cm^3s^-1"),
          alphaN, 1.e-15);
      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Np1, T),
              "cm^3s^-1"),
          alphaNp1, 1.e-15);
      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Np2, T),
              "cm^3s^-1"),
          alphaNp2, 1.e-15);

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_O, T),
              "cm^3s^-1"),
          alphaO, 1.e-15);
      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Op1, T),
              "cm^3s^-1"),
          alphaOp1, 1.e-15);

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Ne, T),
              "cm^3s^-1"),
          alphaNe, 1.e-15);
      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Nep1, T),
              "cm^3s^-1"),
          alphaNep1, 1.e-15);

      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Sp1, T),
              "cm^3s^-1"),
          alphaSp1, 1.e-15);
      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Sp2, T),
              "cm^3s^-1"),
          alphaSp2, 1.e-15);
      assert_values_equal_tol(
          UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
              recombination_rates.get_recombination_rate(ELEMENT_Sp3, T),
              "cm^3s^-1"),
          alphaSp3, 1.e-15);
    }
  }

  return 0;
}
