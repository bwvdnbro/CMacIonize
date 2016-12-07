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

/**
 * @brief Unit test for the VernerRecombinationRates class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  VernerRecombinationRates recombination_rates;

  std::ifstream file("verner_rec_testdata.txt");
  std::string line;
  while (getline(file, line)) {
    if (line[0] != '#') {

      std::istringstream linestream(line);

      double T, alphaH, alphaHe, alphaCp1, alphaCp2, alphaN, alphaNp1, alphaNp2,
          alphaO, alphaOp1, alphaNe, alphaNep1, alphaSp1, alphaSp2, alphaSp3;

      linestream >> T >> alphaH >> alphaHe >> alphaCp1 >> alphaCp2 >> alphaN >>
          alphaNp1 >> alphaNp2 >> alphaO >> alphaOp1 >> alphaNe >> alphaNep1 >>
          alphaSp1 >> alphaSp2 >> alphaSp3;

      double tolerance = 1.e-15;

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_H_n, T),
              "cm^3s^-1"),
          alphaH, tolerance);

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_He_n, T),
              "cm^3s^-1"),
          alphaHe, tolerance);

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_C_p1, T),
              "cm^3s^-1"),
          alphaCp1, tolerance);

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_C_p2, T),
              "cm^3s^-1"),
          alphaCp2, tolerance);

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_N_n, T),
              "cm^3s^-1"),
          alphaN, tolerance);
      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_N_p1, T),
              "cm^3s^-1"),
          alphaNp1, tolerance);
      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_N_p2, T),
              "cm^3s^-1"),
          alphaNp2, tolerance);

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_O_n, T),
              "cm^3s^-1"),
          alphaO, tolerance);
      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_O_p1, T),
              "cm^3s^-1"),
          alphaOp1, tolerance);

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_Ne_n, T),
              "cm^3s^-1"),
          alphaNe, tolerance);
      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_Ne_p1, T),
              "cm^3s^-1"),
          alphaNep1, tolerance);

      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_S_p1, T),
              "cm^3s^-1"),
          alphaSp1, tolerance);
      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_S_p2, T),
              "cm^3s^-1"),
          alphaSp2, tolerance);
      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
              recombination_rates.get_recombination_rate(ION_S_p3, T),
              "cm^3s^-1"),
          alphaSp3, tolerance);
    }
  }

  return 0;
}
