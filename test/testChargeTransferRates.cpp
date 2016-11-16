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
 * @file testChargeTransferRates.cpp
 *
 * @brief Unit test for the ChargeTransferRates class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "ChargeTransferRates.hpp"
#include "UnitConverter.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the ChargeTransferRates class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  ChargeTransferRates rates;

  std::ifstream file("KingdonFerland_testdata.txt");
  std::string line;
  while (getline(file, line)) {
    if (line[0] != '#') {

      std::istringstream linestream(line);

      unsigned int stage, atom;
      double temperature, recombination_rate, ionization_rate;

      linestream >> stage >> atom >> temperature >> recombination_rate >>
          ionization_rate;

      assert_values_equal_tol(recombination_rate,
                              UnitConverter< QUANTITY_REACTION_RATE >::to_unit(
                                  rates.get_charge_transfer_recombination_rate(
                                      stage, atom, temperature),
                                  "cm^3s^-1"),
                              1.e-14);
      assert_values_equal_tol(
          ionization_rate,
          rates.get_charge_transfer_ionization_rate(stage, atom, temperature),
          1.e-14);
    }
  }

  return 0;
}
