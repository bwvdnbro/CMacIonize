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
 * We test our own implementation of the Kingdon & Ferland (1996) fitting
 * functions against output from Kingdon's ct1.f script that was part of Kenny's
 * code.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  ChargeTransferRates rates;

  // table that links atom and stage indices to IonNames
  IonName ion[17][6];
  ion[6][4] = ION_C_p2;
  ion[7][1] = ION_N_n;
  ion[7][2] = ION_N_n;
  ion[7][3] = ION_N_p1;
  ion[7][4] = ION_N_p2;
  ion[8][1] = ION_O_n;
  ion[8][2] = ION_O_n;
  ion[8][3] = ION_O_p1;
  ion[10][3] = ION_Ne_p1;
  ion[16][3] = ION_S_p1;
  ion[16][4] = ION_S_p2;
  ion[16][5] = ION_S_p3;
  std::ifstream file("KingdonFerland_testdata.txt");
  std::string line;
  while (getline(file, line)) {
    if (line[0] != '#') {

      std::istringstream linestream(line);

      unsigned int stage, atom;
      double temperature, recombination_rate, ionization_rate;

      linestream >> stage >> atom >> temperature >> recombination_rate >>
          ionization_rate;

      // since we test both ionization and recombination, there are some stages
      // in the test file that do not have a recombination rate
      // we don't test these
      if (stage > 1) {
        assert_values_equal_rel(
            recombination_rate,
            UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
                rates.get_charge_transfer_recombination_rate_H(
                    ion[atom][stage], temperature * 1.e-4),
                "cm^3s^-1"),
            1.e-6);
      }

      // there are only two ions that have ionization rates: N_n and O_n
      // we only test these two
      if ((atom == 7 && stage == 1) || (atom == 8 && stage == 1)) {

        assert_values_equal_rel(
            ionization_rate, UnitConverter::to_unit< QUANTITY_REACTION_RATE >(
                                 rates.get_charge_transfer_ionization_rate_H(
                                     ion[atom][stage], temperature * 1.e-4),
                                 "cm^3s^-1"),
            1.e-6);
      }
    }
  }

  return 0;
}
