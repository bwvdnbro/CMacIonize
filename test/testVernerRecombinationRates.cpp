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
      double T, alphaH, alphaHe;
      stringstream linestream(line);
      linestream >> T >> alphaH >> alphaHe;
      assert_values_equal_tol(
          recombination_rates.get_recombination_rate(ELEMENT_H, T), alphaH,
          1.e-15);
      assert_values_equal_tol(
          recombination_rates.get_recombination_rate(ELEMENT_He, T), alphaHe,
          1.e-15);
    }
  }

  return 0;
}
