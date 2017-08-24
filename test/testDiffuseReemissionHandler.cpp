/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testDiffuseReemissionHandler.cpp
 *
 * @brief Unit test for the DiffuseReemissionHandler class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "DiffuseReemissionHandler.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the DiffuseReemissionHandler class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  IonizationVariables ionization_variables;

  const double tolerance = 1.e-15;
  std::ifstream ifile("probset_testdata.txt");
  std::string line;
  while (std::getline(ifile, line)) {
    std::stringstream sstream(line);

    double T, pH, pHe[4];
    sstream >> T >> pH >> pHe[0] >> pHe[1] >> pHe[2] >> pHe[3];

    ionization_variables.set_temperature(T);
    DiffuseReemissionHandler::set_reemission_probabilities(
        ionization_variables);

    assert_values_equal_rel(ionization_variables.get_reemission_probability(
                                REEMISSIONPROBABILITY_HYDROGEN),
                            pH, tolerance);
    assert_values_equal_rel(ionization_variables.get_reemission_probability(
                                REEMISSIONPROBABILITY_HELIUM_LYC),
                            pHe[0], tolerance);
    assert_values_equal_rel(ionization_variables.get_reemission_probability(
                                REEMISSIONPROBABILITY_HELIUM_NPEEV),
                            pHe[1], tolerance);
    assert_values_equal_rel(ionization_variables.get_reemission_probability(
                                REEMISSIONPROBABILITY_HELIUM_TPC),
                            pHe[2], tolerance);
    assert_values_equal_rel(ionization_variables.get_reemission_probability(
                                REEMISSIONPROBABILITY_HELIUM_LYA),
                            pHe[3], tolerance);
  }

  return 0;
}
