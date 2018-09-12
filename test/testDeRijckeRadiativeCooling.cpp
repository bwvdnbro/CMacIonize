/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testDeRijckeRadiativeCooling.cpp
 *
 * @brief Unit test for the DeRijckeRadiativeCooling class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DeRijckeRadiativeCooling.hpp"
#include <fstream>

/**
 * @brief Unit test for the DeRijckeRadiativeCooling class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  DeRijckeRadiativeCooling cooling;

  std::ofstream ofile("test_derijckeradiativecooling.txt");
  ofile << "# T (K)\tLambda (J m^-3 s^-1)\n";
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    const double T = 100. + i * 1000.;
    ofile << T << "\t" << cooling.get_cooling_rate(T) << "\n";
  }

  return 0;
}
