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
 * @file testLambertW.cpp
 *
 * @brief Unit test for the Lambert W function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "LambertW.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the Lambert W function.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  std::ifstream ifile("test_lambertw_input.txt");
  std::string line;
  // skip two comment lines
  std::getline(ifile, line);
  std::getline(ifile, line);
  double xref[1000], w0ref[1000], wm1ref[1000];
  unsigned int i = 0;
  while (std::getline(ifile, line)) {
    std::stringstream lstream(line);
    lstream >> xref[i] >> w0ref[i] >> wm1ref[i];
    ++i;
  }
  assert_condition(i == 1000);

  std::ofstream ofile("test_lambertw_output.txt");
  for (unsigned int i = 0; i < 1000; ++i) {
    const double x = (-1. + i * 0.001) / M_E;
    const double w0 = LambertW::lambert_w(x, 0);
    const double wm1 = LambertW::lambert_w(x, -1);
    assert_values_equal_rel(x, xref[i], 1.e-11);
    assert_values_equal_rel(w0, w0ref[i], 1.e-7);
    assert_values_equal_rel(wm1, wm1ref[i], 1.e-8);
    ofile << x << "\t" << w0 << "\t" << wm1 << "\t" << w0ref[i] << "\t"
          << wm1ref[i] << "\n";
  }

  return 0;
}
