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
 * @file testCMILibrary.cpp
 *
 * @brief Unit test for the CMILibrary.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "CMILibrary.hpp"

#include <fstream>
#include <vector>

/**
 * @brief Unit test for the CMILibrary.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // set up a test
  const double pc = 3.086e16;
  std::vector< double > x(1000, 0.);
  std::vector< double > y(1000, 0.);
  std::vector< double > z(1000, 0.);
  std::vector< double > h(1000, 0.);
  std::vector< double > m(1000, 0.);
  const double box_anchor[3] = {-5. * pc, -5. * pc, -5. * pc};
  const double box_sides[3] = {10. * pc, 10. * pc, 10. * pc};
  for (uint_fast8_t ix = 0; ix < 10; ++ix) {
    for (uint_fast8_t iy = 0; iy < 10; ++iy) {
      for (uint_fast8_t iz = 0; iz < 10; ++iz) {
        const uint_fast32_t i = ix * 100 + iy * 10 + iz;
        x[i] = box_anchor[0] + 0.1 * (ix + 0.5) * box_sides[0];
        y[i] = box_anchor[1] + 0.1 * (iy + 0.5) * box_sides[1];
        z[i] = box_anchor[2] + 0.1 * (iz + 0.5) * box_sides[2];
        h[i] = 0.2 * box_sides[0];
        // 100. cm^-3 * (10.pc)^3 * 1.67*10^{-27} kg / 1000
        m[i] = 4.9e30;
      }
    }
  }

  // initialize the library
  cmi_init_periodic_dp("test_CMI_library.param", 1, 1., 1., box_anchor,
                       box_sides, "Petkova", false);

  // run the simulation
  std::vector< double > nH(1000, 0.);
  cmi_compute_neutral_fraction_dp(x.data(), y.data(), z.data(), h.data(),
                                  m.data(), nH.data(), 1000);

  // write an output file for visual checking
  std::ofstream ofile("test_CMI_library.txt");
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    ofile << x[i] << "\t" << y[i] << "\t" << z[i] << "\t" << nH[i] << "\n";
  }
  ofile.close();

  // clean up the library
  cmi_destroy();

  return 0;
}
