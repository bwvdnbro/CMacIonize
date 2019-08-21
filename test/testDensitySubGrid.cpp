/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testDensitySubGrid.cpp
 *
 * @brief Unit test for the DensitySubGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "DensitySubGrid.hpp"
#include "RandomGenerator.hpp"

#include <fstream>

/**
 * @brief Unit test for the DensitySubGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double box[6] = {-1.543e17, -1.543e17, -1.543e17,
                         3.086e17,  3.086e17,  3.086e17};
  const CoordinateVector< int_fast32_t > ncell(16, 16, 16);
  DensitySubGrid grid(box, ncell);
  RandomGenerator random_generator(42);

  for (uint_fast32_t iloop = 0; iloop < 10; ++iloop) {
    for (uint_fast32_t i = 0; i < 1e5; ++i) {
      PhotonPacket photon;

      const double cost =
          2. * random_generator.get_uniform_random_double() - 1.;
      const double phi =
          2. * M_PI * random_generator.get_uniform_random_double();
      const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
      const double cosp = std::cos(phi);
      const double sinp = std::sin(phi);
      const double dx = sint * cosp;
      const double dy = sint * sinp;
      const double dz = cost;

      const double tau =
          -std::log(random_generator.get_uniform_random_double());

      photon.set_position(0.5, 0.5, 0.5);
      photon.set_direction(dx, dy, dz);
      photon.set_photoionization_cross_section(6.3e-22);
      photon.set_weight(1.);
      photon.set_target_optical_depth(tau);

      grid.interact(photon, TRAVELDIRECTION_INSIDE);
    }

    grid.compute_neutral_fraction(1.e49, 1e5);
  }

  std::ofstream ofile("testDensitySubGrid_output.txt");
  ofile << "# x (m)\ty (m)\tz (m)\txH\tnH (m^-3)\n";
  grid.print_intensities(ofile);

  return 0;
}
