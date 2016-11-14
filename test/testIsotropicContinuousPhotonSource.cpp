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
 * @file testIsotropicContinuousPhotonSource.cpp
 *
 * @brief Unit test for the IsotropicContinuousPhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "IsotropicContinuousPhotonSource.hpp"
#include "RandomGenerator.hpp"
#include <fstream>

/**
 * @brief Unit test for the IsotropicContinuousPhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box box(CoordinateVector<>(), CoordinateVector<>(1.));
  RandomGenerator random_generator(42);
  IsotropicContinuousPhotonSource source(box, random_generator, 1.);

  const unsigned int numphoton = 10000000;

  unsigned int bins[6][101];
  for (unsigned int i = 0; i < 6; ++i) {
    for (unsigned int j = 0; j < 101; ++j) {
      bins[i][j] = 0;
    }
  }

  CoordinateVector<> average_position;
  CoordinateVector<> average_direction;
  for (unsigned int i = 0; i < numphoton; ++i) {
    std::pair< CoordinateVector<>, CoordinateVector<> > posdir =
        source.get_random_incoming_direction();
    average_position += posdir.first;
    average_direction += posdir.second;

    for (unsigned int ip = 0; ip < 3; ++ip) {
      unsigned int ibin = 100 * posdir.first[ip];
      ++bins[ip][ibin];
    }
    for (unsigned int id = 0; id < 3; ++id) {
      unsigned int ibin = 50 * (posdir.second[id] + 1.);
      ++bins[id + 3][ibin];
    }
  }
  average_position /= numphoton;
  average_direction /= numphoton;

  std::ofstream ofile("isotropic_bins.txt");
  unsigned int binsum[6] = {0};
  for (unsigned int i = 0; i < 101; ++i) {
    ofile << i;
    for (unsigned int j = 0; j < 6; ++j) {
      ofile << "\t" << bins[j][i];
      binsum[j] += bins[j][i];
    }
    ofile << "\n";
  }

  return 0;
}
