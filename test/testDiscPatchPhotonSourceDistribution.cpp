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
 * @file testDiscPatchPhotonSourceDistribution.cpp
 *
 * @brief Unit test for the DiscPatchPhotonSourceDistribution class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DiscPatchPhotonSourceDistribution.hpp"
#include "Error.hpp"
#include <fstream>

/**
 * @brief Unit test for the DiscPatchPhotonSourceDistribution class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double Myr_in_s = 1.e6 * 365.25 * 24. * 3600.;
  DiscPatchPhotonSourceDistribution distribution(20. * Myr_in_s, 4.26e49, 6, 0.,
                                                 1., 0., 1., 0., 1., 42,
                                                 0.1 * Myr_in_s, true);

  std::ofstream ofile("test_disc_patch_photon_source_distribution.txt");
  ofile << "# t (Myr)\tnumber of sources\n";
  const uint_fast32_t num_steps = 1000;
  const double dt = 0.2 * Myr_in_s;
  double average_number_of_sources = 0.;
  for (uint_fast32_t i = 0; i < num_steps; ++i) {
    const double t = i * dt;
    ofile << t << "\t" << distribution.get_number_of_sources() << "\n";
    average_number_of_sources += distribution.get_number_of_sources();
    distribution.update(t);
  }
  average_number_of_sources /= num_steps;
  cmac_status("Average number of sources: %g", average_number_of_sources);

  return 0;
}
