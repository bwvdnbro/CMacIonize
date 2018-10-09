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
 * @file testCaproniStellarFeedback.cpp
 *
 * @brief Unit test for the CaproniStellarFeedback class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CaproniStellarFeedback.hpp"
#include "fstream"

/**
 * @brief Unit test for the CaproniStellarFeedback class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double Myr_in_s = 1.e6 * 365.25 * 24. * 3600.;

  /// test SN rate function
  {
    std::ofstream ofile("test_caproni_rate.txt");
    ofile << "# t (s)\tsupernova rate (s^-1)\n";
    const double dt = 2. * Myr_in_s;
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      const double t = i * dt;
      const double rate = CaproniStellarFeedback::get_SN_rate(t);
      ofile << t << "\t" << rate << "\n";
    }
  }

  CaproniStellarFeedback stellar_feedback;

  return 0;
}
