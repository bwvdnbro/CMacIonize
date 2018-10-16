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
#include "CartesianDensityGrid.hpp"
#include "HomogeneousDensityFunction.hpp"
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
  const double kpc_in_m = 3.086e19;

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

  CaproniStellarFeedback stellar_feedback(42, 1., true);

  HomogeneousDensityFunction testfunction(1., 2000.);
  testfunction.initialize();
  CartesianDensityGrid grid(Box<>(-1.5 * kpc_in_m, 3. * kpc_in_m), 32, false,
                            true);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, testfunction);

  stellar_feedback.add_stellar_feedback(grid, 200. * Myr_in_s);
  double N_SN = 0.;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    N_SN += 1.e-44 * it.get_hydro_variables().get_energy_term();
  }
  cmac_status("N_SN: %.0f", N_SN);

  return 0;
}
