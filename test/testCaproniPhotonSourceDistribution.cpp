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
 * @file testCaproniPhotonSourceDistribution.cpp
 *
 * @brief Unit test for the CaproniPhotonSourceDistribution class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CaproniPhotonSourceDistribution.hpp"
#include "CartesianDensityGrid.hpp"
#include "Error.hpp"
#include "HomogeneousDensityFunction.hpp"
#include <fstream>

/**
 * @brief Unit test for the CaproniPhotonSourceDistribution class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double Myr_in_s = 1.e6 * 365.25 * 24. * 3600.;

  /// test number of stars function
  {
    std::ofstream ofile("test_caproni_numstars.txt");
    ofile << "# t (Myr)\tnumber of stars\n";
    const uint_fast32_t num_steps = 1000;
    const double dt = 2. * Myr_in_s;
    double average_number_of_sources = 0.;
    for (uint_fast32_t i = 0; i < num_steps; ++i) {
      const double t = i * dt;
      const uint_fast32_t nstar =
          CaproniPhotonSourceDistribution::get_number_of_stars(t);
      ofile << t << "\t" << nstar << "\n";
      average_number_of_sources += nstar;
    }
    average_number_of_sources /= num_steps;
    cmac_status("Average number of stars: %g", average_number_of_sources);
  }

  /// test random star function
  {
    RandomGenerator random_generator;
    std::ofstream ofile("test_caproni_star.txt");
    ofile << "# mass (Msol)\tluminosity (s^-1)\tlifetime (s)\n";
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      double luminosity, lifetime;
      const double mass = CaproniPhotonSourceDistribution::get_random_star(
          luminosity, lifetime, random_generator);
      ofile << mass << "\t" << luminosity << "\t" << lifetime << "\n";
    }
  }

  /// test random galactic radius function
  {
    RandomGenerator random_generator;
    std::ofstream ofile("test_caproni_radius.txt");
    ofile << "#t (s)\travg (m)\trstd (m)\n";
    const uint_fast32_t num_steps = 1000;
    const double dt = 2. * Myr_in_s;
    for (uint_fast32_t i = 0; i < num_steps; ++i) {
      const double t = i * dt;
      double ravg = 0.;
      double ravg2 = 0.;
      for (uint_fast32_t j = 0; j < 100; ++j) {
        const double r = CaproniPhotonSourceDistribution::get_galactic_radius(
            t, random_generator);
        ravg += r;
        ravg2 += r * r;
      }
      ravg *= 0.01;
      ravg2 *= 0.01;
      const double rstd = std::sqrt(ravg2 - ravg * ravg);
      ofile << t << "\t" << ravg << "\t" << rstd << "\n";
    }
  }

  CaproniPhotonSourceDistribution distribution(0.01, 1., 42, 10. * Myr_in_s, 0.,
                                               1., true);

  std::ofstream ofile("test_caproniphotonsourcedistribution.txt");
  ofile << "# t (Myr)\tnumber of sources\n";
  const uint_fast32_t num_steps = 100;
  const double dt = 20. * Myr_in_s;
  double average_number_of_sources = 0.;
  for (uint_fast32_t i = 0; i < num_steps; ++i) {
    const double t = i * dt;
    ofile << t << "\t" << distribution.get_number_of_sources() << "\n";
    average_number_of_sources += distribution.get_number_of_sources();
    distribution.update(t);
  }
  average_number_of_sources /= num_steps;
  cmac_status("Average number of sources: %g", average_number_of_sources);

  const double kpc_in_m = 3.086e19;

  /// test SN rate function
  {
    std::ofstream ofile("test_caproni_rate.txt");
    ofile << "# t (s)\tsupernova rate (s^-1)\n";
    const double dt = 2. * Myr_in_s;
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      const double t = i * dt;
      const double rate = CaproniPhotonSourceDistribution::get_SN_rate(t);
      ofile << t << "\t" << rate << "\n";
    }
  }

  HomogeneousDensityFunction testfunction(1., 2000.);
  testfunction.initialize();
  CartesianDensityGrid grid(Box<>(-1.5 * kpc_in_m, 3. * kpc_in_m), 32, false,
                            true);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, testfunction);

  distribution.add_stellar_feedback(grid, 200. * Myr_in_s);
  double N_SN = 0.;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    N_SN += 1.e-44 * it.get_hydro_variables().get_energy_term();
  }
  cmac_status("N_SN: %.0f", N_SN);

  return 0;
}
