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
 * @file testPhotonSourceSpectrum.cpp
 *
 * @brief Unit test for the PhotonSourceSpectrum interface and its
 * implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Error.hpp"
#include "HydrogenLymanContinuumSpectrum.hpp"
#include "PlanckPhotonSourceSpectrum.hpp"
#include "VernerCrossSections.hpp"
#include <cmath>
#include <fstream>
using namespace std;

/**
 * @brief Get the Planck black body luminosity for a given frequency.
 *
 * @param frequency Frequency value.
 * @return Planck luminosity.
 */
double planck_luminosity(double frequency) {
  double min_frequency = 3.289e15;
  double planck_constant = 6.626e-27;
  double boltzmann_constant = 1.38e-16;
  double temperature_star = 40000.;
  return frequency * frequency * frequency /
         (exp(planck_constant * frequency * min_frequency /
              (boltzmann_constant * temperature_star)) -
          1.);
}

/**
 * @brief Get the hydrogen Lyman continuum luminosity at the given temperature
 * and for the given frequency.
 *
 * @param cross_sections Photoionization cross sections.
 * @param T Temperature.
 * @param frequency Frequency.
 * @return Hydrogen Lyman continuum luminosity.
 */
double HLyc_luminosity(CrossSections &cross_sections, double T,
                       double frequency) {
  double xsecH = cross_sections.get_cross_section(ELEMENT_H, frequency * 13.6);
  return frequency * frequency * frequency * xsecH *
         exp(-157919.667 * (frequency - 1.) / T);
}

/**
 * @brief Unit test for the PhotonSourceSpectrum interface and its
 * implementations.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // PlanckPhotonSourceSpectrum
  {
    PlanckPhotonSourceSpectrum spectrum;

    unsigned int counts[100];
    for (unsigned int i = 0; i < 100; ++i) {
      counts[i] = 0;
    }
    unsigned int numsample = 1000000;
    for (unsigned int i = 0; i < numsample; ++i) {
      double rand_freq = spectrum.get_random_frequency();
      unsigned int index = (rand_freq - 1.) * 100. / 3.;
      ++counts[index];
    }

    double enorm = planck_luminosity(1.);
    if (counts[0]) {
      enorm /= counts[0];
    }
    for (unsigned int i = 0; i < 100; ++i) {
      double nu = 1. + i * 0.03;
      assert_values_equal_tol(planck_luminosity(nu), counts[i] * enorm, 1.e-2);
    }
  }

  // HydrogenLymanContinuumSpectrum
  {
    VernerCrossSections cross_sections;
    HydrogenLymanContinuumSpectrum spectrum(cross_sections);
    for (unsigned int iT = 0; iT < 10; ++iT) {
      double T = 1500. + (iT + 0.5) * 13500. / 10;
      spectrum.set_temperature(T);

      unsigned int counts[100];
      for (unsigned int i = 0; i < 100; ++i) {
        counts[i] = 0;
      }
      unsigned int numsample = 1000000;
      for (unsigned int i = 0; i < numsample; ++i) {
        double rand_freq = spectrum.get_random_frequency();
        unsigned int index = (rand_freq - 1.) * 100. / 3.;
        ++counts[index];
      }

      double enorm = HLyc_luminosity(cross_sections, T, 1.03);
      if (counts[1]) {
        enorm /= counts[1];
      }
      for (unsigned int i = 0; i < 100; ++i) {
        double nu = 1. + i * 0.03;
        assert_values_equal_tol(HLyc_luminosity(cross_sections, T, nu),
                                counts[i] * enorm, 1.e-1);
      }
    }
  }

  return 0;
}
