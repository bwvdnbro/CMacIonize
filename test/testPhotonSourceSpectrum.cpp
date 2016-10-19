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
#include "HeliumLymanContinuumSpectrum.hpp"
#include "HeliumTwoPhotonContinuumSpectrum.hpp"
#include "HydrogenLymanContinuumSpectrum.hpp"
#include "PlanckPhotonSourceSpectrum.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"
#include "VernerCrossSections.hpp"
#include <cmath>
#include <fstream>
#include <vector>
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
 * @brief Get the helium Lyman continuum luminosity at the given temperature
 * and for the given frequency.
 *
 * @param cross_sections Photoionization cross sections.
 * @param T Temperature.
 * @param frequency Frequency.
 * @return Helium Lyman continuum luminosity.
 */
double HeLyc_luminosity(CrossSections &cross_sections, double T,
                        double frequency) {
  if (frequency >= 1.81) {
    double xsecHe =
        cross_sections.get_cross_section(ELEMENT_He, frequency * 13.6);
    return frequency * frequency * frequency * xsecHe *
           exp(-157919.667 * (frequency - 1.81) / T);
  } else {
    return 0.;
  }
}

/**
 * @brief Get the helium 2-photon continuum luminosity.
 *
 * @param yHe2q y values of the spectrum.
 * @param AHe2q A values of the spectrum.
 * @param frequency Frequency value.
 * @return Helium 2-photon continuum luminosity.
 */
double He2pc_luminosity(vector< double > &yHe2q, vector< double > &AHe2q,
                        double frequency) {
  double y = frequency * 3.289e15 / 4.98e15;
  if (y < 1.) {
    unsigned int i = Utilities::locate(y, &yHe2q[0], 41);
    double f = (y - yHe2q[i]) / (yHe2q[i + 1] - yHe2q[i]);
    return AHe2q[i] + f * (AHe2q[i + 1] - AHe2q[i]);
  } else {
    return 0.;
  }
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
      double rand_freq = UnitConverter< QUANTITY_FREQUENCY >::to_unit(
                             spectrum.get_random_frequency(), "eV") /
                         13.6;
      unsigned int index = (rand_freq - 1.) * 100. / 3.;
      ++counts[index];
    }

    double enorm = planck_luminosity(1.015);
    if (counts[0]) {
      enorm /= counts[0];
    }
    for (unsigned int i = 0; i < 100; ++i) {
      double nu = 1. + (i + 0.5) * 0.03;
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
        double rand_freq = UnitConverter< QUANTITY_FREQUENCY >::to_unit(
                               spectrum.get_random_frequency(), "eV") /
                           13.6;
        unsigned int index = (rand_freq - 1.) * 100. / 3.;
        ++counts[index];
      }

      double enorm = HLyc_luminosity(cross_sections, T, 1.045);
      if (counts[1]) {
        enorm /= counts[1];
      }
      for (unsigned int i = 0; i < 100; ++i) {
        double nu = 1. + (i + 0.5) * 0.03;
        assert_values_equal_tol(HLyc_luminosity(cross_sections, T, nu),
                                counts[i] * enorm, 0.1);
      }
    }
  }

  // HeliumLymanContinuumSpectrum
  {
    VernerCrossSections cross_sections;
    HeliumLymanContinuumSpectrum spectrum(cross_sections);
    for (unsigned int iT = 0; iT < 10; ++iT) {
      double T = 1500. + (iT + 0.5) * 13500. / 10;
      spectrum.set_temperature(T);

      unsigned int counts[100];
      for (unsigned int i = 0; i < 100; ++i) {
        counts[i] = 0;
      }
      unsigned int numsample = 1000000;
      for (unsigned int i = 0; i < numsample; ++i) {
        double rand_freq = UnitConverter< QUANTITY_FREQUENCY >::to_unit(
                               spectrum.get_random_frequency(), "eV") /
                           13.6;
        unsigned int index = (rand_freq - 1.) * 100. / 3.;
        ++counts[index];
      }

      double enorm = HeLyc_luminosity(cross_sections, T, 1.825);
      // we can obviously not normalize on the lowest frequency, as the spectrum
      // is zero below frequencies of 1.81
      if (counts[27]) {
        enorm /= counts[27];
      }
      for (unsigned int i = 0; i < 100; ++i) {
        double nu = 1. + (i + 0.5) * 0.03;
        assert_values_equal_tol(HeLyc_luminosity(cross_sections, T, nu),
                                counts[i] * enorm, 0.1);
      }
    }
  }

  // HeliumTwoPhotonContinuumSpectrum
  {
    HeliumTwoPhotonContinuumSpectrum spectrum;
    vector< double > yHe2q;
    vector< double > AHe2q;
    spectrum.get_spectrum(yHe2q, AHe2q);

    unsigned int counts[100];
    for (unsigned int i = 0; i < 100; ++i) {
      counts[i] = 0;
    }
    unsigned int numsample = 1000000;
    for (unsigned int i = 0; i < numsample; ++i) {
      double rand_freq = UnitConverter< QUANTITY_FREQUENCY >::to_unit(
                             spectrum.get_random_frequency(), "eV") /
                         13.6;
      unsigned int index = (rand_freq - 1.) * 100. / 0.6;
      ++counts[index];
    }

    double enorm = spectrum.get_integral(yHe2q, AHe2q) / numsample / 0.006;
    for (unsigned int i = 0; i < 100; ++i) {
      double nu = 1. + (i + 0.5) * 0.006;
      assert_values_equal_tol(He2pc_luminosity(yHe2q, AHe2q, nu),
                              counts[i] * enorm, 0.1);
    }
  }

  return 0;
}
