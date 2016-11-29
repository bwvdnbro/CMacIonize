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
 * @file HeliumTwoPhotonContinuumSpectrum.cpp
 *
 * @brief HeliumTwoPhotonContinuumSpectrum implementation.
 *
 * We do not care about units inside this class, and just use the internal units
 * as they were used in Kenny's code. However, we do convert the frequency when
 * it leaves the class, so that the outer world does not need to know about our
 * strange unit system.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "HeliumTwoPhotonContinuumSpectrum.hpp"
#include "HeliumTwoPhotonContinuumDataLocation.hpp"
#include "Utilities.hpp"
#include <fstream>
#include <vector>
using namespace std;

/**
 * @brief Constructor.
 *
 * Reads in the data file containing the spectrum and pre-calculates the
 * cumulative distribution used for random sampling.
 *
 * @param random_generator RandomGenerator used to generate random numbers.
 */
HeliumTwoPhotonContinuumSpectrum::HeliumTwoPhotonContinuumSpectrum(
    RandomGenerator &random_generator)
    : _random_generator(random_generator) {
  vector< double > yHe2q;
  vector< double > AHe2q;
  get_spectrum(yHe2q, AHe2q);

  // 13.6 eV in Hz
  const double min_frequency = 3.288465385e15;
  const double max_frequency = 1.6 * min_frequency;
  const double nu0 = 4.98e15;
  for (unsigned int i = 0; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _frequency[i] = min_frequency +
                    i * (max_frequency - min_frequency) /
                        (HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ - 1.);
  }
  _cumulative_distribution[0] = 0.;
  for (unsigned int i = 1; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    double y1 = _frequency[i - 1] / nu0;
    double AHe2q1 = 0.;
    if (y1 < 1.) {
      unsigned int iHe1 = Utilities::locate(y1, &yHe2q[0], 41);
      double f = (y1 - yHe2q[iHe1]) / (yHe2q[iHe1 + 1] - yHe2q[iHe1]);
      AHe2q1 = AHe2q[iHe1] + f * (AHe2q[iHe1 + 1] - AHe2q[iHe1]);
    }
    double AHe2q2 = 0.;
    double y2 = _frequency[i] / nu0;
    if (y2 < 1.) {
      unsigned int iHe2 = Utilities::locate(y2, &yHe2q[0], 41);
      double f = (y2 - yHe2q[iHe2]) / (yHe2q[iHe2 + 1] - yHe2q[iHe2]);
      AHe2q2 = AHe2q[iHe2] + f * (AHe2q[iHe2 + 1] - AHe2q[iHe2]);
    }
    _cumulative_distribution[i] =
        0.5 * (AHe2q1 + AHe2q2) * (_frequency[i] - _frequency[i - 1]);
  }
  for (unsigned int i = 1; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] =
        _cumulative_distribution[i - 1] + _cumulative_distribution[i];
  }
  for (unsigned int i = 0; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] /=
        _cumulative_distribution[HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ - 1];
  }
}

/**
 * @brief Read the spectrum from the data file and store it in the given
 * vectors.
 *
 * We use this approach to replicate the spectrum in the unit test, since we do
 * not want to store the actual spectrum in this class.
 *
 * @param yHe2q std::vector to store the y values in.
 * @param AHe2q std::vector to store the A values in.
 */
void HeliumTwoPhotonContinuumSpectrum::get_spectrum(
    std::vector< double > &yHe2q, std::vector< double > &AHe2q) {
  ifstream ifile(HELIUMTWOPHOTONCONTINUUMDATALOCATION);
  yHe2q.resize(41);
  AHe2q.resize(41);
  for (unsigned int i = 0; i < 41; ++i) {
    ifile >> yHe2q[i] >> AHe2q[i];
  }
}

/**
 * @brief Get the integral under the curve of the tabulated spectrum.
 *
 * @param yHe2q y values of the spectrum.
 * @param AHe2q A values of the spectrum.
 * @return Integral under the curve in frequency space.
 */
double
HeliumTwoPhotonContinuumSpectrum::get_integral(std::vector< double > &yHe2q,
                                               std::vector< double > &AHe2q) {
  double integral = 0.;
  for (unsigned int i = 1; i < 41; ++i) {
    if (yHe2q[i - 1] > 3.289e15 / 4.98e15) {
      integral += 0.5 * (AHe2q[i - 1] + AHe2q[i]) * (yHe2q[i] - yHe2q[i - 1]);
    } else {
      if (yHe2q[i] > 3.289e15 / 4.98e15) {
        integral += AHe2q[i] * (yHe2q[i] - 3.289e15 / 4.98e15);
      }
    }
  }
  return integral * 4.98e15 / 3.289e15;
}

/**
 * @brief Get a random frequency distributed according to the spectrum.
 *
 * @return Random frequency (in Hz).
 */
double HeliumTwoPhotonContinuumSpectrum::get_random_frequency() {
  double x = _random_generator.get_uniform_random_double();
  unsigned int inu = Utilities::locate(
      x, _cumulative_distribution, HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ);
  double frequency =
      _frequency[inu] +
      (_frequency[inu + 1] - _frequency[inu]) *
          (x - _cumulative_distribution[inu]) /
          (_cumulative_distribution[inu + 1] - _cumulative_distribution[inu]);
  return frequency;
}
