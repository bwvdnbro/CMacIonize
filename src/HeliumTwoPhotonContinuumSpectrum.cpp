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
#include <string>
#include <vector>

/**
 * @brief Constructor.
 *
 * Reads in the data file containing the spectrum and pre-calculates the
 * cumulative distribution used for random sampling.
 */
HeliumTwoPhotonContinuumSpectrum::HeliumTwoPhotonContinuumSpectrum() {

  // allocate memory for the data tables
  _frequency.resize(HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ, 0.);
  _cumulative_distribution.resize(HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ, 0.);

  std::vector< double > yHe2q;
  std::vector< double > AHe2q;
  get_spectrum(yHe2q, AHe2q);

  // 13.6 eV in Hz
  const double min_frequency = 3.288465385e15;
  const double max_frequency = 1.6 * min_frequency;
  const double nu0 = 4.98e15;
  for (uint_fast32_t i = 0; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _frequency[i] = min_frequency +
                    i * (max_frequency - min_frequency) /
                        (HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ - 1.);
  }
  _cumulative_distribution[0] = 0.;
  // NOTE that we compute every y1/y2 twice (except the first and last bin edge)
  // it would be fairly easy to make this more efficient, but since this routine
  // is only called once at the start of the program, we don't bother
  for (uint_fast32_t i = 1; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    // get the y value for the lower edge of the frequency bin
    const double y1 = _frequency[i - 1] / nu0;
    // find the corresponding rate in the energy distribution by linear
    // interpolatino on the tabulated spectrum
    double AHe2q1 = 0.;
    if (y1 < 1.) {
      const uint_fast32_t iHe1 = Utilities::locate(y1, &yHe2q[0], 41);
      const double f = (y1 - yHe2q[iHe1]) / (yHe2q[iHe1 + 1] - yHe2q[iHe1]);
      AHe2q1 = AHe2q[iHe1] + f * (AHe2q[iHe1 + 1] - AHe2q[iHe1]);
    }
    // now do the same for the upper edge of the frequency bin
    const double y2 = _frequency[i] / nu0;
    double AHe2q2 = 0.;
    if (y2 < 1.) {
      const uint_fast32_t iHe2 = Utilities::locate(y2, &yHe2q[0], 41);
      const double f = (y2 - yHe2q[iHe2]) / (yHe2q[iHe2 + 1] - yHe2q[iHe2]);
      AHe2q2 = AHe2q[iHe2] + f * (AHe2q[iHe2 + 1] - AHe2q[iHe2]);
    }
    // set the value for the spectrum in the bin using a simple first order
    // quadrature rule
    _cumulative_distribution[i] =
        0.5 * (AHe2q1 + AHe2q2) * (_frequency[i] - _frequency[i - 1]);
  }
  // now make the spectrum cumulative...
  for (uint_fast32_t i = 1; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] =
        _cumulative_distribution[i - 1] + _cumulative_distribution[i];
  }
  // ...and normalize it (the last entry in the table contains the total sum)
  for (uint_fast32_t i = 0; i < HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ; ++i) {
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
    std::vector< double > &yHe2q, std::vector< double > &AHe2q) const {

  yHe2q.resize(41);
  AHe2q.resize(41);
  std::ifstream ifile(HELIUMTWOPHOTONCONTINUUMDATALOCATION);
  // skip comment line
  std::string line;
  std::getline(ifile, line);
  for (uint_fast8_t i = 0; i < 41; ++i) {
    ifile >> yHe2q[i] >> AHe2q[i];
  }
}

/**
 * @brief Get the integral under the curve of the tabulated spectrum.
 *
 * This is used to normalize the randomly sampled distribution in the unit test.
 *
 * @param yHe2q y values of the spectrum.
 * @param AHe2q A values of the spectrum.
 * @return Integral under the curve in frequency space.
 */
double HeliumTwoPhotonContinuumSpectrum::get_integral(
    std::vector< double > &yHe2q, std::vector< double > &AHe2q) const {

  // we use a simple linear quadrature rule
  const double miny = 3.289e15 / 4.98e15;
  double integral = 0.;
  for (uint_fast8_t i = 1; i < 41; ++i) {
    // the spectrum is cut off at miny; we do not take into account the part of
    // the spectrum below that frequency
    if (yHe2q[i - 1] > miny) {
      integral += 0.5 * (AHe2q[i - 1] + AHe2q[i]) * (yHe2q[i] - yHe2q[i - 1]);
    } else {
      if (yHe2q[i] > miny) {
        integral += AHe2q[i] * (yHe2q[i] - miny);
      }
    }
  }
  return integral * 4.98e15 / 3.289e15;
}

/**
 * @brief Get a random frequency distributed according to the spectrum.
 *
 * We first sample a random uniform number, and locate the position of that
 * number in the normalized cumulative distribution table. We then linearly
 * interpolate within the bin that contains the random number.
 *
 * @param random_generator RandomGenerator to use.
 * @param temperature Temperature of the cell that reemits the photon (in K).
 * @return Random frequency (in Hz).
 */
double HeliumTwoPhotonContinuumSpectrum::get_random_frequency(
    RandomGenerator &random_generator, double temperature) const {

  const double x = random_generator.get_uniform_random_double();
  const uint_fast32_t inu =
      Utilities::locate(x, _cumulative_distribution.data(),
                        HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ);
  const double frequency =
      _frequency[inu] +
      (_frequency[inu + 1] - _frequency[inu]) *
          (x - _cumulative_distribution[inu]) /
          (_cumulative_distribution[inu + 1] - _cumulative_distribution[inu]);
  return frequency;
}

/**
 * @brief Get the total ionizing flux of the spectrum.
 *
 * @warning This method is currently not used and therefore not implemented.
 *
 * @return Total ionizing flux (in m^-2 s^-1).
 */
double HeliumTwoPhotonContinuumSpectrum::get_total_flux() const {
  cmac_error(
      "HeliumTwoPhotonContinuumSpectrum::get_total_flux() is not implemented!");
  return 0.;
}
