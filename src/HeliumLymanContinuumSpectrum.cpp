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
 * @file HeliumLymanContinuumSpectrum.cpp
 *
 * @brief HeliumLymanContinuumSpectrum implementation.
 *
 * We do not care about units inside this class, and just use the internal units
 * as they were used in Kenny's code. However, we do convert the frequency when
 * it leaves the class, so that the outer world does not need to know about our
 * strange unit system.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "HeliumLymanContinuumSpectrum.hpp"
#include "CrossSections.hpp"
#include "ElementNames.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Fills the precalculated tables.
 *
 * @param cross_sections Photoionization cross sections.
 * @param random_generator RandomGenerator used to generate random numbers.
 */
HeliumLymanContinuumSpectrum::HeliumLymanContinuumSpectrum(
    CrossSections &cross_sections, RandomGenerator &random_generator)
    : _random_generator(random_generator) {
  // 24.6 eV in Hz (1.81 x 13.6 eV)
  const double min_frequency = 1.81 * 3.288465385e15;
  // 54.4 eV in Hz
  const double max_frequency = 4. * 3.288465385e15;
  // Planck constant (in J s)
  const double planck_constant = 6.626e-34;
  // Boltzmann constant (in J s^-1)
  const double boltzmann_constant = 1.38e-23;
  // set up the frequency bins
  for (unsigned int i = 0; i < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _frequency[i] = min_frequency +
                    i * (max_frequency - min_frequency) /
                        (HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ - 1.);
  }
  for (unsigned int iT = 0; iT < HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP; ++iT) {
    _cumulative_distribution[iT][0] = 0.;
    _temperature[iT] =
        1500. + (iT + 0.5) * 13500. / HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP;
    for (unsigned int inu = 1; inu < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      double xsecHe =
          cross_sections.get_cross_section(ION_He_n, _frequency[inu - 1]);
      double jHeIi1 =
          _frequency[inu - 1] * _frequency[inu - 1] * _frequency[inu - 1] *
          xsecHe *
          std::exp(-(planck_constant * (_frequency[inu - 1] - min_frequency)) /
                   (boltzmann_constant * _temperature[iT]));
      xsecHe = cross_sections.get_cross_section(ION_He_n, _frequency[inu]);
      double jHeIi2 =
          _frequency[inu] * _frequency[inu] * _frequency[inu] * xsecHe *
          std::exp(-(planck_constant * (_frequency[inu] - min_frequency)) /
                   (boltzmann_constant * _temperature[iT]));
      _cumulative_distribution[iT][inu] =
          0.5 * (jHeIi1 / _frequency[inu] + jHeIi2 / _frequency[inu - 1]) *
          (_frequency[inu] - _frequency[inu - 1]);
    }
    // make cumulative
    for (unsigned int inu = 1; inu < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] =
          _cumulative_distribution[iT][inu - 1] +
          _cumulative_distribution[iT][inu];
    }
    // normalize
    for (unsigned int inu = 0; inu < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] /=
          _cumulative_distribution[iT]
                                  [HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ - 1];
    }
  }

  _current_T = 0.;
  _current_iT = 0;
}

/**
 * @brief Set the current temperature for the interpolation.
 *
 * @param T New value for the temperature (in K).
 */
void HeliumLymanContinuumSpectrum::set_temperature(double T) {
  _current_T = T;
  _current_iT = Utilities::locate(_current_T, _temperature,
                                  HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP);
}

/**
 * @brief Sample a random frequency from the spectrum.
 *
 * @return Random frequency (in Hz).
 */
double HeliumLymanContinuumSpectrum::get_random_frequency() {
  double x = _random_generator.get_uniform_random_double();
  unsigned int inu1 =
      Utilities::locate(x, _cumulative_distribution[_current_iT],
                        HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ);
  unsigned int inu2 =
      Utilities::locate(x, _cumulative_distribution[_current_iT + 1],
                        HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ);
  double frequency =
      _frequency[inu1] +
      (_current_T - _temperature[_current_iT]) *
          (_frequency[inu2] - _frequency[inu1]) /
          (_temperature[_current_iT + 1] - _temperature[_current_iT]);
  return frequency;
}
