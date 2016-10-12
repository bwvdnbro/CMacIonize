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
 * @file HydrogenLymanContinuumSpectrum.cpp
 *
 * @brief HydrogenLymanContinuumSpectrum implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "HydrogenLymanContinuumSpectrum.hpp"
#include "CrossSections.hpp"
#include "ElementNames.hpp"
#include "Utilities.hpp"
#include <cmath>
using namespace std;

/**
 * @brief Constructor.
 *
 * Fills the precalculated tables.
 *
 * @param cross_sections Photoionization cross sections.
 */
HydrogenLymanContinuumSpectrum::HydrogenLymanContinuumSpectrum(
    CrossSections &cross_sections) {
  double max_frequency = 4.;
  // set up the frequency bins
  for (unsigned int i = 0; i < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _frequency[i] = 1. +
                    i * (max_frequency - 1.) /
                        (HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ - 1.);
  }
  for (unsigned int iT = 0; iT < HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP; ++iT) {
    _cumulative_distribution[iT][0] = 0.;
    _temperature[iT] =
        1500. + (iT + 0.5) * 13500. / HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP;
    for (unsigned int inu = 1; inu < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      double xsecH = cross_sections.get_cross_section(
          ELEMENT_H, _frequency[inu - 1] * 13.6);
      double jHIi1 =
          _frequency[inu - 1] * _frequency[inu - 1] * _frequency[inu - 1] *
          xsecH *
          exp(-157919.667 * (_frequency[inu - 1] - 1.) / _temperature[iT]);
      xsecH =
          cross_sections.get_cross_section(ELEMENT_H, _frequency[inu] * 13.6);
      double jHIi2 =
          _frequency[inu] * _frequency[inu] * _frequency[inu] * xsecH *
          exp(157919.667 * (_frequency[inu] - 1.) / _temperature[iT]);
      _cumulative_distribution[iT][inu] =
          _cumulative_distribution[iT][inu - 1] +
          0.5 * (jHIi1 / _frequency[inu] + jHIi2 / _frequency[inu - 1]) *
              (_frequency[inu] - _frequency[inu - 1]);
    }
    // normalize
    for (unsigned int inu = 0; inu < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] /=
          _cumulative_distribution[HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP - 1]
                                  [inu];
    }
  }
}

/**
 * @brief Set the current temperature for the interpolation.
 *
 * @param New value for the temperature.
 */
void HydrogenLymanContinuumSpectrum::set_temperature(double T) {
  _current_T = T;
}

/**
 * @brief Get a random frequency from the spectrum.
 *
 * @return Random frequency.
 */
double HydrogenLymanContinuumSpectrum::get_random_frequency() {
  unsigned int iT = Utilities::locate(_current_T, _temperature,
                                      HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP);
  double x = Utilities::random_double();
  unsigned int inu = Utilities::locate(x, _cumulative_distribution[iT],
                                       HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ);
  return _frequency[inu];
}
