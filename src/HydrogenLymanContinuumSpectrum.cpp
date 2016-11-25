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
 * We do not care about units inside this class, and just use the internal units
 * as they were used in Kenny's code. However, we do convert the frequency when
 * it leaves the class, so that the outer world does not need to know about our
 * strange unit system.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "HydrogenLymanContinuumSpectrum.hpp"
#include "CrossSections.hpp"
#include "ElementNames.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"
#include <cmath>
using namespace std;

/**
 * @brief Constructor.
 *
 * Fills the precalculated tables.
 *
 * @param cross_sections Photoionization cross sections.
 * @param random_generator RandomGenerator used to generate random numbers.
 */
HydrogenLymanContinuumSpectrum::HydrogenLymanContinuumSpectrum(
    CrossSections &cross_sections, RandomGenerator &random_generator)
    : _random_generator(random_generator) {
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
          ION_H_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(
                       _frequency[inu - 1] * 13.6, "eV"));
      xsecH = UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(xsecH, "cm^2") /
              1.e-18;
      double jHIi1 =
          _frequency[inu - 1] * _frequency[inu - 1] * _frequency[inu - 1] *
          xsecH *
          std::exp(-157919.667 * (_frequency[inu - 1] - 1.) / _temperature[iT]);
      xsecH = cross_sections.get_cross_section(
          ION_H_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(
                       _frequency[inu] * 13.6, "eV"));
      xsecH = UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(xsecH, "cm^2") /
              1.e-18;
      double jHIi2 =
          _frequency[inu] * _frequency[inu] * _frequency[inu] * xsecH *
          std::exp(-157919.667 * (_frequency[inu] - 1.) / _temperature[iT]);
      _cumulative_distribution[iT][inu] =
          0.5 * (jHIi1 / _frequency[inu] + jHIi2 / _frequency[inu - 1]) *
          (_frequency[inu] - _frequency[inu - 1]);
    }
    // make cumulative
    for (unsigned int inu = 1; inu < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] =
          _cumulative_distribution[iT][inu - 1] +
          _cumulative_distribution[iT][inu] * 1.e25;
    }
    // normalize
    for (unsigned int inu = 0; inu < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] /=
          _cumulative_distribution[iT]
                                  [HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ - 1];
    }
  }
}

/**
 * @brief Set the current temperature for the interpolation.
 *
 * @param T New value for the temperature (in K).
 */
void HydrogenLymanContinuumSpectrum::set_temperature(double T) {
  _current_T = T;
}

/**
 * @brief Get a random frequency from the spectrum.
 *
 * @return Random frequency (in Hz).
 */
double HydrogenLymanContinuumSpectrum::get_random_frequency() {
  unsigned int iT = Utilities::locate(_current_T, _temperature,
                                      HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP);
  double x = _random_generator.get_uniform_random_double();
  unsigned int inu1 = Utilities::locate(x, _cumulative_distribution[iT],
                                        HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ);
  unsigned int inu2 = Utilities::locate(x, _cumulative_distribution[iT + 1],
                                        HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ);
  double frequency = _frequency[inu1] +
                     (_current_T - _temperature[iT]) *
                         (_frequency[inu2] - _frequency[inu1]) /
                         (_temperature[iT + 1] - _temperature[iT]);
  return UnitConverter::to_SI< QUANTITY_FREQUENCY >(13.6 * frequency, "eV");
}
