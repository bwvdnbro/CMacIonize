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
#include "UnitConverter.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Fills the precalculated tables.
 *
 * @param cross_sections Photoionization cross sections.
 */
HeliumLymanContinuumSpectrum::HeliumLymanContinuumSpectrum(
    CrossSections &cross_sections) {
  double max_frequency = 4.;
  // set up the frequency bins
  // there is a hard lower limit on the spectrum, which we manually enforce here
  _frequency[0] = 1.;
  _frequency[1] = 1.81;
  for (unsigned int i = 2; i < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _frequency[i] = 1.81 +
                    i * (max_frequency - 1.81) /
                        (HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ - 3.);
  }
  for (unsigned int iT = 0; iT < HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP; ++iT) {
    _cumulative_distribution[iT][0] = 0.;
    _temperature[iT] =
        1500. + (iT + 0.5) * 13500. / HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP;
    for (unsigned int inu = 1; inu < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      double xsecHe = cross_sections.get_cross_section(
          ELEMENT_He, _frequency[inu - 1] * 13.6);
      double jHeIi1;
      if (_frequency[inu - 1] > 1.81) {
        jHeIi1 = _frequency[inu - 1] * _frequency[inu - 1] *
                 _frequency[inu - 1] * xsecHe *
                 std::exp(-157919.667 * (_frequency[inu - 1] - 1.81) /
                          _temperature[iT]);
      } else {
        jHeIi1 = 0.;
      }
      xsecHe =
          cross_sections.get_cross_section(ELEMENT_He, _frequency[inu] * 13.6);
      double jHeIi2;
      if (_frequency[inu] > 1.81) {
        jHeIi2 =
            _frequency[inu] * _frequency[inu] * _frequency[inu] * xsecHe *
            std::exp(-157919.667 * (_frequency[inu] - 1.81) / _temperature[iT]);
      } else {
        jHeIi2 = 0.;
      }
      _cumulative_distribution[iT][inu] =
          0.5 * (jHeIi1 / _frequency[inu] + jHeIi2 / _frequency[inu - 1]) *
          (_frequency[inu] - _frequency[inu - 1]);
    }
    // make cumulative
    for (unsigned int inu = 1; inu < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] =
          _cumulative_distribution[iT][inu - 1] +
          _cumulative_distribution[iT][inu] * 1.e25;
    }
    // normalize
    for (unsigned int inu = 0; inu < HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] /=
          _cumulative_distribution[iT]
                                  [HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ - 1];
    }
  }
}

/**
 * @brief Set the current temperature for the interpolation.
 *
 * @param T New value for the temperature (in K).
 */
void HeliumLymanContinuumSpectrum::set_temperature(double T) { _current_T = T; }

/**
 * @brief Sample a random frequency from the spectrum.
 *
 * @return Random frequency (in Hz).
 */
double HeliumLymanContinuumSpectrum::get_random_frequency() {
  unsigned int iT = Utilities::locate(_current_T, _temperature,
                                      HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP);
  double x = Utilities::random_double();
  unsigned int inu = Utilities::locate(x, _cumulative_distribution[iT],
                                       HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ);
  if (_frequency[inu] < 1.81) {
    error("This is not right... (%g, %g)", _frequency[inu],
          _frequency[inu + 1]);
  }
  return UnitConverter< QUANTITY_FREQUENCY >::to_SI(13.6 * _frequency[inu],
                                                    "eV");
}
