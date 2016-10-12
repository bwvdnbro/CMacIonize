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
 * @file PlanckPhotonSourceSpectrum.cpp
 *
 * @brief PhotonSourceSpectrum implementation for a Planck blackbody spectrum:
 * implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "PlanckPhotonSourceSpectrum.hpp"
#include "Error.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Sets up the internal arrays used for random sampling.
 */
PlanckPhotonSourceSpectrum::PlanckPhotonSourceSpectrum() {
  // some constants
  double max_frequency = 4.;
  double min_frequency = 3.289e15;
  double planck_constant = 6.626e-27;
  double boltzmann_constant = 1.38e-16;
  // not a constant! Should be replaced by a parameter!
  double temperature_star = 40000.;
  // set up the frequency bins and calculate the Planck luminosities
  for (unsigned int i = 0; i < PLANCKPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _frequency[i] =
        1. +
        i * (max_frequency - 1.) / (PLANCKPHOTONSOURCESPECTRUM_NUMFREQ - 1.);
    _luminosity[i] = _frequency[i] * _frequency[i] * _frequency[i] /
                     (exp(planck_constant * _frequency[i] * min_frequency /
                          (boltzmann_constant * temperature_star)) -
                      1.);
  }

  // convert the Planck luminosities to a cumulative distribution
  _cumulative_distribution[0] = 0.;
  for (unsigned int i = 1; i < PLANCKPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] = _cumulative_distribution[i - 1] +
                                  0.5 *
                                      (_luminosity[i] / _frequency[i] +
                                       _luminosity[i - 1] / _frequency[i - 1]) *
                                      (_frequency[i] - _frequency[i - 1]);
  }

  // normalize the cumulative distribution and calculate the logarithms
  _log_cumulative_distribution[0] = -10.;
  _log_frequency[0] = 0.;
  for (unsigned int i = 1; i < PLANCKPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] /=
        _cumulative_distribution[PLANCKPHOTONSOURCESPECTRUM_NUMFREQ - 1];
    _log_cumulative_distribution[i] = log10(_cumulative_distribution[i]);
    _log_frequency[i] = log10(_frequency[i]);
  }
}

/**
 * @brief Get a random frequency from a Planck blackbody spectrum.
 *
 * @return Random frequency.
 */
double PlanckPhotonSourceSpectrum::get_random_frequency() {
  double x = Utilities::random_double();

  unsigned int ix = Utilities::locate(x, _cumulative_distribution,
                                      PLANCKPHOTONSOURCESPECTRUM_NUMFREQ);
  double log_random_frequency =
      (log10(x) - _log_cumulative_distribution[ix]) /
          (_log_cumulative_distribution[ix + 1] -
           _log_cumulative_distribution[ix]) *
          (_log_frequency[ix + 1] - _log_frequency[ix]) +
      _log_frequency[ix];
  return pow(10., log_random_frequency);
}
