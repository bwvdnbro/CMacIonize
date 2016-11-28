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
 * We do not care about units inside this class, and just use the internal units
 * as they were used in Kenny's code. However, we do convert the frequency when
 * it leaves the class, so that the outer world does not need to know about our
 * strange unit system.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "PlanckPhotonSourceSpectrum.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Sets up the internal arrays used for random sampling.
 *
 * @param random_generator RandomGenerator used to generate random numbers.
 * @param temperature Temperature of the black body (in K).
 * @param log Log to write logging info to.
 */
PlanckPhotonSourceSpectrum::PlanckPhotonSourceSpectrum(
    RandomGenerator &random_generator, double temperature, Log *log)
    : _random_generator(random_generator) {
  // some constants
  const double max_frequency = 4.;
  const double min_frequency = 3.289e15;
  const double planck_constant = 6.626e-27;
  const double boltzmann_constant = 1.38e-16;
  // set up the frequency bins and calculate the Planck luminosities
  for (unsigned int i = 0; i < PLANCKPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _frequency[i] =
        1. +
        i * (max_frequency - 1.) / (PLANCKPHOTONSOURCESPECTRUM_NUMFREQ - 1.);
    _luminosity[i] = _frequency[i] * _frequency[i] * _frequency[i] /
                     (exp(planck_constant * _frequency[i] * min_frequency /
                          (boltzmann_constant * temperature)) -
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

  if (log) {
    log->write_status("Set up a Planck black body spectrum with temperature ",
                      temperature, " K.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param random_generator RandomGenerator.
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
PlanckPhotonSourceSpectrum::PlanckPhotonSourceSpectrum(
    RandomGenerator &random_generator, ParameterFile &params, Log *log)
    : PlanckPhotonSourceSpectrum(
          random_generator, params.get_physical_value< QUANTITY_TEMPERATURE >(
                                "photonsourcespectrum.temperature", "40000 K"),
          log) {}

/**
 * @brief Get a random frequency from a Planck blackbody spectrum.
 *
 * @return Random frequency (in Hz).
 */
double PlanckPhotonSourceSpectrum::get_random_frequency() {
  double x = _random_generator.get_uniform_random_double();

  unsigned int ix = Utilities::locate(x, _cumulative_distribution,
                                      PLANCKPHOTONSOURCESPECTRUM_NUMFREQ);
  double log_random_frequency =
      (log10(x) - _log_cumulative_distribution[ix]) /
          (_log_cumulative_distribution[ix + 1] -
           _log_cumulative_distribution[ix]) *
          (_log_frequency[ix + 1] - _log_frequency[ix]) +
      _log_frequency[ix];
  double frequency = pow(10., log_random_frequency);
  // we manually convert from 13.6 eV to Hz, since the UnitConverter is too
  // slow (and binning the actual frequencies in Hz yields a bad interpolation)
  return frequency * 3.288465385e15;
}
