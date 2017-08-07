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
#include "RandomGenerator.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Sets up the internal arrays used for random sampling. Note that we use
 * logarithmic binning for this spectrum due to its high dynamic range in
 * values. This means we need to store both the cumulative distribution function
 * and the logarithm of the cumulative distribution function.
 *
 * @param temperature Temperature of the black body (in K).
 * @param ionizing_flux Ionizing flux of the spectrum (in m^-2 s^-1).
 * @param log Log to write logging info to.
 */
PlanckPhotonSourceSpectrum::PlanckPhotonSourceSpectrum(double temperature,
                                                       double ionizing_flux,
                                                       Log *log)
    : _ionizing_flux(ionizing_flux) {

  _log_frequency.resize(PLANCKPHOTONSOURCESPECTRUM_NUMFREQ, 0.);
  _cumulative_distribution.resize(PLANCKPHOTONSOURCESPECTRUM_NUMFREQ, 0.);
  _log_cumulative_distribution.resize(PLANCKPHOTONSOURCESPECTRUM_NUMFREQ, 0.);

  // some constants
  // in units 13.6 eV (corresponds to 54.4 eV)
  const double max_frequency = 4.;
  // 13.6 eV in Hz
  const double min_frequency = 3.289e15;
  // Planck constant (in J s)
  const double planck_constant = 6.626e-34;
  // Boltzmann constant (in J s^-1)
  const double boltzmann_constant = 1.38e-23;
  // set up the frequency bins and calculate the Planck luminosities
  std::vector< double > frequency(PLANCKPHOTONSOURCESPECTRUM_NUMFREQ, 0.);
  std::vector< double > luminosity(PLANCKPHOTONSOURCESPECTRUM_NUMFREQ, 0.);
  for (unsigned int i = 0; i < PLANCKPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    frequency[i] =
        1. +
        i * (max_frequency - 1.) / (PLANCKPHOTONSOURCESPECTRUM_NUMFREQ - 1.);
    luminosity[i] = frequency[i] * frequency[i] * frequency[i] /
                    (std::exp(planck_constant * frequency[i] * min_frequency /
                              (boltzmann_constant * temperature)) -
                     1.);
  }

  // convert the Planck luminosities to a cumulative distribution
  _cumulative_distribution[0] = 0.;
  for (unsigned int i = 1; i < PLANCKPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] = _cumulative_distribution[i - 1] +
                                  0.5 * (luminosity[i] / frequency[i] +
                                         luminosity[i - 1] / frequency[i - 1]) *
                                      (frequency[i] - frequency[i - 1]);
  }

  // normalize the cumulative distribution and calculate the logarithms
  _log_cumulative_distribution[0] = -10.;
  _log_frequency[0] = 0.;
  for (unsigned int i = 1; i < PLANCKPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] /=
        _cumulative_distribution[PLANCKPHOTONSOURCESPECTRUM_NUMFREQ - 1];
    _log_cumulative_distribution[i] = std::log10(_cumulative_distribution[i]);
    _log_frequency[i] = std::log10(frequency[i]);
  }

  if (log) {
    log->write_status("Set up a Planck black body spectrum with temperature ",
                      temperature, " K.");
    if (_ionizing_flux > 1.) {
      log->write_status("Planck spectrum has a non trivial ionizing flux of ",
                        _ionizing_flux, " m^-2 s^-1.");
    }
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param role Role the spectrum will fulfil in the simulation. Parameters are
 * read from the corresponding block in the parameter file.
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
PlanckPhotonSourceSpectrum::PlanckPhotonSourceSpectrum(std::string role,
                                                       ParameterFile &params,
                                                       Log *log)
    : PlanckPhotonSourceSpectrum(
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              role + ":temperature", "40000. K"),
          params.get_physical_value< QUANTITY_FLUX >(role + ":ionizing_flux",
                                                     "-1. m^-2 s^-1"),
          log) {}

/**
 * @brief Get a random frequency from a Planck blackbody spectrum.
 *
 * We sample a random uniform number and find its location in the cumulative
 * distribution array. We then use the logarithmic arrays to convert this into a
 * frequency.
 *
 * @param random_generator RandomGenerator to use.
 * @param temperature Not used for this spectrum.
 * @return Random frequency (in Hz).
 */
double PlanckPhotonSourceSpectrum::get_random_frequency(
    RandomGenerator &random_generator, double temperature) const {
  double x = random_generator.get_uniform_random_double();

  unsigned int ix = Utilities::locate(x, _cumulative_distribution.data(),
                                      PLANCKPHOTONSOURCESPECTRUM_NUMFREQ);
  double log_random_frequency =
      (std::log10(x) - _log_cumulative_distribution[ix]) /
          (_log_cumulative_distribution[ix + 1] -
           _log_cumulative_distribution[ix]) *
          (_log_frequency[ix + 1] - _log_frequency[ix]) +
      _log_frequency[ix];
  double frequency = std::pow(10., log_random_frequency);
  // we manually convert from 13.6 eV to Hz, since the UnitConverter is too
  // slow (and binning the actual frequencies in Hz yields a bad interpolation)
  return frequency * 3.288465385e15;
}

/**
 * @brief Get the total ionizing flux of the spectrum.
 *
 * @return Total ionizing flux (in m^-2 s^-1).
 */
double PlanckPhotonSourceSpectrum::get_total_flux() const {
  if (_ionizing_flux < 0.) {
    cmac_error("PlanckPhotonSourceSpectrum is used as external spectrum, but "
               "no ionizing flux was provided in the parameter file!");
  }
  return _ionizing_flux;
}
