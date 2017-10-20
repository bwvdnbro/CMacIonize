/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file MaskedPhotonSourceSpectrum.cpp
 *
 * @brief MaskedPhotonSourceSpectrum implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "MaskedPhotonSourceSpectrum.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "PhotonSourceSpectrumMaskFactory.hpp"
#include "RandomGenerator.hpp"
#include "Utilities.hpp"

/**
 * @brief Constructor.
 *
 * @param unmasked_spectrum Unmasked PhotonSpectrum.
 * @param mask PhotonSourceSpectrumMask to apply.
 * @param number_of_bins Number of bins to use to bin the spectrum.
 * @param number_of_samples Number of random samples to take to sample the
 * unmasked spectrum.
 */
MaskedPhotonSourceSpectrum::MaskedPhotonSourceSpectrum(
    PhotonSourceSpectrum *unmasked_spectrum, PhotonSourceSpectrumMask *mask,
    const uint_fast16_t number_of_bins, const uint_fast32_t number_of_samples) {

  // allocate memory for frequency arrays
  _frequency_bins.resize(number_of_bins, 0.);
  _cumulative_distribution.resize(number_of_bins, 0.);

  // construct the frequency bins
  // 13.6 eV in Hz
  const double min_frequency = 3.289e15;
  const double max_frequency = 4. * min_frequency;
  const double frequency_bin_size =
      (max_frequency - min_frequency) / (number_of_bins - 1.);
  for (uint_fast16_t i = 0; i < number_of_bins; ++i) {
    _frequency_bins[i] = min_frequency + i * frequency_bin_size;
  }

  // now sample the masked distribution to get the masked spectrum
  RandomGenerator random_generator;
  for (uint_fast32_t i = 0; i < number_of_samples; ++i) {
    const double random_frequency =
        unmasked_spectrum->get_random_frequency(random_generator);
    const uint_fast16_t index =
        (random_frequency - min_frequency) / frequency_bin_size;
    _cumulative_distribution[index] += 1.;
  }

  // apply mask
  for (uint_fast16_t i = 0; i < number_of_bins; ++i) {
    _cumulative_distribution[i] *= mask->get_bin_fraction(_frequency_bins[i]);
  }

  // make cumulative
  for (uint_fast16_t i = 1; i < number_of_bins; ++i) {
    _cumulative_distribution[i] += _cumulative_distribution[i - 1];
  }

  // normalize
  const double norm = _cumulative_distribution.back();
  const double norm_inv = 1. / norm;
  for (uint_fast16_t i = 0; i < number_of_bins; ++i) {
    _cumulative_distribution[i] *= norm_inv;
  }

  // apply the mask to the total ionizing flux
  _ionizing_flux =
      norm * unmasked_spectrum->get_total_flux() / number_of_samples;

  delete unmasked_spectrum;
  delete mask;
}

/**
 * @brief ParameterFile constructor.
 *
 * @param role Role the masked spectrum will fulfil in the simulation.
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
MaskedPhotonSourceSpectrum::MaskedPhotonSourceSpectrum(std::string role,
                                                       ParameterFile &params,
                                                       Log *log)
    : MaskedPhotonSourceSpectrum(
          PhotonSourceSpectrumFactory::generate_from_type(
              params.get_value< std::string >(role + ":masked type", "Planck"),
              role, params, log),
          PhotonSourceSpectrumMaskFactory::generate(role, params, log),
          params.get_value< uint_fast16_t >(role + ":mask number of bins",
                                            1000),
          params.get_value< uint_fast32_t >(role + ":mask number of samples",
                                            1e7)) {}

/**
 * @brief Get a random frequency from the masked spectrum.
 *
 * @param random_generator RandomGenerator to use to generate pseudo random
 * numbers.
 * @param temperature Temperature at which to evaluate the spectrum (in K,
 * ignored).
 * @return Random frequency distributed according to the spectrum (in Hz).
 */
double MaskedPhotonSourceSpectrum::get_random_frequency(
    RandomGenerator &random_generator, double temperature) const {

  const double x = random_generator.get_uniform_random_double();
  const uint_fast32_t inu = Utilities::locate(
      x, _cumulative_distribution.data(), _frequency_bins.size());
  const double frequency =
      _frequency_bins[inu] +
      (_frequency_bins[inu + 1] - _frequency_bins[inu]) *
          (x - _cumulative_distribution[inu]) /
          (_cumulative_distribution[inu + 1] - _cumulative_distribution[inu]);
  return frequency;
}

/**
 * @brief Get the total ionizing flux of the masked spectrum.
 *
 * @return Total ionizing flux of the spectrum (in m^-2 s^-1).
 */
double MaskedPhotonSourceSpectrum::get_total_flux() const {
  return _ionizing_flux;
}

/**
 * @brief Get the masked spectrum.
 *
 * @return Frequency bins and cumulative spectrum.
 */
std::pair< std::vector< double >, std::vector< double > >
MaskedPhotonSourceSpectrum::get_spectrum() const {

  std::vector< double > spectrum = _cumulative_distribution;
  // retrieve original (non cumulative, non normalized) spectrum
  const uint_fast16_t number_of_bins = _frequency_bins.size();
  for (uint_fast16_t i = number_of_bins; i > 0; --i) {
    spectrum[i - 1] -= spectrum[i - 2];
    spectrum[i - 1] *= _ionizing_flux;
  }
  return std::make_pair(_frequency_bins, spectrum);
}
