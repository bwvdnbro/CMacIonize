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
#include "PhysicalConstants.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Fills the precalculated tables.
 *
 * @param cross_sections Photoionization cross sections.
 */
HydrogenLymanContinuumSpectrum::HydrogenLymanContinuumSpectrum(
    const CrossSections &cross_sections) {

  // allocate memory for the data tables
  _frequency.resize(HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ, 0.);
  _temperature.resize(HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP, 0.);
  _cumulative_distribution.resize(
      HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP,
      std::vector< double >(HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ, 0.));

  // some constants
  // 13.6 eV in Hz
  const double min_frequency = 3.289e15;
  // 54.4 eV in Hz
  const double max_frequency = 4. * min_frequency;
  // Planck constant (in J s)
  const double planck_constant =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  // Boltzmann constant (in J s^-1)
  const double boltzmann_constant =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

  // set up the frequency bins
  for (uint_fast32_t i = 0; i < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ; ++i) {
    _frequency[i] = min_frequency +
                    i * (max_frequency - min_frequency) /
                        (HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ - 1.);
  }

  // set up the temperature bins and precompute the spectrum
  for (uint_fast32_t iT = 0; iT < HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP;
       ++iT) {
    _cumulative_distribution[iT][0] = 0.;
    _temperature[iT] =
        1500. + (iT + 0.5) * 13500. / HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP;
    // precompute the spectrum for this temperature
    for (uint_fast32_t inu = 1; inu < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      // first do the lower edge of the frequency interval
      double xsecH =
          cross_sections.get_cross_section(ION_H_n, _frequency[inu - 1]);
      // Wood, Mathis & Ercolano (2004), equation (8)
      // note that we ignore all constant prefactors, since we normalize the
      // spectrum afterwards
      // note that the temperature prefactor is also constant within a
      // temperature table and since we normalize temperature tables, we can
      // also ignore it here
      const double jHIi1 =
          _frequency[inu - 1] * _frequency[inu - 1] * _frequency[inu - 1] *
          xsecH *
          std::exp(-(planck_constant * (_frequency[inu - 1] - min_frequency)) /
                   (boltzmann_constant * _temperature[iT]));
      // now do the upper edge of the interval
      xsecH = cross_sections.get_cross_section(ION_H_n, _frequency[inu]);
      const double jHIi2 =
          _frequency[inu] * _frequency[inu] * _frequency[inu] * xsecH *
          std::exp(-(planck_constant * (_frequency[inu] - min_frequency)) /
                   (boltzmann_constant * _temperature[iT]));
      // the spectrum in the bin is computed using a simple linear quadrature
      // rule
      // we convert the energy spectrum to a number spectrum by dividing by the
      // frequency
      _cumulative_distribution[iT][inu] =
          0.5 * (jHIi1 / _frequency[inu] + jHIi2 / _frequency[inu - 1]) *
          (_frequency[inu] - _frequency[inu - 1]);
    }

    // make cumulative
    for (uint_fast32_t inu = 1; inu < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] =
          _cumulative_distribution[iT][inu - 1] +
          _cumulative_distribution[iT][inu];
    }
    // normalize
    for (uint_fast32_t inu = 0; inu < HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ;
         ++inu) {
      _cumulative_distribution[iT][inu] /=
          _cumulative_distribution[iT]
                                  [HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ - 1];
    }
  }
}

/**
 * @brief Get a random frequency from the spectrum.
 *
 * We first locate the given temperature value in the temperature table. Then we
 * generate a random uniform number and locate it in the two cumulative
 * distributions that border the temperature value. The sampled frequency is
 * then given by linear interpolation on the temperature and frequency tables.
 *
 * @param random_generator RandomGenerator to use.
 * @param temperature Temperature of the cell that reemits the photon (in K).
 * @return Random frequency (in Hz).
 */
double HydrogenLymanContinuumSpectrum::get_random_frequency(
    RandomGenerator &random_generator, double temperature) const {

  const uint_fast32_t iT = Utilities::locate(
      temperature, _temperature.data(), HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP);
  const double x = random_generator.get_uniform_random_double();
  const uint_fast32_t inu1 =
      Utilities::locate(x, _cumulative_distribution[iT].data(),
                        HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ);
  const uint_fast32_t inu2 =
      Utilities::locate(x, _cumulative_distribution[iT + 1].data(),
                        HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ);
  const double frequency = _frequency[inu1] +
                           (temperature - _temperature[iT]) *
                               (_frequency[inu2] - _frequency[inu1]) /
                               (_temperature[iT + 1] - _temperature[iT]);
  return frequency;
}

/**
 * @brief Get the total ionizing flux of the spectrum.
 *
 * @warning This method is currently not used and therefore not implemented.
 *
 * @return Total ionizing flux (in m^-2 s^-1).
 */
double HydrogenLymanContinuumSpectrum::get_total_flux() const {
  cmac_error(
      "HydrogenLymanContinuumSpectrum::get_total_flux() is not implemented!");
  return 0.;
}
