/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file HolmesPhotonSourceSpectrum.cpp
 *
 * @brief HolmesPhotonSourceSpectrum implementation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "HolmesPhotonSourceSpectrum.hpp"
#include "HolmesDataLocation.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

/**
 * @brief Constructor.
 *
 * Reads in the correct model spectrum and resamples it on the 1000 bin internal
 * frequency grid.
 *
 * @param log Log to write logging info to.
 */
HolmesPhotonSourceSpectrum::HolmesPhotonSourceSpectrum(Log *log) {

  std::ifstream file(HOLMESDATALOCATION);
  if (!file) {
    cmac_error("Error while opening file '%s'. Spectrum does not exist!",
               HOLMESDATALOCATION);
  }

  if (log) {
    log->write_status("Reading spectrum from ", HOLMESDATALOCATION, "...");
  }
  std::string line;
  // skip the first five lines from the file, they contain comments
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  // read the data from the remainder of the file
  std::vector< double > file_frequencies;
  std::vector< double > file_eddington_fluxes;
  uint_fast32_t i = 0;
  const double lightspeed =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED);
  while (std::getline(file, line)) {
    std::stringstream linestream(line);
    double file_frequency, file_flux;
    linestream >> file_frequency >> file_flux;

    // file contains wavelength (in Angstrom), convert to frequency
    file_frequencies.push_back(lightspeed * 1.e10 / file_frequency);
    file_eddington_fluxes.push_back(file_flux);

    ++i;
  }

  // the file contains the inverse of the spectrum
  std::reverse(file_frequencies.begin(), file_frequencies.end());
  std::reverse(file_eddington_fluxes.begin(), file_eddington_fluxes.end());

  // allocate memory for the data tables
  _frequencies.resize(HOLMESPHOTONSOURCESPECTRUM_NUMFREQ, 0.);
  _cumulative_distribution.resize(HOLMESPHOTONSOURCESPECTRUM_NUMFREQ, 0.);

  // construct the frequency bins
  // 13.6 eV in Hz
  const double min_frequency = 3.289e15;
  const double max_frequency = 4. * min_frequency;
  for (uint_fast32_t i = 0; i < HOLMESPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _frequencies[i] =
        min_frequency + i * (max_frequency - min_frequency) /
                            (HOLMESPHOTONSOURCESPECTRUM_NUMFREQ - 1.);
  }

  // create the spectrum in the bin range of interest
  _cumulative_distribution[0] = 0.;
  const size_t num_frequency = file_frequencies.size();
  for (uint_fast32_t i = 1; i < HOLMESPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    double y1 = _frequencies[i - 1];
    uint_fast32_t i1 =
        Utilities::locate(y1, &file_frequencies[0], num_frequency);
    double f = (y1 - file_frequencies[i1]) /
               (file_frequencies[i1 + 1] - file_frequencies[i1]);
    double e1 = file_eddington_fluxes[i1] +
                f * (file_eddington_fluxes[i1 + 1] - file_eddington_fluxes[i1]);
    double y2 = _frequencies[i];
    uint_fast32_t i2 =
        Utilities::locate(y2, &file_frequencies[0], num_frequency);
    f = (y2 - file_frequencies[i2]) /
        (file_frequencies[i2 + 1] - file_frequencies[i2]);
    double e2 = file_eddington_fluxes[i2] +
                f * (file_eddington_fluxes[i2 + 1] - file_eddington_fluxes[i2]);
    _cumulative_distribution[i] = 0.5 * (e1 / y2 + e2 / y1) * (y2 - y1);
  }

  // _cumulative_distribution now contains the actual ionizing spectrum
  // make it cumulative (and at the same time get the total ionizing
  // luminosity)
  for (uint_fast32_t i = 1; i < HOLMESPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] += _cumulative_distribution[i - 1];
  }

  // the total ionizing luminosity is the last element of
  // _cumulative_distribution (using a zeroth order quadrature)
  // its value is in erg Hz^-1 s^-1 cm^-2 sr^-1
  // we convert to s^-1 cm^-2 sr^-1 by dividing by Planck's constant
  // (in J s = J Hz^-1) and multiplying with 10^-7
  _total_flux =
      1.e-7 * _cumulative_distribution[HOLMESPHOTONSOURCESPECTRUM_NUMFREQ - 1] /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  // we integrate out over all solid angles
  _total_flux *= 4. * M_PI;
  // and convert from cm^-2 to m^-2
  _total_flux *= 1.e4;
  // normalize the spectrum
  for (uint_fast32_t i = 0; i < HOLMESPHOTONSOURCESPECTRUM_NUMFREQ; ++i) {
    _cumulative_distribution[i] /=
        _cumulative_distribution[HOLMESPHOTONSOURCESPECTRUM_NUMFREQ - 1];
  }

  if (log) {
    log->write_status(
        "Constructed HolmesPhotonSourceSpectrum with total ionizing flux ",
        _total_flux, " m^-2 s^-1.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * No parameters are read from the parameter file.
 *
 * @param role Role the spectrum will fulfil in the simulation. Parameters are
 * read from the corresponding block in the parameter file.
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
HolmesPhotonSourceSpectrum::HolmesPhotonSourceSpectrum(std::string role,
                                                       ParameterFile &params,
                                                       Log *log)
    : HolmesPhotonSourceSpectrum(log) {}

/**
 * @brief Get a random frequency from a stellar model spectrum.
 *
 * We draw a random uniform number in the range [0, 1] and find the bin in the
 * tabulated cumulative distribution that contains that number. The sampled
 * frequency is then the linearly interpolated value corresponding to that bin.
 *
 * @param random_generator RandomGenerator to use.
 * @param temperature Not used for this spectrum.
 * @return Random frequency (in Hz).
 */
double HolmesPhotonSourceSpectrum::get_random_frequency(
    RandomGenerator &random_generator, double temperature) const {

  const double x = random_generator.get_uniform_random_double();
  const uint_fast32_t inu = Utilities::locate(
      x, _cumulative_distribution.data(), HOLMESPHOTONSOURCESPECTRUM_NUMFREQ);
  const double frequency =
      _frequencies[inu] +
      (_frequencies[inu + 1] - _frequencies[inu]) *
          (x - _cumulative_distribution[inu]) /
          (_cumulative_distribution[inu + 1] - _cumulative_distribution[inu]);
  return frequency;
}

/**
 * @brief Get the total ionizing flux of the spectrum.
 *
 * The total ionizing flux depends on the model spectrum that is used and is
 * computed in the constructor, after the spectrum has been regridded.
 *
 * @return Total ionizing flux (in m^-2 s^-1).
 */
double HolmesPhotonSourceSpectrum::get_total_flux() const {
  cmac_error("Should not be used!");
  return _total_flux;
}
