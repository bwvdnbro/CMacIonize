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
 * @file FaucherGiguerePhotonSourceSpectrum.cpp
 *
 * @brief FaucherGiguerePhotonSourceSpectrum implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "FaucherGiguerePhotonSourceSpectrum.hpp"
#include "Error.hpp"
#include "FaucherGiguereDataLocation.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"
#include "Utilities.hpp"
#include <cmath>
#include <fstream>
#include <sstream>

/**
 * @brief Constructor.
 *
 * @param redshift Redshift for which we want the spectrum.
 * @param log Log to write logging info to.
 */
FaucherGiguerePhotonSourceSpectrum::FaucherGiguerePhotonSourceSpectrum(
    double redshift, Log *log) {

  // allocate memory for the data tables
  _frequencies.resize(FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ, 0.);
  _cumulative_distribution.resize(FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ,
                                  0.);

  // construct the frequency bins
  // 13.6 eV in Hz
  const double min_frequency = 3.289e15;
  const double max_frequency = 4. * min_frequency;
  for (uint_fast32_t i = 0; i < FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ;
       ++i) {
    _frequencies[i] = min_frequency +
                      i * (max_frequency - min_frequency) /
                          (FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ - 1.);
  }

  // find the redshift bin that contains the requested redshift (we have values
  // for all redshift in [0, 10.65], for every delta z = 0.05
  if (redshift <= 10.65) {
    double spectrum_freq[261];
    double spectrum_ener[261];
    const unsigned int izlo = redshift / 0.05;
    const unsigned int izhi = izlo + 1;
    const double zlo = izlo * 0.05;
    const double zhi = izhi * 0.05;
    const double zlo_fac = 20. * (zhi - redshift);
    // open the first file and read in the spectrum
    std::ifstream zlofile(get_filename(zlo));
    std::string line;
    // skip the first two comment lines
    getline(zlofile, line);
    getline(zlofile, line);
    // now read the spectrum
    for (uint_fast16_t i = 0; i < 261; ++i) {
      getline(zlofile, line);
      std::istringstream linestream(line);
      double nu, e;
      linestream >> nu >> e;
      // convert 13.6 eV to Hz
      spectrum_freq[i] = nu * 3.289e15;
      spectrum_ener[i] = zlo_fac * e;
    }
    if (zhi <= 10.65) {
      // we do not execute this part if zhi > 10.65, which can only happen if
      // zlo == z == 10.65
      const double zhi_fac = 20. * (redshift - zlo);
      std::ifstream zhifile(get_filename(zhi));
      // skip the first two comment lines
      getline(zhifile, line);
      getline(zhifile, line);
      // now read the spectrum
      for (uint_fast16_t i = 0; i < 261; ++i) {
        getline(zlofile, line);
        std::istringstream linestream(line);
        double nu, e;
        linestream >> nu >> e;
        spectrum_ener[i] += zhi_fac * e;
      }
    }

    // we now have the interpolated full spectrum
    // create the spectrum in our bin range of interest by linear interpolation
    _cumulative_distribution[0] = 0.;
    for (uint_fast32_t i = 1; i < FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ;
         ++i) {
      const double y1 = _frequencies[i - 1];
      const uint_fast16_t i1 = Utilities::locate(y1, spectrum_freq, 261);
      double f = (y1 - spectrum_freq[i1]) /
                 (spectrum_freq[i1 + 1] - spectrum_freq[i1]);
      const double e1 =
          spectrum_ener[i1] + f * (spectrum_ener[i1 + 1] - spectrum_ener[i1]);
      const double y2 = _frequencies[i];
      const uint_fast16_t i2 = Utilities::locate(y2, spectrum_freq, 261);
      f = (y2 - spectrum_freq[i2]) /
          (spectrum_freq[i2 + 1] - spectrum_freq[i2]);
      const double e2 =
          spectrum_ener[i2] + f * (spectrum_ener[i2 + 1] - spectrum_ener[i2]);
      _cumulative_distribution[i] = 0.5 * (e1 / y2 + e2 / y1) * (y2 - y1);
    }
    // _cumulative_distribution now contains the actual ionizing spectrum
    // make it cumulative (and at the same time get the total ionizing
    // luminosity)
    for (uint_fast32_t i = 1; i < FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ;
         ++i) {
      _cumulative_distribution[i] += _cumulative_distribution[i - 1];
    }
    // the total ionizing luminosity is the last element of
    // _cumulative_distribution (using a zeroth order quadrature)
    // its value is in 10^-21 erg Hz^-1 s^-1 cm^-2 sr^-1
    // we convert to s^-1 cm^-2 sr^-1 by dividing by Planck's constant
    // (in J s = J Hz^-1) and multiplying by 10^-28
    _total_flux =
        1.e-28 *
        _cumulative_distribution[FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ -
                                 1] /
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
    // we integrate out over all solid angles
    _total_flux *= 4. * M_PI;
    // and convert from cm^-2 to m^-2
    _total_flux *= 1.e4;
    // normalize the spectrum
    for (uint_fast32_t i = 0; i < FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ;
         ++i) {
      _cumulative_distribution[i] /=
          _cumulative_distribution[FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ -
                                   1];
    }
  } else {
    // no UVB. We set the cumulative distribution and the total luminosity to 0
    for (uint_fast32_t i = 0; i < FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ;
         ++i) {
      _cumulative_distribution[i] = 0.;
    }
    _total_flux = 0.;
  }

  if (log) {
    log->write_status(
        "Constructed FaucherGiguerePhotonSourceSpectrum at redshift ", redshift,
        ", with total ionizing flux ", _total_flux, " m^-2 s^-1.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - redshift: Redshift for which we want to use the UVB spectrum (default: 0.)
 *
 * @param role Role the spectrum will fulfil in the simulation. Parameters are
 * read from the corresponding block in the parameter file.
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
FaucherGiguerePhotonSourceSpectrum::FaucherGiguerePhotonSourceSpectrum(
    std::string role, ParameterFile &params, Log *log)
    : FaucherGiguerePhotonSourceSpectrum(
          params.get_value< double >(role + ":redshift", 0.), log) {}

/**
 * @brief Get the name of the file containing the spectrum for the given
 * redshift.
 *
 * @param z Redshift.
 * @return Name of the file containing the spectrum at that redshift.
 */
std::string FaucherGiguerePhotonSourceSpectrum::get_filename(double z) {
  // due to the strange way Faucher-Giguere's file names are formatted, and
  // since we need to be sure round off does not affect the names we generate,
  // we find the name using integer arithmetics rather than floating points
  uint_fast32_t iz = std::round(z / 0.05) * 5;
  const uint_fast32_t iz100 = iz / 100;
  iz -= iz100 * 100;
  const uint_fast32_t iz10 = iz / 10;
  iz -= iz10 * 10;
  std::stringstream namestream;
  namestream << FAUCHERGIGUEREDATALOCATION << "fg_uvb_dec11_z_" << iz100 << "."
             << iz10;
  if (iz > 0) {
    namestream << iz;
  }
  namestream << ".dat";
  std::string name = namestream.str();
  // check if the file exists
  std::ifstream file(name);
  if (!file) {
    cmac_error("File not found: %s!", name.c_str());
  }
  return name;
}

/**
 * @brief Get the total ionizing flux.
 *
 * @return Total ionizing flux (in m^-2 s^-1).
 */
double FaucherGiguerePhotonSourceSpectrum::get_total_flux() const {
  return _total_flux;
}

/**
 * @brief Get a random frequency from the spectrum.
 *
 * @param random_generator RandomGenerator to use.
 * @param temperature Temperature (in K). Not used for this spectrum.
 * @return Random frequency (in Hz).
 */
double FaucherGiguerePhotonSourceSpectrum::get_random_frequency(
    RandomGenerator &random_generator, double temperature) const {

  const double x = random_generator.get_uniform_random_double();
  const uint_fast32_t inu =
      Utilities::locate(x, _cumulative_distribution.data(),
                        FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ);
  const double frequency =
      _frequencies[inu] +
      (_frequencies[inu + 1] - _frequencies[inu]) *
          (x - _cumulative_distribution[inu]) /
          (_cumulative_distribution[inu + 1] - _cumulative_distribution[inu]);
  return frequency;
}
