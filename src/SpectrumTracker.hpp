/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SpectrumTracker.hpp
 *
 * @brief Instrument that keeps track of the spectrum of photons that pass
 * through a cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPECTRUMTRACKER_HPP
#define SPECTRUMTRACKER_HPP

#include "Photon.hpp"

#include <cinttypes>
#include <fstream>
#include <string>
#include <vector>

/**
 * @brief Instrument that keeps track of the spectrum of photons that pass
 * through a cell.
 */
class SpectrumTracker {
private:
  /*! @brief Minimum frequency for the spectral bins (in Hz). */
  const double _minimum_frequency;

  /*! @brief Width of a single frequency bin (in Hz). */
  const double _frequency_width;

  /*! @brief Inverse width of a single frequency bin (in Hz^-1). */
  const double _inverse_frequency_width;

  /*! @brief Number counts per bin for direct radiation from a source. */
  std::vector< uint_fast64_t > _number_counts_primary;

  /*! @brief Number counts per bin for diffuse hydrogen reemission. */
  std::vector< uint_fast64_t > _number_counts_diffuse_H;

  /*! @brief Number counts per bin for diffuse helium reemission. */
  std::vector< uint_fast64_t > _number_counts_diffuse_He;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_bins Number of bins to use.
   */
  SpectrumTracker(const uint_fast32_t number_of_bins)
      : _minimum_frequency(3.289e15),
        _frequency_width(3. * 3.289e15 / number_of_bins),
        _inverse_frequency_width(1. / _frequency_width),
        _number_counts_primary(number_of_bins, 0),
        _number_counts_diffuse_H(number_of_bins, 0),
        _number_counts_diffuse_He(number_of_bins, 0) {}

  /**
   * @brief Add the contribution of the given photon packet to the bins.
   *
   * @param photon Photon to add.
   */
  inline void count_photon(const Photon &photon) {

    const double frequency = photon.get_energy();
    const uint_fast32_t index =
        (frequency - _minimum_frequency) * _inverse_frequency_width;
    if (index < _number_counts_primary.size()) {
      if (photon.get_type() == PHOTONTYPE_PRIMARY) {
        ++_number_counts_primary[index];
      }
      if (photon.get_type() == PHOTONTYPE_DIFFUSE_HI) {
        ++_number_counts_diffuse_H[index];
      }
      if (photon.get_type() == PHOTONTYPE_DIFFUSE_HeI) {
        ++_number_counts_diffuse_He[index];
      }
    }
  }

  /**
   * @brief Output the spectrum to the file with the given name.
   *
   * @param filename Name of the output file.
   */
  inline void output_spectrum(const std::string filename) const {

    std::ofstream ofile(filename);
    ofile << "# frequency (Hz)\tprimary count\tdiffuse H count\tdiffuse He "
             "count\n";
    for (uint_fast32_t i = 0; i < _number_counts_primary.size(); ++i) {
      const double nu = _minimum_frequency + (i + 0.5) * _frequency_width;
      ofile << nu << "\t" << _number_counts_primary[i] << "\t"
            << _number_counts_diffuse_H[i] << "\t"
            << _number_counts_diffuse_He[i] << "\n";
    }
  }
};

#endif // SPECTRUMTRACKER_HPP
