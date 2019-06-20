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
 * @file MaskedPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum that masks another PhotonSourceSpectrum.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MASKEDPHOTONSOURCESPECTRUM_HPP
#define MASKEDPHOTONSOURCESPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

#include <cstdint>
#include <string>
#include <vector>

class Log;
class ParameterFile;
class PhotonSourceSpectrumMask;

/**
 * @brief PhotonSourceSpectrum that masks another PhotonSourceSpectrum.
 */
class MaskedPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins (in Hz). */
  std::vector< double > _frequency_bins;

  /*! @brief Cumulative distribution in each bin. */
  std::vector< double > _cumulative_distribution;

  /*! @brief Total ionizing flux of the spectrum (in m^-2 s^-1). */
  double _ionizing_flux;

public:
  MaskedPhotonSourceSpectrum(PhotonSourceSpectrum *unmasked_spectrum,
                             PhotonSourceSpectrumMask *mask,
                             const uint_fast16_t number_of_bins,
                             const uint_fast32_t number_of_samples);

  MaskedPhotonSourceSpectrum(std::string role, ParameterFile &params,
                             Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~MaskedPhotonSourceSpectrum() {}

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature) const;

  virtual double get_total_flux() const;

  // unit testing routines

  std::pair< std::vector< double >, std::vector< double > >
  get_spectrum() const;
};

#endif // MASKEDPHOTONSOURCESPECTRUM_HPP
