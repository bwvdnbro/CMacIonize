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
 * @file CastelliKuruczPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for the Castelli & Kurucz (2003)
 * stellar atmosphere models.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef CASTELLIKURUCZPHOTONSOURCESPECTRUM_HPP
#define CASTELLIKURUCZPHOTONSOURCESPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;
class RandomGenerator;
namespace HDF5Tools {
template < typename _datatype_, uint_fast8_t _size_ > class HDF5DataBlock;
}

/**
 * @brief Number of frequency bins used in the internal table.
 */
#define CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ 1000

/**
 * @brief PhotonSourceSpectrum implementation for the Castelli & Kurucz (2003)
 * (https://ui.adsabs.harvard.edu/abs/2003IAUS..210P.A20C) stellar
 * atmosphere models.
 *
 * These spectra are stored in an HDF5 file that was constructed from the
 * SKIRT stored table data file distributed along with SKIRT
 * (https://skirt.ugent.be), which in itself was based on the online tables
 * provided by Castelli & Kurucz, which are no longer available from the
 * original URL.
 */
class CastelliKuruczPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  std::vector< double > _frequencies;

  /*! @brief Cumulative distribution of the spectrum. */
  std::vector< double > _cumulative_distribution;

  /*! @brief Total ionizing flux of the spectrum (in m^-2 s^-1). */
  double _total_flux;

  static double interpolate(const double fZ, const double fT, const double fg,
                            const double fl, const uint_fast32_t iZ,
                            const uint_fast32_t iT, const uint_fast32_t ig,
                            const uint_fast32_t il,
                            const HDF5Tools::HDF5DataBlock< double, 4 > &F);

public:
  CastelliKuruczPhotonSourceSpectrum(const double temperature,
                                     const double surface_gravity,
                                     const double metallicity,
                                     Log *log = nullptr);

  CastelliKuruczPhotonSourceSpectrum(std::string role, ParameterFile &params,
                                     Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~CastelliKuruczPhotonSourceSpectrum() {}

  static bool is_valid(const double temperature, const double surface_gravity,
                       const double metallicity);

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const;

  virtual double get_total_flux() const;
};

#endif // CASTELLIKURUCZPHOTONSOURCESPECTRUM_HPP
