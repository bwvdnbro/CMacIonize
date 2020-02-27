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
 * @file PopStarPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for PopStar stellar spectra.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef POPSTARPHOTONSOURCESPECTRUM_HPP
#define POPSTARPHOTONSOURCESPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;
class RandomGenerator;

/**
 * @brief Number of frequency bins used in the internal table.
 */
#define POPSTARPHOTONSOURCESPECTRUM_NUMFREQ 1000

/**
 * @brief PhotonSourceSpectrum implementation for PopStar stellar spectra.
 *
 * These spectra are taken from https://www.fractal-es.com/PopStar/#download
 *
 * We resample the spectra on a frequency grid of 1000 bins in the range
 * [13.6 eV, 54.4 eV].
 */
class PopStarPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  std::vector< double > _frequencies;

  /*! @brief Cumulative distribution of the spectrum. */
  std::vector< double > _cumulative_distribution;

  /*! @brief Total ionizing flux of the spectrum (in m^-2 s^-1). */
  double _total_flux;

public:
  static std::string get_filename(const double log_age_in_yr,
                                  const double metallicity);

  PopStarPhotonSourceSpectrum(const double log_age_in_yr,
                              const double metallicity, Log *log = nullptr);

  PopStarPhotonSourceSpectrum(std::string role, ParameterFile &params,
                              Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~PopStarPhotonSourceSpectrum() {}

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const;

  virtual double get_total_flux() const;
};

#endif // POPSTARPHOTONSOURCESPECTRUM_HPP
