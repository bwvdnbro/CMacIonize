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
 * @file WMBasicPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for the Pauldrach, Hoffmann &
 * Lennon (2001) stellar model spectra.
 *
 * The data used comes from Sternberg, Hoffmann & Pauldrach (2003), and was
 * downloaded from the link in that paper (ftp://wise3.tau.ac.il/pub/stars).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef WMBASICPHOTONSOURCESPECTRUM_HPP
#define WMBASICPHOTONSOURCESPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;
class RandomGenerator;

/**
 * @brief Number of frequency bins used in the internal table.
 */
#define WMBASICPHOTONSOURCESPECTRUM_NUMFREQ 1000

/**
 * @brief PhotonSourceSpectrum implementation for the Pauldrach, Hoffmann &
 * Lennon (2001) stellar model spectra.
 */
class WMBasicPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  std::vector< double > _frequencies;

  /*! @brief Cumulative distribution of the spectrum. */
  std::vector< double > _cumulative_distribution;

  /*! @brief Total ionizing flux of the spectrum (in m^-2 s^-1). */
  double _total_flux;

public:
  WMBasicPhotonSourceSpectrum(double temperature, double surface_gravity,
                              Log *log = nullptr);

  WMBasicPhotonSourceSpectrum(std::string role, ParameterFile &params,
                              Log *log = nullptr);

  static std::string get_log_g_name(double surface_gravity);
  static std::string get_filename(double temperature, double surface_gravity);

  /**
   * @brief Virtual destructor.
   */
  virtual ~WMBasicPhotonSourceSpectrum() {}

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const;

  virtual double get_total_flux() const;
};

#endif // WMBASICPHOTONSOURCESPECTRUM_HPP
