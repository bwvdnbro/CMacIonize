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
 * @file PlanckPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for a Planck blackbody spectrum:
 * header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PLANCKPHOTONSOURCESPECTRUM_HPP
#define PLANCKPHOTONSOURCESPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;
class RandomGenerator;

/**
 * @brief Number of frequency bins used in the internal table.
 */
#define PLANCKPHOTONSOURCESPECTRUM_NUMFREQ 1000

/**
 * @brief PhotonSourceSpectrum implementation for a Planck blackbody spectrum.
 *
 * We used the spectrum found on the Wikipedia page about black body radiation:
 * https://en.wikipedia.org/wiki/Black-body_radiation.
 */
class PlanckPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Base 10 logarithm of the frequency bins (in log(frequency/13.6 eV).
   */
  std::vector< double > _log_frequency;

  /*! @brief Cumulative distribution in each bin. */
  std::vector< double > _cumulative_distribution;

  /*! @brief Base 10 logarithm of the cumulative distribution in each bin. */
  std::vector< double > _log_cumulative_distribution;

  /*! @brief Ionizing flux of the spectrum (in m^-2 s^-1). */
  const double _ionizing_flux;

public:
  PlanckPhotonSourceSpectrum(double temperature, double ionizing_flux = -1.,
                             Log *log = nullptr);

  PlanckPhotonSourceSpectrum(std::string role, ParameterFile &params,
                             Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~PlanckPhotonSourceSpectrum() {}

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const;

  virtual double get_total_flux() const;
};

#endif // PLANCKPHOTONSOURCESPECTRUM_HPP
