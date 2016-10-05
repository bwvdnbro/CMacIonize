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

/**
 * @brief Number of frequency bins used in the internal table.
 */
#define PLANCKPHOTONSOURCESPECTRUM_NUMFREQ 1000

/**
 * @brief PhotonSourceSpectrum implementation for a Planck blackbody spectrum.
 */
class PlanckPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  double _frequency[PLANCKPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Base 10 logarithm of the frequency bins. */
  double _log_frequency[PLANCKPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Luminosity in each bin. */
  double _luminosity[PLANCKPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Cumulative distribution in each bin. */
  double _cumulative_distribution[PLANCKPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Base 10 logarithm of the cumulative distribution in each bin. */
  double _log_cumulative_distribution[PLANCKPHOTONSOURCESPECTRUM_NUMFREQ];

public:
  PlanckPhotonSourceSpectrum();

  virtual double get_random_frequency();
};

#endif // PLANCKPHOTONSOURCESPECTRUM_HPP
