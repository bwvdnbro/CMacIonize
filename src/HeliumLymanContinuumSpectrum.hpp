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
 * @file HeliumLymanContinuumSpectrum.hpp
 *
 * @brief Helium Lyman continuum PhotonSourceSpectrum implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HELIUMLYMANCONTINUUMSPECTRUM_HPP
#define HELIUMLYMANCONTINUUMSPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

class CrossSections;

#define HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ 1000
#define HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP 100

class HeliumLymanContinuumSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  double _frequency[HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ];

  /*! @brief Temperature bins. */
  double _temperature[HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP];

  /*! @brief Cumulative distribution function. */
  double _cumulative_distribution[HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP]
                                 [HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ];

  /*! @brief Current temperature. */
  double _current_T;

public:
  HeliumLymanContinuumSpectrum(CrossSections &cross_sections);

  void set_temperature(double T);

  virtual double get_random_frequency();
};

#endif // HELIUMLYMANCONTINUUMSPECTRUM_HPP
