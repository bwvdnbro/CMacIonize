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
#include "RandomGenerator.hpp"

class CrossSections;

/*! @brief Number of frequencies in the internal table. */
#define HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ 1000

/*! @brief Number of temperatures in the internal table. */
#define HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP 100

/**
 * @brief Helium Lyman continuum photoionization spectrum.
 */
class HeliumLymanContinuumSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins (in 13.6 eV). */
  double _frequency[HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ];

  /*! @brief Temperature bins (in K). */
  double _temperature[HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP];

  /*! @brief Cumulative distribution function. */
  double _cumulative_distribution[HELIUMLYMANCONTINUUMSPECTRUM_NUMTEMP]
                                 [HELIUMLYMANCONTINUUMSPECTRUM_NUMFREQ];

  /*! @brief Current temperature (in K). */
  double _current_T;

  /*! @brief Index of the current temperature bin. */
  unsigned int _current_iT;

  /*! @brief RandomGenerator used to generate random numbers. */
  RandomGenerator &_random_generator;

public:
  HeliumLymanContinuumSpectrum(CrossSections &cross_sections,
                               RandomGenerator &random_generator);

  void set_temperature(double T);

  virtual double get_random_frequency();

  virtual double get_total_flux();
};

#endif // HELIUMLYMANCONTINUUMSPECTRUM_HPP
