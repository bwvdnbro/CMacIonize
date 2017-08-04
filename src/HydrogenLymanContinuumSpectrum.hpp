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
 * @file HydrogenLymanContinuumSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for the hydrogen Lyman continuum
 * spectrum.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROGENLYMANCONTINUUMSPECTRUM_HPP
#define HYDROGENLYMANCONTINUUMSPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"

#include <vector>

class CrossSections;

/*! @brief Number of frequencies in the internal table. */
#define HYDROGENLYMANCONTINUUMSPECTRUM_NUMFREQ 1000

/*! @brief Number of temperatures in the internal table. */
#define HYDROGENLYMANCONTINUUMSPECTRUM_NUMTEMP 100

/**
 * @brief Hydrogen Lyman continuum photoionization spectrum.
 *
 * We use the spectrum given by Wood, K., Mathis, J. S. & Ercolano, B. 2004,
 * MNRAS, 348, 1337 (http://adsabs.harvard.edu/abs/2004MNRAS.348.1337W),
 * equation (8), which uses the same ionization cross sections that are used in
 * other parts of the program. We pretabulate values in a 2D temperature
 * frequency space in the range [1,500 K; 15,000 K] (for temperature values
 * outside this range, extrapolation is safe).
 */
class HydrogenLymanContinuumSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins (in 13.6 eV). */
  std::vector< double > _frequency;

  /*! @brief Temperature bins (in K). */
  std::vector< double > _temperature;

  /*! @brief Cumulative distribution function. */
  std::vector< std::vector< double > > _cumulative_distribution;

public:
  HydrogenLymanContinuumSpectrum(CrossSections &cross_sections);

  /**
   * @brief Virtual destructor.
   */
  virtual ~HydrogenLymanContinuumSpectrum() {}

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature) const;

  virtual double get_total_flux() const;
};

#endif // HYDROGENLYMANCONTINUUMSPECTRUM_HPP
