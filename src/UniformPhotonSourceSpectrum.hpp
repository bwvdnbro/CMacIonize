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
 * @file UniformPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for a uniform spectrum.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UNIFORMPHOTONSOURCESPECTRUM_HPP
#define UNIFORMPHOTONSOURCESPECTRUM_HPP

#include "Error.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief PhotonSourceSpectrum implementation for a uniform spectrum.
 */
class UniformPhotonSourceSpectrum : public PhotonSourceSpectrum {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~UniformPhotonSourceSpectrum() {}

  /**
   * @brief Get a random frequency from the spectrum.
   *
   * @param random_generator RandomGenerator to use.
   * @param temperature Temperature of the gas (for reemission spectra) (in K).
   * @return Random frequency, distributed according to the spectrum (in Hz).
   */
  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const {
    return (1. + 3. * random_generator.get_uniform_random_double()) * 3.289e15;
  }

  /**
   * @brief Get the total ionizing flux emitted by the spectrum.
   *
   * @return Total ionizing flux (in m^-2 s^-1).
   */
  virtual double get_total_flux() const {
    cmac_error("This function should not be used!");
    return 0.;
  }
};

#endif // UNIFORMPHOTONSOURCESPECTRUM_HPP
