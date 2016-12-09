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
 * @file PhotonSourceSpectrum.hpp
 *
 * @brief Interface for a photon source spectrum.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCESPECTRUM_HPP
#define PHOTONSOURCESPECTRUM_HPP

/**
 * @brief Interface for a photon source spectrum.
 */
class PhotonSourceSpectrum {
public:
  /**
   * @brief Get a random frequency from the spectrum.
   *
   * @param temperature Temperature of the gas (for reemission spectra) (in K).
   * @return Random frequency, distributed according to the spectrum (in Hz).
   */
  virtual double get_random_frequency(double temperature = 0.) = 0;

  /**
   * @brief Get the total ionizing flux emitted by the spectrum.
   *
   * @return Total ionizing flux (in m^-2 s^-1).
   */
  virtual double get_total_flux() = 0;
};

#endif // PHOTONSOURCESPECTRUM_HPP
