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
 * @file PhotonSourceSpectrumMask.hpp
 *
 * @brief General interface for PhotonSourceSpectrum masks.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCESPECTRUMMASK_HPP
#define PHOTONSOURCESPECTRUMMASK_HPP

/**
 * @brief General interface for PhotonSourceSpectrum masks.
 */
class PhotonSourceSpectrumMask {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~PhotonSourceSpectrumMask() {}

  /**
   * @brief Get the fraction of the spectrum at the given frequency that should
   * be retained in the masked spectrum.
   *
   * @param frequency Frequency of the bin (in Hz).
   * @return Fraction of the spectrum that should be retained.
   */
  virtual double get_bin_fraction(double frequency) const = 0;
};

#endif // PHOTONSOURCESPECTRUMMASK_HPP
