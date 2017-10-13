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
 * @file LinearPhotonSourceSpectrumMask.hpp
 *
 * @brief Linear PhotonSourceSpecrumMask that linearly masks out the high
 * frequency part of the spectrum.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LINEARPHOTONSOURCESPECTRUMMASK_HPP
#define LINEARPHOTONSOURCESPECTRUMMASK_HPP

#include "PhotonSourceSpectrumMask.hpp"

/**
 * @brief Linear PhotonSourceSpecrumMask that linearly masks out the high
 * frequency part of the spectrum.
 */
class LinearPhotonSourceSpectrumMask : public PhotonSourceSpectrumMask {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~LinearPhotonSourceSpectrumMask() {}

  /**
   * @brief Get the fraction of the spectrum at the given frequency that should
   * be retained in the masked spectrum.
   *
   * @param frequency Frequency of the bin (in Hz).
   * @return Fraction of the spectrum that should be retained.
   */
  virtual double get_bin_fraction(double frequency) const {

    const double min_frequency = 3.289e15;
    const double max_frequency = 4. * min_frequency;
    const double frequency_range = max_frequency - min_frequency;
    return 1. - (frequency - min_frequency) / frequency_range;
  }
};

#endif // LINEARPHOTONSOURCESPECTRUMMASK_HPP
