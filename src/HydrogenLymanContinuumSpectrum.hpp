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

class HydrogenLymanContinuumSpectrum : public PhotonSourceSpectrum {
public:
  virtual double get_random_frequency();
};

#endif // HYDROGENLYMANCONTINUUMSPECTRUM_HPP
