/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file CaproniStellarFeedback.hpp
 *
 * @brief Dwarf galaxy stellar feedback prescription based on the SN rates in
 * Caproni et al. (2017).
 *
 * For more details, see Vandenbroucke et al., in prep.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CAPRONISTELLARFEEDBACK_HPP
#define CAPRONISTELLARFEEDBACK_HPP

#include "DensityGrid.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief Dwarf galaxy stellar feedback prescription based on the SN rates in
 * Caproni et al. (2017).
 */
class CaproniStellarFeedback {
public:
  /**
   * @brief Get the instanteneous supernova rate at the given time.
   *
   * We use a 8th order polynomial fit to the Caproni et al. (2017) data.
   *
   * @param time Current simulation time (in s).
   * @return Instanteneous supernova rate (in s^-1).
   */
  inline static double get_SN_rate(const double time) {

    if (time > 6.4e16) {
      cmac_error("Time value outside fit validity range!");
    }

    // fit coefficients
    const double a[9] = {
        6.89799700195e-143, -1.91715448814e-125, 2.19302246787e-108,
        -1.32201654133e-91, 4.43164838505e-75,   -7.80981517111e-59,
        5.39845077482e-43,  9.70796661139e-28,   -8.44606535214e-14};

    double rate = a[0] * time + a[1];
    rate = rate * time + a[2];
    rate = rate * time + a[3];
    rate = rate * time + a[4];
    rate = rate * time + a[5];
    rate = rate * time + a[6];
    rate = rate * time + a[7];
    rate = rate * time + a[8];

    return rate;
  }

  /**
   * @brief Constructor.
   *
   * @param log Log to write logging info to.
   */
  inline CaproniStellarFeedback(Log *log = nullptr) {}

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline CaproniStellarFeedback(ParameterFile &params, Log *log = nullptr)
      : CaproniStellarFeedback(log) {}

  /**
   * @brief Add stellar feedback at the given time.
   *
   * @param grid DensityGrid to operate on.
   * @param time Current simulation time (in s).
   */
  inline void add_stellar_feedback(DensityGrid &grid, const double time) {}
};

#endif // CAPRONISTELLARFEEDBACK_HPP
