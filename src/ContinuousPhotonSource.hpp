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
 * @file ContinuousPhotonSource.hpp
 *
 * @brief General interface for external photon sources.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CONTINUOUSPHOTONSOURCE_HPP
#define CONTINUOUSPHOTONSOURCE_HPP

#include "CoordinateVector.hpp"
#include "Error.hpp"

class RandomGenerator;

/**
 * @brief General interface for external photon sources.
 */
class ContinuousPhotonSource {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~ContinuousPhotonSource() {}

  /**
   * @brief Get the entrance position and direction of a random external photon.
   *
   * @param random_generator RandomGenerator to use.
   * @return std::pair of CoordinateVector instances, specifying a starting
   * position (in m) and direction for an incoming photon.
   */
  virtual std::pair< CoordinateVector<>, CoordinateVector<> >
  get_random_incoming_direction(RandomGenerator &random_generator) const = 0;

  /**
   * @brief Get the total surface area through which the radiation enters the
   * simulation box.
   *
   * This is used to calculate the luminosity of the incoming radiation.
   *
   * @return Total surface area through which radiation enters the box (in m^2).
   */
  virtual double get_total_surface_area() const = 0;

  /**
   * @brief Does this ContinuousPhotonSource have a total luminosity value?
   *
   * @return False, as the total luminosity by default depends on the spectrum.
   */
  virtual bool has_total_luminosity() const { return false; }

  /**
   * @brief Get the total luminosity for the source.
   *
   * Only provided if has_total_luminosity() returns true.
   *
   * @return Total ionizing luminosity of the source (in s^-1).
   */
  virtual double get_total_luminosity() const {
    cmac_error(
        "This function should always be implemented for "
        "ContinuousPhotonSources that provide a total ionizing luminosity!");
    return 0;
  }
};

#endif // CONTINUOUSPHOTONSOURCE_HPP
