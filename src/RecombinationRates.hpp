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
 * @file RecombinationRates.hpp
 *
 * @brief General interface for recombination rates.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RECOMBINATIONRATES_HPP
#define RECOMBINATIONRATES_HPP

#include "ElementNames.hpp"

/**
 * @brief General interface for recombination rates.
 */
class RecombinationRates {
public:
  /**
   * @brief Get the recombination rate for the given ion at the given
   * temperature.
   *
   * @param ion IonName for a valid ion.
   * @param temperature Temperature (in K).
   * @return Recombination rate (in m^3s^-1).
   */
  virtual double get_recombination_rate(IonName ion,
                                        double temperature) const = 0;
};

#endif // RECOMBINATIONRATES_HPP
