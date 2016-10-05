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
 * @file VernerRecombinationRates.cpp
 *
 * @brief RecombinationRates implementation with Verner's recombination rates:
 * implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "VernerRecombinationRates.hpp"
#include "Error.hpp"
#include <cmath>
using namespace std;

/**
 * @brief Get the recombination rate of the given element at the given
 * temperature.
 *
 * @param element ElementName for an element.
 * @param temperature Temperature.
 * @return Recombination rate.
 */
double VernerRecombinationRates::get_recombination_rate(ElementName element,
                                                        double temperature) {
  switch (element) {
  case ELEMENT_H:
    return 7.982e-11 / (sqrt(temperature / 3.148) *
                        pow(1. + sqrt(temperature / 3.148), 0.252) *
                        pow(1. + sqrt(temperature / 7.036e5), 1.748));
  case ELEMENT_He:
    return 3.294e-11 / (sqrt(temperature / 15.54) *
                        pow(1. + sqrt(temperature / 15.54), 0.309) *
                        pow(1. + sqrt(temperature / 3.676e7), 1.691));
  default:
    error("Unknown element: %i", element);
  }
  return 0.;
}
