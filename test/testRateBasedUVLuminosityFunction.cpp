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
 * @file testRateBasedUVLuminosityFunction.cpp
 *
 * @brief Unit test for the RateBasedUVLuminosityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "RateBasedUVLuminosityFunction.hpp"

/**
 * @brief Unit test for the RateBasedUVLuminosityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  RateBasedUVLuminosityFunction luminosity_function(2., 2.);

  assert_condition(luminosity_function(1., 1.) == 2.);

  return 0;
}
