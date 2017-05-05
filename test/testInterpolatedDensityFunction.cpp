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
 * @file testInterpolatedDensityFunction.cpp
 *
 * @brief Unit test for the InterpolatedDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "InterpolatedDensityFunction.hpp"

/**
 * @brief Get the number density at the given vertical position.
 *
 * The test number density is @f$n = (z-1)^2 @f$.
 *
 * @param z Vertical position (in m).
 * @return Number density (in m^-3).
 */
double get_number_density(double z) { return 1. + z * z; }

/**
 * @brief Unit test for the InterpolateDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  InterpolatedDensityFunction density_function("test_interpolated_density.txt",
                                               4000.);

  for (unsigned int i = 0; i < 1000; ++i) {
    CoordinateVector<> p(0.5, 0.5, (i + 0.5) * 0.001);
    assert_values_equal_rel(density_function(p).get_number_density(),
                            get_number_density(p.z()), 1.e-4);
    assert_condition(density_function(p).get_temperature() == 4000.);
  }

  return 0;
}
