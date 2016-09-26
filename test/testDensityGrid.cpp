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
 * @file testDensityGrid.cpp
 *
 * @brief Unit test for the DensityGrid class
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "DensityFunction.hpp"
#include "DensityGrid.hpp"

/**
 * @brief Test implementation of DensityFunction.
 */
class TestDensityFunction : public DensityFunction {
  /**
   * @brief Get the density at the given coordinate.
   *
   * @param position CoordinateVector specifying a coordinate position.
   * @return A constant density 1.
   */
  virtual double operator()(CoordinateVector position) { return 1.; }
};

/**
 * @brief Unit test for the DensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit value: 0 on success.
 */
int main(int argc, char **argv) {
  TestDensityFunction testfunction;
  DensityGrid grid(0., 0., 0., 1., 1., 1., 64, testfunction);

  assert_values_equal(1., grid.get_total_mass());

  return 0;
}
