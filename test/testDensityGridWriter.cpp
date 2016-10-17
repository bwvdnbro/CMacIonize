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
 * @file testDensityGridWriter.cpp
 *
 * @brief Unit test for the DensityGridWriter class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "DensityGrid.hpp"
#include "DensityGridWriter.hpp"
#include "RecombinationRates.hpp"

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
  virtual double operator()(CoordinateVector<> position) { return 1.; }
};

/**
 * @brief Test implementation of RecombinationRates.
 */
class TestRecombinationRates : public RecombinationRates {
public:
  /**
   * @brief Get the recombination rate for the given element at the given
   * temperature.
   *
   * @param element ElementName for an element.
   * @param temperature Temperature.
   * @return Recombination rate.
   */
  virtual double get_recombination_rate(ElementName element,
                                        double temperature) {
    return 1.;
  }
};

/**
 * @brief Unit test for the DensityGridWriter class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  CoordinateVector<> origin;
  CoordinateVector<> side(1.);
  Box box(origin, side);
  CoordinateVector< unsigned char > ncell(8);
  TestDensityFunction density_function;
  TestRecombinationRates recombination_rates;
  DensityGrid grid(box, ncell, 0.1, 8000., density_function,
                   recombination_rates);

  DensityGridWriter writer("testgrid.hdf5", grid);
  writer.write();

  return 0;
}
