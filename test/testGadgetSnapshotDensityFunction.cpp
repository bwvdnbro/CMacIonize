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
 * @file testGadgetSnapshotDensityFunction.cpp
 *
 * @brief Unit test for the GadgetSnapshotDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CoordinateVector.hpp"
#include "DensityGrid.hpp"
#include "Error.hpp"
#include "GadgetSnapshotDensityFunction.hpp"
#include "RecombinationRates.hpp"
#include "TerminalLog.hpp"
using namespace std;

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
 * @brief Unit test for the GadgetSnapshotDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // before we can test this, we need to make sure we can open and read a
  // Gadget2 snapshot file.
  TerminalLog tlog(LOGLEVEL_INFO);
  GadgetSnapshotDensityFunction density("test.hdf5", false, 0., 0., &tlog);
  TestRecombinationRates testrecombinationrates;

  CoordinateVector<> anchor;
  CoordinateVector<> sides(1., 1., 1.);
  Box box(anchor, sides);
  DensityGrid grid(box, 32, 0.1, 8000., density, testrecombinationrates);
  assert_values_equal(grid.get_total_hydrogen_number(),
                      density.get_total_hydrogen_number());

  return 0;
}
