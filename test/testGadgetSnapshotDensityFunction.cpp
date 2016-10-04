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
#include "CrossSections.hpp"
#include "DensityGrid.hpp"
#include "Error.hpp"
#include "GadgetSnapshotDensityFunction.hpp"
using namespace std;

/**
 * @brief Test implementation of CrossSections.
 */
class TestCrossSections : public CrossSections {
public:
  /**
   * @brief Get the photoionization cross section for the given element at the
   * given photon energy.
   *
   * @param element CrossSectionElements index for an element.
   * @param energy Photon energy.
   * @return Photoionization cross section.
   */
  virtual double get_cross_section(int element, double energy) { return 1.; }
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
  GadgetSnapshotDensityFunction density("test.hdf5");
  TestCrossSections testcrosssections;

  CoordinateVector<> anchor;
  CoordinateVector<> sides(1., 1., 1.);
  Box box(anchor, sides);
  DensityGrid grid(box, 32, density, testcrosssections);
  assert_values_equal(grid.get_total_mass(), density.get_total_mass());

  return 0;
}
