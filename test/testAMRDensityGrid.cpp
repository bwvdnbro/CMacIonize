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
 * @file testAMRDensityGrid.cpp
 *
 * @brief Unit test for the AMRDensityGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AMRDensityGrid.hpp"
#include "Assert.hpp"
#include "TerminalLog.hpp"

/**
 * @brief Test implementation of DensityFunction.
 */
class TestDensityFunction : public DensityFunction {
public:
  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  DensityValues operator()(const Cell &cell) const {
    DensityValues values;

    const CoordinateVector<> position = cell.get_cell_midpoint();

    double density;
    if (position.z() < 0.5) {
      density = 1.;
    } else {
      density = 2.;
    }

    values.set_number_density(density);
    values.set_temperature(4000.);
    return values;
  }
};

/**
 * @brief Test implementation of AMRRefinementScheme.
 */
class TestAMRRefinementScheme : public AMRRefinementScheme {
public:
  /**
   * @brief Decide if the given cell should be refine or not.
   *
   * @param level Current refinement level of the cell.
   * @param cell DensityGrid::iterator pointing to a cell.
   * @return True if the density is larger than 1.
   */
  virtual bool refine(uint_fast8_t level, DensityGrid::iterator &cell) const {
    return cell.get_ionization_variables().get_number_density() > 1. &&
           level < 6;
  }

  /**
   * @brief Decide if the given cells should be replaced by a single cell or
   * not.
   *
   * @param level Current refinement level of the cells.
   * @param cells DensityGrid::iterators pointing to the cells.
   * @return True if the cells can be replaced by a single cell on a coarser
   * level.
   */
  virtual bool coarsen(uint_fast8_t level,
                       const DensityGrid::iterator *cells) const {
    double avg_density = 0.;
    for (uint_fast8_t i = 0; i < 8; ++i) {
      avg_density += cells[i].get_ionization_variables().get_number_density();
    }
    avg_density *= 0.125;
    return avg_density < 1.;
  }
};

/**
 * @brief Unit test for the AMRDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  TestDensityFunction density_function;
  density_function.initialize();
  AMRRefinementScheme *scheme = new TestAMRRefinementScheme();
  TerminalLog log(LOGLEVEL_INFO);
  AMRDensityGrid grid(Box<>(CoordinateVector<>(0.), CoordinateVector<>(1.)), 32,
                      scheme, 5, false, false, &log);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  assert_values_equal(1.5, grid.get_total_hydrogen_number());
  assert_values_equal(4000., grid.get_average_temperature());

  assert_condition(grid.get_number_of_cells() == 32 * 32 * 16 + 32 * 64 * 64);
  // 32768 is the grid in the lower left front corner, which has indices 000 on
  // all levels
  // the 32768 bit is set to indicate its level: 5
  amrkey_t key = grid.get_cell_index(CoordinateVector<>(0.01));
  assert_condition(key == 0);

  // pick the first cell to check the volume and midpoint calculation
  // due to the order in which the cell list is constructed, this should be the
  // cell with key 32768
  assert_condition(grid.get_cell_volume(0) ==
                   (1. / 32) * (1. / 32) * (1. / 32));
  CoordinateVector<> midpoint = grid.get_cell_midpoint(0);
  assert_condition(midpoint.x() == 0.015625);
  assert_condition(midpoint.y() == 0.015625);
  assert_condition(midpoint.z() == 0.015625);

  uint_fast32_t ncell = 0;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    ++ncell;
  }
  assert_condition(ncell == grid.get_number_of_cells());

  return 0;
}
