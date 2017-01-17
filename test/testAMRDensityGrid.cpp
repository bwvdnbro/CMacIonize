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
   * @brief Get the density at the given position.
   *
   * @param position CoordinateVector specifying a position.
   * @return Density at that position.
   */
  DensityValues operator()(CoordinateVector<> position) const {
    DensityValues cell;

    double density;
    if (position.z() < 0.5) {
      density = 1.;
    } else {
      density = 2.;
    }

    cell.set_total_density(density);
    cell.set_temperature(4000.);
    return cell;
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
   * @param midpoint Midpoint of the cell (in m).
   * @param volume Volume of the cell (in m^3).
   * @param cell DensityValues of a cell.
   * @return True if the density is larger than 1.
   */
  virtual bool refine(unsigned char level, CoordinateVector<> midpoint,
                      double volume, DensityValues &cell) const {
    return cell.get_total_density() > 1. && level < 6;
  }

  /**
   * @brief Decide if the given cells should be replaced by a single cell or
   * not.
   *
   * @param level Current refinement level of the cells.
   * @param midpoints Midpoints of the cells (in m).
   * @param volumes  Volumes of the cells (in m^3).
   * @param cells DensityValues of the cells.
   * @return True if the cells can be replaced by a single cell on a coarser
   * level.
   */
  virtual bool coarsen(unsigned char level, CoordinateVector<> *midpoints,
                       double *volumes, DensityValues *cells) const {
    double avg_density = 0.;
    for (unsigned int i = 0; i < 8; ++i) {
      avg_density += cells->get_total_density();
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
  AMRRefinementScheme *scheme = new TestAMRRefinementScheme();
  TerminalLog log(LOGLEVEL_INFO);
  AMRDensityGrid grid(Box(CoordinateVector<>(0.), CoordinateVector<>(1.)), 32,
                      density_function, scheme, false, &log);
  grid.initialize();

  assert_values_equal(1.5, grid.get_total_hydrogen_number());
  assert_values_equal(4000., grid.get_average_temperature());

  assert_condition(grid.get_number_of_cells() == 32 * 32 * 16 + 32 * 64 * 64);
  // 32768 is the grid in the lower left front corner, which has indices 000 on
  // all levels
  // the 32768 bit is set to indicate its level: 5
  unsigned long key = grid.get_cell_index(CoordinateVector<>(0.01));
  assert_condition(key == 32768);

  // pick the first cell to check the volume and midpoint calculation
  // due to the order in which the cell list is constructed, this should be the
  // cell with key 32768
  assert_condition(grid.get_cell_volume(0) ==
                   (1. / 32) * (1. / 32) * (1. / 32));
  CoordinateVector<> midpoint = grid.get_cell_midpoint(0);
  assert_condition(midpoint.x() == 0.015625);
  assert_condition(midpoint.y() == 0.015625);
  assert_condition(midpoint.z() == 0.015625);

  unsigned int ncell = 0;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    ++ncell;
  }
  assert_condition(ncell == grid.get_number_of_cells());

  return 0;
}
