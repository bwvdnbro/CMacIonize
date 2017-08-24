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
 * @file testSPHArrayDensityGridWriter.cpp
 *
 * @brief Unit test for the SPHArrayDensityGridWriter class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "Octree.hpp"
#include "SPHArrayDensityGridWriter.hpp"

class TestDensityFunction : public DensityFunction {
public:
  /**
   * @brief Get the neutral fraction for the given cell.
   *
   * @param cell Cell.
   * @return Neutral fraction.
   */
  virtual DensityValues operator()(const Cell &cell) const {
    const CoordinateVector<> p = cell.get_cell_midpoint();
    DensityValues values;
    if ((p - CoordinateVector<>(0.5)).norm() < 0.2) {
      values.set_ionic_fraction(ION_H_n, 1.);
    } else {
      values.set_ionic_fraction(ION_H_n, 0.);
    }
    return values;
  }
};

/**
 * @brief Unit test for the SPHArrayDensityGridWriter class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CartesianDensityGrid grid(box, 10);
  TestDensityFunction density_function;
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  std::vector< CoordinateVector<> > sph_positions(1000);
  for (unsigned int ix = 0; ix < 10; ++ix) {
    for (unsigned int iy = 0; iy < 10; ++iy) {
      for (unsigned int iz = 0; iz < 10; ++iz) {
        sph_positions[ix * 100 + iy * 10 + iz][0] = (ix + 0.5) * 0.1;
        sph_positions[ix * 100 + iy * 10 + iz][1] = (iy + 0.5) * 0.1;
        sph_positions[ix * 100 + iy * 10 + iz][2] = (iz + 0.5) * 0.1;
      }
    }
  }
  Octree octree(sph_positions, box);
  SPHArrayDensityGridWriter density_grid_writer;
  density_grid_writer.reset(1000, &octree);
  ParameterFile params;
  density_grid_writer.write(grid, 0, params, 0.);

  std::vector< double > neutral_fractions(1000, -1.);
  density_grid_writer.fill_array(neutral_fractions.data());

  for (unsigned int i = 0; i < 1000; ++i) {
    const CoordinateVector<> p = sph_positions[i];
    if ((p - CoordinateVector<>(0.5)).norm() < 0.2) {
      assert_condition(neutral_fractions[i] == 1.);
    } else {
      assert_condition(neutral_fractions[i] == 0.);
    }
  }

  return 0;
}
