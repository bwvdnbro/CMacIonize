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
 * @file testVoronoiDensityGrid.cpp
 *
 * @brief Unit test for the VoronoiDensityGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "UniformRandomVoronoiGeneratorDistribution.hpp"
#include "UniformRegularVoronoiGeneratorDistribution.hpp"
#include "VoronoiDensityGrid.hpp"
#include "VoronoiGeneratorDistribution.hpp"
/**
 * @brief Unit test for the VoronoiDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  /// random uniform generators
  {
    HomogeneousDensityFunction density_function(1., 2000.);
    density_function.initialize();
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    UniformRandomVoronoiGeneratorDistribution *test_positions =
        new UniformRandomVoronoiGeneratorDistribution(box, 100, 42);
    VoronoiDensityGrid grid(test_positions, box, "Old", 0, false, false, 0.,
                            5. / 3., nullptr);
    std::pair< unsigned long, unsigned long > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, density_function);

    assert_values_equal(1., grid.get_total_hydrogen_number());
    assert_values_equal(2000., grid.get_average_temperature());
  }

  /// regular generators
  {
    HomogeneousDensityFunction density_function(1., 2000.);
    density_function.initialize();
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    UniformRegularVoronoiGeneratorDistribution *test_positions =
        new UniformRegularVoronoiGeneratorDistribution(
            box, CoordinateVector< unsigned int >(5));
    VoronoiDensityGrid grid(test_positions, box, "Old", 0, false, false, 0.,
                            5. / 3., nullptr);
    std::pair< unsigned long, unsigned long > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, density_function);

    assert_values_equal(1., grid.get_total_hydrogen_number());
    assert_values_equal(2000., grid.get_average_temperature());
  }

  return 0;
}
