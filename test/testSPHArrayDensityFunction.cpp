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
 * @file testSPHArrayDensityFunction.cpp
 *
 * @brief Unit test for the SPHArrayDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AsciiFileDensityGridWriter.hpp"
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "SPHArrayDensityFunction.hpp"
#include "Utilities.hpp"

/**
 * @brief Unit test for the SPHArrayDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CartesianDensityGrid grid(box, 32);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  ParameterFile params;
  SPHArrayDensityFunction density_function(1., 1.);

  /// double precision arrays
  {
    std::vector< double > x(1000, 0.);
    std::vector< double > y(1000, 0.);
    std::vector< double > z(1000, 0.);
    std::vector< double > h(1000, 0.);
    std::vector< double > m(1000, 0.);
    for (size_t i = 0; i < 1000; ++i) {
      x[i] = Utilities::random_double();
      y[i] = Utilities::random_double();
      z[i] = Utilities::random_double();
      h[i] = 0.2;
      m[i] = 0.001;
    }

    density_function.reset(x.data(), y.data(), z.data(), h.data(), m.data(),
                           1000);
    density_function.initialize();

    grid.initialize(block, density_function);

    cmac_status("Ntot: %g.", grid.get_total_hydrogen_number());
    assert_values_equal_rel(grid.get_total_hydrogen_number(), 5.25878e26,
                            1.e-6);
    AsciiFileDensityGridWriter writer("test_SPH_array_density_function_double",
                                      grid, ".");
    writer.write(0, params);
  }

  /// double precision positions, single precision smoothing lengths and masses
  {
    std::vector< double > x(1000, 0.);
    std::vector< double > y(1000, 0.);
    std::vector< double > z(1000, 0.);
    std::vector< float > h(1000, 0.);
    std::vector< float > m(1000, 0.);
    for (size_t i = 0; i < 1000; ++i) {
      x[i] = Utilities::random_double();
      y[i] = Utilities::random_double();
      z[i] = Utilities::random_double();
      h[i] = 0.2;
      m[i] = 0.001;
    }

    density_function.reset(x.data(), y.data(), z.data(), h.data(), m.data(),
                           1000);
    density_function.initialize();

    grid.initialize(block, density_function);

    cmac_status("Ntot: %g.", grid.get_total_hydrogen_number());
    assert_values_equal_rel(grid.get_total_hydrogen_number(), 5.19135e26,
                            1.e-6);
    AsciiFileDensityGridWriter writer("test_SPH_array_density_function_dfloat",
                                      grid, ".");
    writer.write(0, params);
  }

  /// single precision arrays
  {
    std::vector< float > x(1000, 0.);
    std::vector< float > y(1000, 0.);
    std::vector< float > z(1000, 0.);
    std::vector< float > h(1000, 0.);
    std::vector< float > m(1000, 0.);
    for (size_t i = 0; i < 1000; ++i) {
      x[i] = Utilities::random_double();
      y[i] = Utilities::random_double();
      z[i] = Utilities::random_double();
      h[i] = 0.2;
      m[i] = 0.001;
    }

    density_function.reset(x.data(), y.data(), z.data(), h.data(), m.data(),
                           1000);
    density_function.initialize();

    grid.initialize(block, density_function);

    cmac_status("Ntot: %g.", grid.get_total_hydrogen_number());
    assert_values_equal_rel(grid.get_total_hydrogen_number(), 5.23194e26,
                            1.e-7);
    AsciiFileDensityGridWriter writer("test_SPH_array_density_function_float",
                                      grid, ".");
    writer.write(0, params);
  }

  return 0;
}
