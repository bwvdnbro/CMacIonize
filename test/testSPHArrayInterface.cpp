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
 * @file testSPHArrayInterface.cpp
 *
 * @brief Unit test for the SPHArrayInterface class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AsciiFileDensityGridWriter.hpp"
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "Octree.hpp"
#include "SPHArrayInterface.hpp"
#include "Utilities.hpp"

/**
 * @brief Test DensityFunction.
 */
class TestDensityFunction : public DensityFunction {
public:
  /**
   * @brief Get the neutral fraction for the given cell.
   *
   * @param cell Cell.
   * @return Neutral fraction.
   */
  virtual DensityValues operator()(const Cell &cell) {
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
 * @brief Unit test for the SPHArrayInterface class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CartesianDensityGrid grid(box, 32);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  ParameterFile params;
  SPHArrayInterface interface(1., 1., "Petkova");

  /// DensityFunction functionality

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

    interface.reset(x.data(), y.data(), z.data(), h.data(), m.data(), 1000);
    interface.initialize();

    grid.initialize(block, interface);

    cmac_status("Ntot: %g.", grid.get_total_hydrogen_number());
    assert_values_equal_rel(grid.get_total_hydrogen_number(), 1.68285e+27,
                            1.e-6);
    AsciiFileDensityGridWriter writer("test_SPH_array_density_function_double",
                                      ".");
    writer.write(grid, 0, params);
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

    interface.reset(x.data(), y.data(), z.data(), h.data(), m.data(), 1000);
    interface.initialize();

    grid.initialize(block, interface);

    cmac_status("Ntot: %g.", grid.get_total_hydrogen_number());
    assert_values_equal_rel(grid.get_total_hydrogen_number(), 1.6572e+27,
                            1.e-5);
    AsciiFileDensityGridWriter writer("test_SPH_array_density_function_dfloat",
                                      ".");
    writer.write(grid, 0, params);
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

    interface.reset(x.data(), y.data(), z.data(), h.data(), m.data(), 1000);
    interface.initialize();

    grid.initialize(block, interface);

    cmac_status("Ntot: %g.", grid.get_total_hydrogen_number());
    assert_values_equal_rel(grid.get_total_hydrogen_number(), 1.68544e+27,
                            1.e-6);
    AsciiFileDensityGridWriter writer("test_SPH_array_density_function_float",
                                      ".");
    writer.write(grid, 0, params);
  }

  /// DensityGridWriter functionality
  {
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    CartesianDensityGrid grid(box, 10);
    TestDensityFunction density_function;
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, density_function);

    std::vector< double > x(1000, 0.);
    std::vector< double > y(1000, 0.);
    std::vector< double > z(1000, 0.);
    std::vector< float > h(1000, 0.);
    std::vector< float > m(1000, 0.);
    for (uint_fast8_t ix = 0; ix < 10; ++ix) {
      for (uint_fast8_t iy = 0; iy < 10; ++iy) {
        for (uint_fast8_t iz = 0; iz < 10; ++iz) {
          const uint_fast32_t i = ix * 100 + iy * 10 + iz;
          x[i] = (ix + 0.5) * 0.1;
          y[i] = (iy + 0.5) * 0.1;
          z[i] = (iz + 0.5) * 0.1;
          h[i] = 0.2;
          m[i] = 0.001;
        }
      }
    }
    interface.reset(x.data(), y.data(), z.data(), h.data(), m.data(), 1000);
    interface.initialize();

    ParameterFile params;
    interface.write(grid, 0, params, 0.);

    std::vector< double > neutral_fractions(1000, -1.);
    interface.fill_array(neutral_fractions.data());

    // this does not work for the new inverse mapping algorithm
    //    for (uint_fast32_t i = 0; i < 1000; ++i) {
    //      const CoordinateVector<> p(x[i], y[i], z[i]);
    //      if ((p - CoordinateVector<>(0.5)).norm() < 0.2) {
    //        assert_condition(neutral_fractions[i] == 1.);
    //      } else {
    //        assert_condition(neutral_fractions[i] == 0.);
    //      }
    //    }
  }

  return 0;
}
