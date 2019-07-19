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
 * @file testFLASHSnapshotDensityFunction.cpp
 *
 * @brief Unit test for the FLASHSnapshotDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "FLASHSnapshotDensityFunction.hpp"

/**
 * @brief Expected density function.
 *
 * @param x CoordinateVector<> specifying a position.
 * @return Expected density at that position.
 */
double expected_density(CoordinateVector<> x) {
  return (1. + 100. * x.x() + 100. * x.y() + 100. * x.z()) * 1.e3 /
         1.6737236e-27;
}

/**
 * @brief Unit test for the FLASHSnapshotDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  FLASHSnapshotDensityFunction density("FLASHtest.hdf5");
  density.initialize();

  const uint_fast32_t np = 128;
  double xi2 = 0.;
  for (uint_fast32_t i = 0; i < np; ++i) {
    CoordinateVector<> p((i + 0.5) * 0.02 / np, (i + 0.5) * 0.01 / np,
                         (i + 0.5) * 0.01 / np);
    DummyCell cell(p.x(), p.y(), p.z());
    DensityValues vals = density(cell);
    double rho = vals.get_number_density();
    double rho_ex = expected_density(p);
    double diff = (rho - rho_ex) / (rho + rho_ex);
    xi2 += diff * diff;

    assert_condition(vals.get_temperature() == 4000.);
  }
  assert_values_equal(xi2, 0.0106294);

  return 0;
}
