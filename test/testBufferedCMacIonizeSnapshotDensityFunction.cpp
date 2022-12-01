/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testBufferedCMacIonizeSnapshotDensityFunction.cpp
 *
 * @brief Unit test for the BufferedCMacIonizeSnapshotDensityFunction class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "BufferedCMacIonizeSnapshotDensityFunction.hpp"
#include "TerminalLog.hpp"

#include <fstream>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/**
 * @brief Unit test for the BufferedCMacIonizeSnapshotDensityFunction class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  TerminalLog log(LOGLEVEL_INFO);
  const double pc = UnitConverter::convert(1., "pc", "m");
  const Box<> box(-5. * pc, 10. * pc);
  const CoordinateVector< uint_fast32_t > ncell(8);
  BufferedCMacIonizeSnapshotDensityFunction density_function(
      "taskbased.hdf5", 10, false, box, ncell, &log);

  density_function.initialize();

  // we will test the density function twice
  // the first time, we do a serial pass, output some values to a text file for
  // manual checking, and store the neutral fractions for the second time
  std::ofstream ofile("test_buffered_snapshot.txt");
  ofile << "# x (m)\ty (m)\txH\tT (K)\n";
  std::vector< double > xH(ncell.x() * ncell.y());
  for (uint_fast32_t ix = 0; ix < ncell.x(); ++ix) {
    for (uint_fast32_t iy = 0; iy < ncell.y(); ++iy) {
      const CoordinateVector<> p(
          box.get_anchor().x() + (ix + 0.5) * box.get_sides().x() / ncell.x(),
          box.get_anchor().y() + (iy + 0.5) * box.get_sides().y() / ncell.y(),
          0.);
      const DummyCell cell(p.x(), p.y(), p.z());
      const DensityValues values = density_function(cell);
      assert_condition(values.get_number_density() == 1.e8);
      ofile << p.x() << "\t" << p.y() << "\t"
            << values.get_ionic_fraction(ION_H_n) << "\t"
            << values.get_temperature() << "\n";
      xH[ix * ncell.y() + iy] = values.get_ionic_fraction(ION_H_n);
    }
  }
  ofile.close();

  // for the second pass, we call the density function in parallel and check
  // that the neutral fractions have not changed
#ifdef HAVE_OPENMP
#pragma omp parallel for num_threads(4)
  for (uint_fast32_t ix = 0; ix < ncell.x(); ++ix) {
    for (uint_fast32_t iy = 0; iy < ncell.y(); ++iy) {
      const CoordinateVector<> p(
          box.get_anchor().x() + (ix + 0.5) * box.get_sides().x() / ncell.x(),
          box.get_anchor().y() + (iy + 0.5) * box.get_sides().y() / ncell.y(),
          0.);
      const DummyCell cell(p.x(), p.y(), p.z());
      const DensityValues values = density_function(cell);
      assert_condition(values.get_number_density() == 1.e8);
      assert_condition(xH[ix * ncell.y() + iy] ==
                       values.get_ionic_fraction(ION_H_n));
    }
  }
#endif

  // clean up
  density_function.free();

  return 0;
}
