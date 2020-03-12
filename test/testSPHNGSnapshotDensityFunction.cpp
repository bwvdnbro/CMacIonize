/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testSPHNGSnapshotDensityFunction.cpp
 *
 * @brief Unit test for the SPHNGSnapshotDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "SPHNGSnapshotDensityFunction.hpp"
#include "SPHNGVoronoiGeneratorDistribution.hpp"
#include "UnitConverter.hpp"

#include <fstream>
#include <sstream>

/**
 * @brief Unit test for the SPHNGSnapshotDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// tagged file
  {
    cmac_status("Start reading tagged file...");
    SPHNGSnapshotDensityFunction density_function("SPHNGtest.dat", 8000., false,
                                                  0, 0., 0., "");
    density_function.initialize();

    std::vector< CoordinateVector<> > positions;
    std::ifstream file("SPHNG_data.txt");
    std::string line;
    uint_fast32_t index = 0;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double x, y, z, m, h;
      lstream >> x >> y >> z >> m >> h;

      // convert units
      x = UnitConverter::to_SI< QUANTITY_LENGTH >(x, "cm");
      y = UnitConverter::to_SI< QUANTITY_LENGTH >(y, "cm");
      z = UnitConverter::to_SI< QUANTITY_LENGTH >(z, "cm");
      m = UnitConverter::to_SI< QUANTITY_MASS >(m, "g");
      h = UnitConverter::to_SI< QUANTITY_LENGTH >(h, "cm");

      const double tolerance = 1.e-14;

      const CoordinateVector<> p = density_function.get_position(index);
      assert_values_equal_rel(x, p.x(), tolerance);
      assert_values_equal_rel(y, p.y(), tolerance);
      assert_values_equal_rel(z, p.z(), tolerance);

      assert_values_equal_rel(m, density_function.get_mass(index), tolerance);
      assert_values_equal_rel(h, density_function.get_smoothing_length(index),
                              tolerance);

      positions.push_back(p);

      ++index;
    }

    const Box<> box(0., 1.);
    std::sort(positions.begin(), positions.end(),
              SPHNGVoronoiGeneratorDistribution::position_smaller_than);
    SPHNGVoronoiGeneratorDistribution generator_distribution(
        box, "SPHNGtest.dat", nullptr);
    for (uint_fast32_t i = 0; i < positions.size(); ++i) {
      const CoordinateVector<> p = generator_distribution.get_position();
      assert_condition(positions[i].x() == p.x());
      assert_condition(positions[i].y() == p.y());
      assert_condition(positions[i].z() == p.z());
    }
    cmac_status("Done reading tagged file.");
  }

  /// untagged file
  {
    cmac_status("Start reading untagged file...");
    SPHNGSnapshotDensityFunction density_function("SPHNGtest_notags.dat", 8000.,
                                                  false, 0, 0., 0., "");
    density_function.initialize();

    std::vector< CoordinateVector<> > positions;
    std::ifstream file("SPHNG_data.txt");
    std::string line;
    uint_fast32_t index = 0;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double x, y, z, m, h;
      lstream >> x >> y >> z >> m >> h;

      // convert units
      x = UnitConverter::to_SI< QUANTITY_LENGTH >(x, "cm");
      y = UnitConverter::to_SI< QUANTITY_LENGTH >(y, "cm");
      z = UnitConverter::to_SI< QUANTITY_LENGTH >(z, "cm");
      m = UnitConverter::to_SI< QUANTITY_MASS >(m, "g");
      h = UnitConverter::to_SI< QUANTITY_LENGTH >(h, "cm");

      const double tolerance = 1.e-14;

      const CoordinateVector<> p = density_function.get_position(index);
      assert_values_equal_rel(x, p.x(), tolerance);
      assert_values_equal_rel(y, p.y(), tolerance);
      assert_values_equal_rel(z, p.z(), tolerance);

      assert_values_equal_rel(m, density_function.get_mass(index), tolerance);
      assert_values_equal_rel(h, density_function.get_smoothing_length(index),
                              tolerance);

      positions.push_back(p);

      ++index;
    }

    const Box<> box(0., 1.);
    std::sort(positions.begin(), positions.end(),
              SPHNGVoronoiGeneratorDistribution::position_smaller_than);
    SPHNGVoronoiGeneratorDistribution generator_distribution(
        box, "SPHNGtest_notags.dat", nullptr);
    for (uint_fast32_t i = 0; i < positions.size(); ++i) {
      const CoordinateVector<> p = generator_distribution.get_position();
      assert_condition(positions[i].x() == p.x());
      assert_condition(positions[i].y() == p.y());
      assert_condition(positions[i].z() == p.z());
    }
    cmac_status("Done reading untagged file.");
  }

  return 0;
}
