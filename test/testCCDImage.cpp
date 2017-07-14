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
 * @file testCCDImage.cpp
 *
 * @brief Unit test for the CCDImage class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CCDImage.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Unit test for the CCDImage class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  RandomGenerator random_generator(42);
  CCDImage image(0., 0., 1000, 1000, 0., 0., 1., 1., "PGM", "test_ccdimage",
                 ".");

  for (unsigned int i = 0; i < 1000000; ++i) {
    const double r =
        0.5 * std::sqrt(random_generator.get_uniform_random_double());
    const double cost = 2. * random_generator.get_uniform_random_double() - 1.;
    const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
    const double phi = 2. * M_PI * random_generator.get_uniform_random_double();
    const double cosp = std::cos(phi);
    const double sinp = std::sin(phi);
    const CoordinateVector<> position(r * sint * cosp + 0.5,
                                      r * sint * sinp + 0.5, r * cost + 0.5);
    image.add_photon(position, 1., 0., 0.);
  }

  CCDImage image2(0., 0., 1000, 1000, 0., 0., 1., 1., "PGM", "test_ccdimage",
                  ".");
  for (unsigned int i = 0; i < 100000; ++i) {
    const double r =
        0.2 * std::sqrt(random_generator.get_uniform_random_double());
    const double cost = 2. * random_generator.get_uniform_random_double() - 1.;
    const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
    const double phi = 2. * M_PI * random_generator.get_uniform_random_double();
    const double cosp = std::cos(phi);
    const double sinp = std::sin(phi);
    const CoordinateVector<> position(r * sint * cosp + 0.3,
                                      r * sint * sinp + 0.24, r * cost + 0.54);
    image2.add_photon(position, 1., 0., 0.);
  }

  image += image2;

  image.save();

  return 0;
}
