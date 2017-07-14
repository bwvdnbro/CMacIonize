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
 * @file testSpiralGalaxyContinuousPhotonSource.cpp
 *
 * @brief Unit test for the SpiralGalaxyContinuousPhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CCDImage.hpp"
#include "SpiralGalaxyContinuousPhotonSource.hpp"
#include <fstream>

/**
 * @brief Unit test for the SpiralGalaxyContinuousPhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  const double kpc = 3.086e19;
  const CoordinateVector<> box_anchor(-12. * kpc);
  const CoordinateVector<> box_sides(24. * kpc);
  const Box<> box(box_anchor, box_sides);
  const CoordinateVector< unsigned int > ncell(64);
  SpiralGalaxyContinuousPhotonSource source(box, 5. * kpc, 0.6 * kpc, 0.2);
  RandomGenerator random_generator(42);

  CCDImage image(0.5 * M_PI, 0., 200, 200, -12. * kpc, -12. * kpc, 24. * kpc,
                 24. * kpc, "BinaryArray",
                 "test_spiralgalaxycontinuousphotonsource_image", ".");
  const unsigned int ntot = ncell.x() * ncell.y() * ncell.z();
  std::vector< double > dens(ntot, 0.);
  for (unsigned int i = 0; i < 500000; ++i) {
    auto photon = source.get_random_incoming_direction(random_generator);
    const CoordinateVector<> p = photon.first;
    image.add_photon(p, 1., 0., 0.);
    const unsigned int ix =
        ncell.x() * (p.x() - box_anchor.x()) / box_sides.x();
    const unsigned int iy =
        ncell.y() * (p.y() - box_anchor.y()) / box_sides.y();
    const unsigned int iz =
        ncell.z() * (p.z() - box_anchor.z()) / box_sides.z();

    const unsigned int index = ix * ncell.y() * ncell.z() + iy * ncell.z() + iz;
    dens[index] += 1.;
  }

  image.save();

  std::ofstream ofile("test_spiralgalaxycontinuousphotonsource.txt");
  ofile << "#x (m)\ty (m)\tz (m)\tnumber of photons\n";
  for (unsigned int i = 0; i < ntot; ++i) {
    const unsigned int ix = i / ncell.y() / ncell.z();
    const unsigned int iy = (i - ix * ncell.y() * ncell.z()) / ncell.z();
    const unsigned int iz = i - ix * ncell.y() * ncell.z() - iy * ncell.z();
    const CoordinateVector<> p(
        box_anchor.x() + (ix + 0.5) * box_sides.x() / ncell.x(),
        box_anchor.y() + (iy + 0.5) * box_sides.y() / ncell.y(),
        box_anchor.z() + (iz + 0.5) * box_sides.z() / ncell.z());
    ofile << p.x() << "\t" << p.y() << "\t" << p.z() << "\t" << dens[i] << "\n";
  }
  ofile.close();

  return 0;
}
