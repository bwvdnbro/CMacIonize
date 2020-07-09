/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 * Copyright (C) 2020 Nina Sartorio (sartorio.ninae@gmail.com)
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
 * @file testAbsorptionTracker.cpp
 *
 * @brief Unit test for the AbsorptionTracker class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Nina Sartorio (sartorio.nina@gmail.com)
 */

#include "AbsorptionTracker.hpp"
#include "Assert.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Unit test for the AbsorptionTracker class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  AbsorptionTracker tracker;

  RandomGenerator rg;

  PhotonPacket photon;
  photon.set_type(PHOTONTYPE_PRIMARY);
  double absorption[NUMBER_OF_IONNAMES];
  for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    absorption[i] = rg.get_uniform_random_double();
  }
  tracker.count_photon(photon, absorption);

  tracker.output_tracker("test_absorption_tracker.txt");

#ifdef HAVE_HDF5
  HDF5Tools::HDF5File file = HDF5Tools::open_file(
      "test_absorption_tracker.hdf5", HDF5Tools::HDF5FILEMODE_WRITE);
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "AbsorptionGroup");
  tracker.create_group(group, 1);
  tracker.append_to_group(group, 0);
  HDF5Tools::close_group(group);
  HDF5Tools::close_file(file);
#endif

  return 0;
}
