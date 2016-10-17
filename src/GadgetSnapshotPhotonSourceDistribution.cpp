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
 * @file GadgetSnapshotPhotonSourceDistribution.cpp
 *
 * @brief GadgetSnapshotPhotonSourceDistribution implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "GadgetSnapshotPhotonSourceDistribution.hpp"
#include "HDF5Tools.hpp"

/**
 * @brief Constructor.
 *
 * Reads in the sources from the file and stores them in internal arrays.
 *
 * @param filename Name of the snapshot file to read.
 */
GadgetSnapshotPhotonSourceDistribution::GadgetSnapshotPhotonSourceDistribution(
    std::string filename) {
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);
  // open the group containing the star particle data
  HDF5Tools::HDF5Group starparticles =
      HDF5Tools::open_group(file, "/PartType4");
  // read the positions
  _positions = HDF5Tools::read_dataset< CoordinateVector<> >(starparticles,
                                                             "Coordinates");
  // close the group
  HDF5Tools::close_group(starparticles);
  // close the file
  HDF5Tools::close_file(file);
}

/**
 * @brief Get the number of sources in the snapshot file.
 *
 * @return Number of sources.
 */
unsigned int GadgetSnapshotPhotonSourceDistribution::get_number_of_sources() {
  return _positions.size();
}

/**
 * @brief Get the position of one of the sources.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Position of the given source.
 */
CoordinateVector<>
GadgetSnapshotPhotonSourceDistribution::get_position(unsigned int index) {
  return _positions[index];
}

/**
 * @brief Get the weight of one of the sources.
 *
 * At the moment, all sources have an equal weight.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Weight of the given source.
 */
double GadgetSnapshotPhotonSourceDistribution::get_weight(unsigned int index) {
  return 1. / _positions.size();
}
