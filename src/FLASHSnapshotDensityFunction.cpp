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
 * @file FLASHSnapshotDensityFunction.cpp
 *
 * @brief FLASHSnapshotDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "FLASHSnapshotDensityFunction.hpp"
#include "HDF5Tools.hpp"

/**
 * @brief Constructor.
 *
 * This reads in the data and stores it internally.
 *
 * @param filename Name of the snapshot file to read.
 */
FlashSnapshotDensityFunction::FlashSnapshotDensityFunction(
    std::string filename) {
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

  HDF5Tools::close_file(file);
}

/**
 * @brief Function that returns the density at the given coordinate position.
 *
 * @param position CoordinateVector<> specifying a position.
 * @return Density at that position.
 */
double FlashSnapshotDensityFunction::operator()(CoordinateVector<> position) {
  return 0.;
}
