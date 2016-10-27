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

  // find out the dimensions of the box
  HDF5Tools::HDF5Dictionary< double > real_runtime_pars =
      HDF5Tools::read_dictionary< double >(file, "real runtime parameters");
  CoordinateVector<> anchor;
  anchor[0] = real_runtime_pars["xmin"];
  anchor[1] = real_runtime_pars["ymin"];
  anchor[2] = real_runtime_pars["zmin"];
  CoordinateVector<> top_anchor;
  top_anchor[0] = real_runtime_pars["xmax"];
  top_anchor[1] = real_runtime_pars["ymax"];
  top_anchor[2] = real_runtime_pars["zmax"];
  CoordinateVector<> sides = top_anchor - anchor;
  Box box(anchor, sides);

  // find out the number of blocks in each dimension
  HDF5Tools::HDF5Dictionary< int > integer_runtime_pars =
      HDF5Tools::read_dictionary< int >(file, "integer runtime parameters");
  CoordinateVector< int > nblock;
  nblock[0] = integer_runtime_pars["nblockx"];
  nblock[1] = integer_runtime_pars["nblocky"];
  nblock[2] = integer_runtime_pars["nblockz"];

  // make the grid
  _grid = AMRGrid< double >(box, nblock);

  // fill the grid with values

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
