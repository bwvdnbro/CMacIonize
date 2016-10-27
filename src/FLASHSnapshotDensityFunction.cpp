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

  // read the block extents
  HDF5Tools::HDF5DataBlock< double, 3 > extents =
      HDF5Tools::read_dataset< double, 3 >(file, "bounding box");
  // read the densities
  HDF5Tools::HDF5DataBlock< double, 4 > densities =
      HDF5Tools::read_dataset< double, 4 >(file, "dens");
  // read the refinement levels
  std::vector< int > levels =
      HDF5Tools::read_dataset< int >(file, "refine level");
  // read the node types
  std::vector< int > nodetypes =
      HDF5Tools::read_dataset< int >(file, "node type");
  // determine the level of each block
  unsigned int level = 0;
  unsigned int dsize = densities.size()[1];
  while(dsize > 1){
    ++level;
    dsize >>= 1;
  }
  // add them to the grid
  for (unsigned int i = 0; i < extents.size()[0]; ++i) {
    if (nodetypes[i] == 1) {
      CoordinateVector<> anchor;
      anchor[0] = extents[{i, 0, 0}];
      anchor[1] = extents[{i, 1, 0}];
      anchor[2] = extents[{i, 2, 0}];
      CoordinateVector<> top_anchor;
      top_anchor[0] = extents[{i, 0, 1}];
      top_anchor[1] = extents[{i, 1, 1}];
      top_anchor[2] = extents[{i, 2, 1}];
      CoordinateVector<> sides = top_anchor - anchor;
      for (unsigned int ix = 0; ix < densities.size()[1]; ++ix) {
        for (unsigned int iy = 0; iy < densities.size()[2]; ++iy) {
          for (unsigned int iz = 0; iz < densities.size()[3]; ++iz) {
            CoordinateVector<> centre;
            centre[0] =
                anchor.x() + (ix + 0.5) * sides.x() / densities.size()[1];
            centre[1] =
                anchor.y() + (iy + 0.5) * sides.y() / densities.size()[2];
            centre[2] =
                anchor.z() + (iz + 0.5) * sides.z() / densities.size()[3];
            // this is the ordering as it is in the file
            double rho = densities[{i, iz, iy, ix}];
            // each block contains level^3 cells, hence levels[i] + level
            // (but levels[i] is 1 larger than in our definition, Fortran counts
            // from 1)
            unsigned long key = _grid.get_key(levels[i] + level -1, centre);
            _grid.create_cell(key) = rho;
          }
        }
      }
    }
  }

  HDF5Tools::close_file(file);
}

/**
 * @brief Function that returns the density at the given coordinate position.
 *
 * @param position CoordinateVector<> specifying a position.
 * @return Density at that position.
 */
double FlashSnapshotDensityFunction::operator()(CoordinateVector<> position) {
  return _grid.get_cell(position);
}
