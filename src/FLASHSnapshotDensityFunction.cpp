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
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"

/**
 * @brief Constructor.
 *
 * This reads in the data and stores it internally.
 *
 * @param filename Name of the snapshot file to read.
 * @param temperature Initial temperature for the ISM (in K).
 * @param log Log to write logging info to.
 */
FLASHSnapshotDensityFunction::FLASHSnapshotDensityFunction(std::string filename,
                                                           double temperature,
                                                           Log *log)
    : _log(log) {
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

  // units
  double unit_length_in_SI = UnitConverter::to_SI< QUANTITY_LENGTH >(1., "cm");
  double unit_density_in_SI =
      UnitConverter::to_SI< QUANTITY_DENSITY >(1., "g cm^-3");
  // temperatures are already in K
  double unit_temperature_in_SI = 1.;

  // find out the dimensions of the box
  HDF5Tools::HDF5Dictionary< double > real_runtime_pars =
      HDF5Tools::read_dictionary< double >(file, "real runtime parameters");
  CoordinateVector<> anchor;
  anchor[0] = real_runtime_pars["xmin"] * unit_length_in_SI;
  anchor[1] = real_runtime_pars["ymin"] * unit_length_in_SI;
  anchor[2] = real_runtime_pars["zmin"] * unit_length_in_SI;
  CoordinateVector<> top_anchor;
  top_anchor[0] = real_runtime_pars["xmax"] * unit_length_in_SI;
  top_anchor[1] = real_runtime_pars["ymax"] * unit_length_in_SI;
  top_anchor[2] = real_runtime_pars["zmax"] * unit_length_in_SI;
  CoordinateVector<> sides = top_anchor - anchor;
  Box<> box(anchor, sides);

  // find out the number of blocks in each dimension
  HDF5Tools::HDF5Dictionary< int > integer_runtime_pars =
      HDF5Tools::read_dictionary< int >(file, "integer runtime parameters");
  CoordinateVector< int > nblock;
  nblock[0] = integer_runtime_pars["nblockx"];
  nblock[1] = integer_runtime_pars["nblocky"];
  nblock[2] = integer_runtime_pars["nblockz"];

  // make the grid
  _grid = AMRGrid< DensityValues >(box, nblock);

  // fill the grid with values

  // read the block extents
  HDF5Tools::HDF5DataBlock< double, 3 > extents =
      HDF5Tools::read_dataset< double, 3 >(file, "bounding box");
  // read the densities
  HDF5Tools::HDF5DataBlock< double, 4 > densities =
      HDF5Tools::read_dataset< double, 4 >(file, "dens");
  // read the temperatures
  HDF5Tools::HDF5DataBlock< double, 4 > temperatures =
      HDF5Tools::read_dataset< double, 4 >(file, "temp");
  // read the refinement levels
  std::vector< int > levels =
      HDF5Tools::read_dataset< int >(file, "refine level");
  // read the node types
  std::vector< int > nodetypes =
      HDF5Tools::read_dataset< int >(file, "node type");
  // determine the level of each block
  unsigned int level = 0;
  unsigned int dsize = densities.size()[1];
  while (dsize > 1) {
    ++level;
    dsize >>= 1;
  }
  // add them to the grid
  for (unsigned int i = 0; i < extents.size()[0]; ++i) {
    if (nodetypes[i] == 1) {
      CoordinateVector<> anchor;
      std::array< unsigned int, 3 > ix0 = {{i, 0, 0}};
      anchor[0] = extents[ix0] * unit_length_in_SI;
      std::array< unsigned int, 3 > iy0 = {{i, 1, 0}};
      anchor[1] = extents[iy0] * unit_length_in_SI;
      std::array< unsigned int, 3 > iz0 = {{i, 2, 0}};
      anchor[2] = extents[iz0] * unit_length_in_SI;
      CoordinateVector<> top_anchor;
      std::array< unsigned int, 3 > ix1 = {{i, 0, 1}};
      top_anchor[0] = extents[ix1] * unit_length_in_SI;
      std::array< unsigned int, 3 > iy1 = {{i, 1, 1}};
      top_anchor[1] = extents[iy1] * unit_length_in_SI;
      std::array< unsigned int, 3 > iz1 = {{i, 2, 1}};
      top_anchor[2] = extents[iz1] * unit_length_in_SI;
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
            std::array< unsigned int, 4 > irho = {{i, iz, iy, ix}};
            double rho = densities[irho];
            // each block contains level^3 cells, hence levels[i] + level
            // (but levels[i] is 1 larger than in our definition, Fortran counts
            // from 1)
            unsigned long key = _grid.get_key(levels[i] + level - 1, centre);
            DensityValues &vals = _grid.create_cell(key);
            vals.set_number_density(rho * unit_density_in_SI);
            if (temperature <= 0.) {
              double temp = temperatures[irho];
              vals.set_temperature(temp * unit_temperature_in_SI);
            } else {
              vals.set_temperature(temperature);
            }
          }
        }
      }
    }
  }

  HDF5Tools::close_file(file);

  if (_log) {
    _log->write_status("Successfully read densities from file \"", filename,
                       "\".");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read.
 * @param log Log to write logging info to.
 */
FLASHSnapshotDensityFunction::FLASHSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : FLASHSnapshotDensityFunction(
          params.get_value< std::string >("DensityFunction:filename"),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "DensityFunction:temperature", "-1. K"),
          log) {}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues FLASHSnapshotDensityFunction::operator()(const Cell &cell) const {

  DensityValues values;

  const CoordinateVector<> position = cell.get_cell_midpoint();

  const DensityValues &vals = _grid.get_cell(position);
  values.set_number_density(vals.get_number_density() / 1.6737236e-27);
  values.set_temperature(vals.get_temperature());
  values.set_ionic_fraction(ION_H_n, 1.e-6);
  values.set_ionic_fraction(ION_He_n, 1.e-6);
  return values;
}
