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
 * @param read_cosmic_ray_heating Read the variables for cosmic ray heating?
 * @param log Log to write logging info to.
 */
FLASHSnapshotDensityFunction::FLASHSnapshotDensityFunction(
    std::string filename, double temperature, bool read_cosmic_ray_heating,
    Log *log)
    : _read_cosmic_ray_heating(read_cosmic_ray_heating), _log(log) {

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
  HDF5Tools::HDF5Dictionary< int32_t > integer_runtime_pars =
      HDF5Tools::read_dictionary< int32_t >(file, "integer runtime parameters");
  CoordinateVector< uint_fast32_t > nblock;
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
  std::vector< int32_t > levels =
      HDF5Tools::read_dataset< int32_t >(file, "refine level");
  // read the node types
  std::vector< int32_t > nodetypes =
      HDF5Tools::read_dataset< int32_t >(file, "node type");
  // determine the level of each block
  uint_fast8_t level = 0;
  uint_fast32_t dsize = densities.size()[1];
  while (dsize > 1) {
    ++level;
    dsize >>= 1;
  }
  // add them to the grid
  for (size_t i = 0; i < extents.size()[0]; ++i) {
    if (nodetypes[i] == 1) {
      CoordinateVector<> anchor;
      std::array< size_t, 3 > ix0 = {{i, 0, 0}};
      anchor[0] = extents[ix0] * unit_length_in_SI;
      std::array< size_t, 3 > iy0 = {{i, 1, 0}};
      anchor[1] = extents[iy0] * unit_length_in_SI;
      std::array< size_t, 3 > iz0 = {{i, 2, 0}};
      anchor[2] = extents[iz0] * unit_length_in_SI;
      CoordinateVector<> top_anchor;
      std::array< size_t, 3 > ix1 = {{i, 0, 1}};
      top_anchor[0] = extents[ix1] * unit_length_in_SI;
      std::array< size_t, 3 > iy1 = {{i, 1, 1}};
      top_anchor[1] = extents[iy1] * unit_length_in_SI;
      std::array< size_t, 3 > iz1 = {{i, 2, 1}};
      top_anchor[2] = extents[iz1] * unit_length_in_SI;
      CoordinateVector<> sides = top_anchor - anchor;
      for (size_t ix = 0; ix < densities.size()[1]; ++ix) {
        for (size_t iy = 0; iy < densities.size()[2]; ++iy) {
          for (size_t iz = 0; iz < densities.size()[3]; ++iz) {
            CoordinateVector<> centre;
            centre[0] =
                anchor.x() + (ix + 0.5) * sides.x() / densities.size()[1];
            centre[1] =
                anchor.y() + (iy + 0.5) * sides.y() / densities.size()[2];
            centre[2] =
                anchor.z() + (iz + 0.5) * sides.z() / densities.size()[3];
            // this is the ordering as it is in the file
            std::array< size_t, 4 > irho = {{i, iz, iy, ix}};
            double rho = densities[irho];
            // each block contains level^3 cells, hence levels[i] + level
            // (but levels[i] is 1 larger than in our definition, Fortran counts
            // from 1)
            amrkey_t key = _grid.get_key(levels[i] + level - 1, centre);
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

  if (_read_cosmic_ray_heating) {
    // units
    // 1 Gauss = 1e-4 Tesla
    double unit_magnetic_field_in_SI = 1.e-4;
    // 1 erg g^-1 = 1e-4 J kg^-1
    double unit_cosmic_ray_energy_in_SI = 1.e-4;

    // read the magnetic field
    HDF5Tools::HDF5DataBlock< double, 4 > magnetic_field_x =
        HDF5Tools::read_dataset< double, 4 >(file, "magx");
    HDF5Tools::HDF5DataBlock< double, 4 > magnetic_field_y =
        HDF5Tools::read_dataset< double, 4 >(file, "magy");
    HDF5Tools::HDF5DataBlock< double, 4 > magnetic_field_z =
        HDF5Tools::read_dataset< double, 4 >(file, "magz");
    // read the cosmic ray energy
    HDF5Tools::HDF5DataBlock< double, 4 > cosmic_ray_energy =
        HDF5Tools::read_dataset< double, 4 >(file, "encr");
    for (size_t i = 0; i < extents.size()[0]; ++i) {
      if (nodetypes[i] == 1) {
        CoordinateVector<> anchor;
        std::array< size_t, 3 > ix0 = {{i, 0, 0}};
        anchor[0] = extents[ix0] * unit_length_in_SI;
        std::array< size_t, 3 > iy0 = {{i, 1, 0}};
        anchor[1] = extents[iy0] * unit_length_in_SI;
        std::array< size_t, 3 > iz0 = {{i, 2, 0}};
        anchor[2] = extents[iz0] * unit_length_in_SI;
        CoordinateVector<> top_anchor;
        std::array< size_t, 3 > ix1 = {{i, 0, 1}};
        top_anchor[0] = extents[ix1] * unit_length_in_SI;
        std::array< size_t, 3 > iy1 = {{i, 1, 1}};
        top_anchor[1] = extents[iy1] * unit_length_in_SI;
        std::array< size_t, 3 > iz1 = {{i, 2, 1}};
        top_anchor[2] = extents[iz1] * unit_length_in_SI;
        CoordinateVector<> sides = top_anchor - anchor;
        for (size_t ix = 0; ix < densities.size()[1]; ++ix) {
          for (size_t iy = 0; iy < densities.size()[2]; ++iy) {
            for (size_t iz = 0; iz < densities.size()[3]; ++iz) {
              CoordinateVector<> centre;
              centre[0] =
                  anchor.x() + (ix + 0.5) * sides.x() / densities.size()[1];
              centre[1] =
                  anchor.y() + (iy + 0.5) * sides.y() / densities.size()[2];
              centre[2] =
                  anchor.z() + (iz + 0.5) * sides.z() / densities.size()[3];
              // this is the ordering as it is in the file
              std::array< size_t, 4 > irho = {{i, iz, iy, ix}};
              // each block contains level^3 cells, hence levels[i] + level
              // (but levels[i] is 1 larger than in our definition, Fortran
              // counts
              // from 1)
              amrkey_t key = _grid.get_key(levels[i] + level - 1, centre);
              DensityValues &vals = _grid[key].value();
              vals.set_magnetic_field(
                  CoordinateVector<>(magnetic_field_x[irho],
                                     magnetic_field_y[irho],
                                     magnetic_field_z[irho]) *
                  unit_magnetic_field_in_SI);
              vals.set_cosmic_ray_energy(cosmic_ray_energy[irho] *
                                         unit_cosmic_ray_energy_in_SI);
            }
          }
        }
      }
    }
    // make sure the neighbour information for the grid is set
    _grid.set_ngbs(CoordinateVector< bool >(true, true, false));
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
 * Parameters are:
 *  - filename: Name of the snapshot file (required)
 *  - temperature: Temperature value used to initialize the cells (default: read
 *    temperature from snapshot file)
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
          params.get_value< bool >("DensityFunction:read cosmic ray heating",
                                   false),
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

  if (_read_cosmic_ray_heating) {
    amrkey_t key = _grid.get_key(position);
    AMRGridCell< DensityValues > &cell = _grid[key];
    AMRGridCell< DensityValues > *cell_left = cell.get_ngb(AMRNGBPOSITION_LEFT);
    while (!cell_left->is_single_cell()) {
      cell_left = cell_left->get_child(position);
    }
    AMRGridCell< DensityValues > *cell_right =
        cell.get_ngb(AMRNGBPOSITION_RIGHT);
    while (!cell_right->is_single_cell()) {
      cell_right = cell_right->get_child(position);
    }
    AMRGridCell< DensityValues > *cell_front =
        cell.get_ngb(AMRNGBPOSITION_FRONT);
    while (!cell_front->is_single_cell()) {
      cell_front = cell_front->get_child(position);
    }
    AMRGridCell< DensityValues > *cell_back = cell.get_ngb(AMRNGBPOSITION_BACK);
    while (!cell_back->is_single_cell()) {
      cell_back = cell_back->get_child(position);
    }
    AMRGridCell< DensityValues > *cell_bottom =
        cell.get_ngb(AMRNGBPOSITION_BOTTOM);
    while (!cell_bottom->is_single_cell()) {
      cell_bottom = cell_bottom->get_child(position);
    }
    AMRGridCell< DensityValues > *cell_top = cell.get_ngb(AMRNGBPOSITION_TOP);
    while (!cell_top->is_single_cell()) {
      cell_top = cell_top->get_child(position);
    }
    CoordinateVector<> gradPc;
    gradPc[0] =
        (cell_left->value().get_cosmic_ray_energy() -
         cell_right->value().get_cosmic_ray_energy()) /
        (cell_left->get_midpoint().x() - cell_right->get_midpoint().x());
    gradPc[1] =
        (cell_front->value().get_cosmic_ray_energy() -
         cell_back->value().get_cosmic_ray_energy()) /
        (cell_front->get_midpoint().y() - cell_back->get_midpoint().y());
    gradPc[2] =
        (cell_bottom->value().get_cosmic_ray_energy() -
         cell_top->value().get_cosmic_ray_energy()) /
        (cell_bottom->get_midpoint().z() - cell_top->get_midpoint().z());
    CoordinateVector<> B = cell.value().get_magnetic_field();
    values.set_cosmic_ray_factor(
        std::abs(CoordinateVector<>::dot_product(B, gradPc)));
  } else {
    values.set_cosmic_ray_factor(-1.);
  }

  return values;
}
