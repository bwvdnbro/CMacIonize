/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AmunSnapshotDensityFunction.cpp
 *
 * @brief DensityFunction that reads a density field from an Amun HDF5 snapshot
 * file: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AmunSnapshotDensityFunction.hpp"
#include "HDF5Tools.hpp"
#include "ParameterFile.hpp"
#include "Utilities.hpp"

/**
 * @brief Constructor.
 *
 * @param folder Folder that contains the snapshots.
 * @param prefix Prefix for the snapshots.
 * @param padding Number of padding digits in snapshot names.
 * @param number_of_files Number of files that make up the snapshot.
 * @param box Desired dimensions of the simulation box (in m).
 * @param number_density Desired average number density (in m^-3).
 * @param shift Position shift (in fractions of the box size).
 */
AmunSnapshotDensityFunction::AmunSnapshotDensityFunction(
    const std::string folder, const std::string prefix,
    const uint_fast32_t padding, const uint_fast32_t number_of_files,
    const Box<> box, const double number_density,
    const CoordinateVector<> shift)
    : _box(box), _shift(shift) {

  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the first file for metadata reading
  std::vector< int32_t > dims;
  std::vector< int32_t > pdims;
  {
    const std::string name =
        Utilities::compose_filename(folder, prefix, "h5", 0, padding);
    HDF5Tools::HDF5File file =
        HDF5Tools::open_file(name, HDF5Tools::HDF5FILEMODE_READ);
    HDF5Tools::HDF5Group attributes =
        HDF5Tools::open_group(file, "/attributes");
    dims = HDF5Tools::read_vector_attribute< int32_t >(attributes, "dims");
    pdims = HDF5Tools::read_vector_attribute< int32_t >(attributes, "pdims");

    HDF5Tools::close_group(attributes);
    HDF5Tools::close_file(file);

    _number_of_cells[0] = static_cast< uint_fast32_t >(dims[0] * pdims[0]);
    _number_of_cells[1] = static_cast< uint_fast32_t >(dims[1] * pdims[1]);
    _number_of_cells[2] = static_cast< uint_fast32_t >(dims[2] * pdims[2]);
  }
  const uint_fast32_t totnumcell =
      _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2];

  _densities.resize(totnumcell);

  // now read all the blocks
  double average_density = 0.;
  for (uint_fast32_t ifile = 0; ifile < number_of_files; ++ifile) {
    uint_fast32_t offset_z = ifile / (pdims[0] * pdims[1]);
    uint_fast32_t offset_x =
        (ifile - offset_z * pdims[0] * pdims[1]) / pdims[1];
    uint_fast32_t offset_y =
        ifile - offset_z * pdims[0] * pdims[1] - offset_x * pdims[1];

    offset_x *= dims[0];
    offset_y *= dims[1];
    offset_z *= dims[2];

    const std::string name =
        Utilities::compose_filename(folder, prefix, "h5", ifile, padding);
    HDF5Tools::HDF5File file =
        HDF5Tools::open_file(name, HDF5Tools::HDF5FILEMODE_READ);

    {
      HDF5Tools::HDF5Group variables =
          HDF5Tools::open_group(file, "/variables");
      HDF5Tools::HDF5DataBlock< float, 3 > dens =
          HDF5Tools::read_dataset< float, 3 >(variables, "dens");
      HDF5Tools::close_group(variables);

      for (int_fast32_t ix = 0; ix < dims[0]; ++ix) {
        for (int_fast32_t iy = 0; iy < dims[1]; ++iy) {
          for (int_fast32_t iz = 0; iz < dims[2]; ++iz) {
            std::array< size_t, 3 > index = {{static_cast< size_t >(iz),
                                              static_cast< size_t >(iy),
                                              static_cast< size_t >(ix)}};
            const double this_density = dens[index];
            _densities[(iz + offset_z) * _number_of_cells[1] *
                           _number_of_cells[0] +
                       (iy + offset_y) * _number_of_cells[0] + ix + offset_x] =
                this_density;
            average_density += this_density;
          }
        }
      }
    }

    HDF5Tools::close_file(file);
  }
  average_density /= totnumcell;

  const double conversion_factor = number_density / average_density;
  for (size_t i = 0; i < _densities.size(); ++i) {
    _densities[i] *= conversion_factor;
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * The following parameters are read from the parameter file:
 *  - folder: Folder containing the snapshot(s) (default: .)
 *  - prefix: Prefix for the snapshot names (required)
 *  - padding: Number of digits of padding in each snapshot name (default: 5)
 *  - number of files: Number of files that make up the snapshot (default: 4)
 *  - box anchor: Anchor of the desired simulation box (default: [0. m, 0. m,
 *    0. m])
 *  - box sides: Sides lengths of the desired simulation box (default: [1. m, 1.
 *    m, 1. m])
 *  - average number density: Desired average number density (default: 100.
 *    cm^-3)
 *  - shift: (Periodic) shift to apply to all positions (in fractions of the box
 *    size, default: [0., 0., 0.])
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
AmunSnapshotDensityFunction::AmunSnapshotDensityFunction(ParameterFile &params,
                                                         Log *log)
    : AmunSnapshotDensityFunction(
          params.get_value< std::string >("DensityFunction:folder", "."),
          params.get_value< std::string >("DensityFunction:prefix"),
          params.get_value< uint_fast32_t >("DensityFunction:padding", 5),
          params.get_value< uint_fast32_t >("DensityFunction:number of files",
                                            4),
          Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                    "DensityFunction:box anchor", "[0. m, 0. m, 0. m]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "DensityFunction:box sides", "[1. m, 1. m, 1. m]")),
          params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
              "DensityFunction:average number density", "100. cm^-3"),
          params.get_value< CoordinateVector<> >("DensityFunction:shift",
                                                 CoordinateVector<>(0.))) {}

/**
 * @brief Virtual destructor.
 */
AmunSnapshotDensityFunction::~AmunSnapshotDensityFunction() {}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues AmunSnapshotDensityFunction::operator()(const Cell &cell) const {

  const CoordinateVector<> midpoint = cell.get_cell_midpoint();
  CoordinateVector<> dx = midpoint - _box.get_anchor();
  for (uint_fast8_t i = 0; i < 3; ++i) {
    dx[i] -= _shift[i] * _box.get_sides()[i];
    while (dx[i] >= _box.get_sides()[i]) {
      dx[i] -= _box.get_sides()[i];
    }
    while (dx[i] < 0.) {
      dx[i] += _box.get_sides()[i];
    }
  }

  const uint_fast32_t ix = dx.x() / _box.get_sides().x() * _number_of_cells.x();
  const uint_fast32_t iy = dx.y() / _box.get_sides().y() * _number_of_cells.y();
  const uint_fast32_t iz = dx.z() / _box.get_sides().z() * _number_of_cells.z();

  const double nH = _densities[iz * _number_of_cells[1] * _number_of_cells[0] +
                               iy * _number_of_cells[0] + ix];

  DensityValues values;
  values.set_number_density(nH);

  return values;
}
