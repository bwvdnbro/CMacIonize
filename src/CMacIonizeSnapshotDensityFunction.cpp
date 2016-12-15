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
 * @file CMacIonizeSnapshotDensityFunction.cpp
 *
 * @brief CMacIonizeSnapshotDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CMacIonizeSnapshotDensityFunction.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief Constructor.
 *
 * @param filename Name of the snapshot file to read.
 * @param log Log to write logging info to.
 */
CMacIonizeSnapshotDensityFunction::CMacIonizeSnapshotDensityFunction(
    std::string filename, Log *log) {
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

  // read grid parameters
  HDF5Tools::HDF5Group group = HDF5Tools::open_group(file, "/Parameters");
  std::vector< std::string > parameternames =
      HDF5Tools::get_attribute_names(group);
  ParameterFile parameters;
  for (auto it = parameternames.begin(); it != parameternames.end(); ++it) {
    std::string attname = *it;
    std::string attvalue =
        HDF5Tools::read_attribute< std::string >(group, attname);
    parameters.add_value(attname, attvalue);
  }
  _box = Box(parameters.get_physical_vector< QUANTITY_LENGTH >(
                 "densitygrid.box_anchor"),
             parameters.get_physical_vector< QUANTITY_LENGTH >(
                 "densitygrid.box_sides"));
  _ncell = parameters.get_value< CoordinateVector< int > >("densitygrid.ncell");
  HDF5Tools::close_group(group);

  // read cell midpoints, densities, and temperatures
  group = HDF5Tools::open_group(file, "/PartType0");
  std::vector< CoordinateVector<> > cell_midpoints =
      HDF5Tools::read_dataset< CoordinateVector<> >(group, "Coordinates");
  std::vector< double > cell_densities =
      HDF5Tools::read_dataset< double >(group, "NumberDensity");
  std::vector< double > cell_temperatures =
      HDF5Tools::read_dataset< double >(group, "Temperature");
  std::vector< std::vector< double > > neutral_fractions(NUMBER_OF_IONNAMES);
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    neutral_fractions[i] = HDF5Tools::read_dataset< double >(
        group, "NeutralFraction" + get_ion_name(i));
  }
  HDF5Tools::close_group(group);

  HDF5Tools::close_file(file);

  if (log) {
    log->write_status(
        "Constructing a CMacIonizeSnapshotDensityFunction containing ",
        _ncell.x(), " x ", _ncell.y(), " x ", _ncell.z(), " cells.");
  }

  _grid = new DensityValues **[_ncell.x()];
  for (int ix = 0; ix < _ncell.x(); ++ix) {
    _grid[ix] = new DensityValues *[_ncell.y()];
    for (int iy = 0; iy < _ncell.y(); ++iy) {
      _grid[ix][iy] = new DensityValues[_ncell.z()];
      for (int iz = 0; iz < _ncell.z(); ++iz) {
        _grid[ix][iy][iz].set_total_density(-1.);
      }
    }
  }

  for (unsigned int i = 0; i < cell_midpoints.size(); ++i) {
    CoordinateVector<> p = cell_midpoints[i];
    // the anchor of the box is always [0., 0., 0.] in the snapshot file
    int ix = _ncell.x() * p.x() / _box.get_sides().x();
    int iy = _ncell.y() * p.y() / _box.get_sides().y();
    int iz = _ncell.z() * p.z() / _box.get_sides().z();
    _grid[ix][iy][iz].set_total_density(cell_densities[i]);
    _grid[ix][iy][iz].set_temperature(cell_temperatures[i]);
    for (int j = 0; j < NUMBER_OF_IONNAMES; ++j) {
      IonName ion = static_cast< IonName >(j);
      _grid[ix][iy][iz].set_ionic_fraction(ion, neutral_fractions[j][i]);
    }
  }

  for (int ix = 0; ix < _ncell.x(); ++ix) {
    for (int iy = 0; iy < _ncell.y(); ++iy) {
      for (int iz = 0; iz < _ncell.z(); ++iz) {
        if (_grid[ix][iy][iz].get_total_density() < 0.) {
          cmac_error("No values found for cell (%i, %i, %i)!", ix, iy, iz);
        }
      }
    }
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
CMacIonizeSnapshotDensityFunction::CMacIonizeSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : CMacIonizeSnapshotDensityFunction(
          params.get_value< std::string >("densityfunction.filename"), log) {}

/**
 * @brief Destructor.
 *
 * Delete the internal density grid.
 */
CMacIonizeSnapshotDensityFunction::~CMacIonizeSnapshotDensityFunction() {
  for (int ix = 0; ix < _ncell.x(); ++ix) {
    for (int iy = 0; iy < _ncell.y(); ++iy) {
      delete[] _grid[ix][iy];
    }
    delete[] _grid[ix];
  }
  delete[] _grid;
}

/**
 * @brief Get the DensityValues at the given position.
 *
 * @param position CoordinateVector specifying a position (in m).
 * @return DensityValues at that position (in SI units).
 */
DensityValues CMacIonizeSnapshotDensityFunction::
operator()(CoordinateVector<> position) const {
  // get the indices of the cell containing the position
  int ix = _ncell.x() * (position.x() - _box.get_anchor().x()) /
           _box.get_sides().x();
  int iy = _ncell.y() * (position.y() - _box.get_anchor().y()) /
           _box.get_sides().y();
  int iz = _ncell.z() * (position.z() - _box.get_anchor().z()) /
           _box.get_sides().z();

  return _grid[ix][iy][iz];
}
