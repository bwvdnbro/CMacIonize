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
    std::string filename, Log *log)
    : _cartesian_grid(nullptr), _amr_grid(nullptr) {
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
                 "densitygrid:box_anchor"),
             parameters.get_physical_vector< QUANTITY_LENGTH >(
                 "densitygrid:box_sides"));
  _ncell = parameters.get_value< CoordinateVector< int > >("densitygrid:ncell");
  std::string type = parameters.get_value< std::string >("densitygrid:type");
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

  if (type == "Cartesian") {
    _cartesian_grid = new DensityValues **[_ncell.x()];
    for (int ix = 0; ix < _ncell.x(); ++ix) {
      _cartesian_grid[ix] = new DensityValues *[_ncell.y()];
      for (int iy = 0; iy < _ncell.y(); ++iy) {
        _cartesian_grid[ix][iy] = new DensityValues[_ncell.z()];
        for (int iz = 0; iz < _ncell.z(); ++iz) {
          _cartesian_grid[ix][iy][iz].set_number_density(-1.);
        }
      }
    }

    for (unsigned int i = 0; i < cell_midpoints.size(); ++i) {
      CoordinateVector<> p = cell_midpoints[i];
      // the anchor of the box is always [0., 0., 0.] in the snapshot file
      int ix = _ncell.x() * p.x() / _box.get_sides().x();
      int iy = _ncell.y() * p.y() / _box.get_sides().y();
      int iz = _ncell.z() * p.z() / _box.get_sides().z();
      _cartesian_grid[ix][iy][iz].set_number_density(cell_densities[i]);
      _cartesian_grid[ix][iy][iz].set_temperature(cell_temperatures[i]);
      for (int j = 0; j < NUMBER_OF_IONNAMES; ++j) {
        IonName ion = static_cast< IonName >(j);
        _cartesian_grid[ix][iy][iz].set_ionic_fraction(ion,
                                                       neutral_fractions[j][i]);
      }
    }

    for (int ix = 0; ix < _ncell.x(); ++ix) {
      for (int iy = 0; iy < _ncell.y(); ++iy) {
        for (int iz = 0; iz < _ncell.z(); ++iz) {
          if (_cartesian_grid[ix][iy][iz].get_number_density() < 0.) {
            cmac_error("No values found for cell (%i, %i, %i)!", ix, iy, iz);
          }
        }
      }
    }
  } else if (type == "AMR") {
    // find the smallest number of blocks that fits the requested top level grid
    int power_of_2_x = get_power_of_two(_ncell.x());
    int power_of_2_y = get_power_of_two(_ncell.y());
    int power_of_2_z = get_power_of_two(_ncell.z());
    int power_of_2 = std::min(power_of_2_x, power_of_2_y);
    power_of_2 = std::min(power_of_2, power_of_2_z);
    CoordinateVector< int > nblock = _ncell / power_of_2;
    // find out how many cells each block should have at the lowest level
    // this is just the power in power_of_2
    unsigned char level = 0;
    while (power_of_2 > 1) {
      power_of_2 >>= 1;
      ++level;
    }
    // this creates all cells on the highest level
    _amr_grid = new AMRGrid< DensityValues >(_box, nblock);
    _amr_grid->create_all_cells(level);
    // now try to fill the cells
    for (unsigned int i = 0; i < cell_midpoints.size(); ++i) {
      CoordinateVector<> p = cell_midpoints[i];
      p += _box.get_anchor();
      unsigned long key = _amr_grid->get_key(p);
      AMRGridCell< DensityValues > *cell = &(*_amr_grid)[key];
      CoordinateVector<> cell_midpoint = cell->get_midpoint();
      // now check if the cell midpoint equals p. If not, we have to refine the
      // cell.
      while (std::abs(cell_midpoint.x() - p.x()) >
                 1.e-9 * std::abs(cell_midpoint.x() + p.x()) &&
             std::abs(cell_midpoint.y() - p.y()) >
                 1.e-9 * std::abs(cell_midpoint.y() + p.y()) &&
             std::abs(cell_midpoint.z() - p.z()) >
                 1.e-9 * std::abs(cell_midpoint.z() + p.z())) {
        cell->create_all_cells(cell->get_level(), cell->get_level() + 1);
        cell = cell->get_child(p);
        cell_midpoint = cell->get_midpoint();
      }
      // cell now is the actual cell corresponding to the midpoint
      // set its contents
      DensityValues &values = cell->value();
      values.set_number_density(cell_densities[i]);
      values.set_temperature(cell_temperatures[i]);
      for (int j = 0; j < NUMBER_OF_IONNAMES; ++j) {
        IonName ion = static_cast< IonName >(j);
        values.set_ionic_fraction(ion, neutral_fractions[j][i]);
      }
    }
  } else {
    cmac_error("Reconstructing a density field from a %sDensityGrid is not yet "
               "supported!",
               type.c_str());
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
          params.get_value< std::string >("densityfunction:filename"), log) {}

/**
 * @brief Destructor.
 *
 * Delete the internal density grid.
 */
CMacIonizeSnapshotDensityFunction::~CMacIonizeSnapshotDensityFunction() {
  if (_cartesian_grid) {
    for (int ix = 0; ix < _ncell.x(); ++ix) {
      for (int iy = 0; iy < _ncell.y(); ++iy) {
        delete[] _cartesian_grid[ix][iy];
      }
      delete[] _cartesian_grid[ix];
    }
    delete[] _cartesian_grid;
  } else if (_amr_grid) {
    delete _amr_grid;
  }
}

/**
 * @brief Get the DensityValues at the given position.
 *
 * @param position CoordinateVector specifying a position (in m).
 * @return DensityValues at that position (in SI units).
 */
DensityValues CMacIonizeSnapshotDensityFunction::
operator()(CoordinateVector<> position) const {
  if (_cartesian_grid) {
    // get the indices of the cell containing the position
    int ix = _ncell.x() * (position.x() - _box.get_anchor().x()) /
             _box.get_sides().x();
    int iy = _ncell.y() * (position.y() - _box.get_anchor().y()) /
             _box.get_sides().y();
    int iz = _ncell.z() * (position.z() - _box.get_anchor().z()) /
             _box.get_sides().z();

    return _cartesian_grid[ix][iy][iz];
  } else if (_amr_grid) {
    unsigned long key = _amr_grid->get_key(position);
    AMRGridCell< DensityValues > &cell = (*_amr_grid)[key];
    DensityValues values = cell.value();
    // we make sure cells are refined to the correct level
    CoordinateVector<> cell_midpoint = cell.get_midpoint();
    if (std::abs(cell_midpoint.x() - position.x()) >
            1.e-9 * std::abs(cell_midpoint.x() + position.x()) &&
        std::abs(cell_midpoint.y() - position.y()) >
            1.e-9 * std::abs(cell_midpoint.y() + position.y()) &&
        std::abs(cell_midpoint.z() - position.z()) >
            1.e-9 * std::abs(cell_midpoint.z() + position.z())) {
      values.set_number_density(-1.);
    }
    return values;
  } else {
    cmac_error("This grid type is not supported (and you should never see this "
               "error)!");
    return DensityValues();
  }
}
