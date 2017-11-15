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
#include <cinttypes>

/**
 * @brief Constructor.
 *
 * @param filename Name of the snapshot file to read.
 * @param log Log to write logging info to.
 */
CMacIonizeSnapshotDensityFunction::CMacIonizeSnapshotDensityFunction(
    std::string filename, Log *log)
    : _cartesian_grid(nullptr), _amr_grid(nullptr),
      _voronoi_pointlocations(nullptr) {

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
  _box = Box<>(
      parameters.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:anchor"),
      parameters.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:sides"));
  _ncell = parameters.get_value< CoordinateVector< uint_fast32_t > >(
      "DensityGrid:number of cells", CoordinateVector< uint_fast32_t >(-1));
  std::string type = parameters.get_value< std::string >("DensityGrid:type");
  HDF5Tools::close_group(group);

  // units
  double unit_length_in_SI = 1.;
  double unit_density_in_SI = 1.;
  double unit_temperature_in_SI = 1.;
  if (HDF5Tools::group_exists(file, "/Units")) {
    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Units");
    double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit length in cgs (U_L)");
    double unit_temperature_in_cgs = HDF5Tools::read_attribute< double >(
        units, "Unit temperature in cgs (U_T)");
    unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    unit_density_in_SI =
        1. / unit_length_in_SI / unit_length_in_SI / unit_length_in_SI;
    // K is K
    unit_temperature_in_SI = unit_temperature_in_cgs;
    HDF5Tools::close_group(units);
  }

  // read cell midpoints, densities, and temperatures
  group = HDF5Tools::open_group(file, "/PartType0");
  std::vector< CoordinateVector<> > cell_midpoints =
      HDF5Tools::read_dataset< CoordinateVector<> >(group, "Coordinates");
  std::vector< double > cell_densities =
      HDF5Tools::read_dataset< double >(group, "NumberDensity");
  std::vector< double > cell_temperatures =
      HDF5Tools::read_dataset< double >(group, "Temperature");
  std::vector< std::vector< double > > neutral_fractions(NUMBER_OF_IONNAMES);
  for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    neutral_fractions[i] = HDF5Tools::read_dataset< double >(
        group, "NeutralFraction" + get_ion_name(i));
  }
  HDF5Tools::close_group(group);

  HDF5Tools::close_file(file);

  // unit conversion
  for (size_t i = 0; i < cell_midpoints.size(); ++i) {
    cell_midpoints[i][0] *= unit_length_in_SI;
    cell_midpoints[i][1] *= unit_length_in_SI;
    cell_midpoints[i][2] *= unit_length_in_SI;
    cell_densities[i] *= unit_density_in_SI;
    cell_temperatures[i] *= unit_temperature_in_SI;
  }

  if (log) {
    log->write_status(
        "Constructing a CMacIonizeSnapshotDensityFunction containing ",
        _ncell.x(), " x ", _ncell.y(), " x ", _ncell.z(), " cells.");
  }

  if (type == "Cartesian") {
    _cartesian_grid = new DensityValues **[_ncell.x()];
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      _cartesian_grid[ix] = new DensityValues *[_ncell.y()];
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        _cartesian_grid[ix][iy] = new DensityValues[_ncell.z()];
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          _cartesian_grid[ix][iy][iz].set_number_density(-1.);
        }
      }
    }

    for (size_t i = 0; i < cell_midpoints.size(); ++i) {
      CoordinateVector<> p = cell_midpoints[i];
      // the anchor of the box is always [0., 0., 0.] in the snapshot file
      uint_fast32_t ix = _ncell.x() * p.x() / _box.get_sides().x();
      uint_fast32_t iy = _ncell.y() * p.y() / _box.get_sides().y();
      uint_fast32_t iz = _ncell.z() * p.z() / _box.get_sides().z();
      _cartesian_grid[ix][iy][iz].set_number_density(cell_densities[i]);
      _cartesian_grid[ix][iy][iz].set_temperature(cell_temperatures[i]);
      for (int_fast32_t j = 0; j < NUMBER_OF_IONNAMES; ++j) {
        IonName ion = static_cast< IonName >(j);
        _cartesian_grid[ix][iy][iz].set_ionic_fraction(ion,
                                                       neutral_fractions[j][i]);
      }
    }

    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          if (_cartesian_grid[ix][iy][iz].get_number_density() < 0.) {
            cmac_error("No values found for cell (%" PRIuFAST32 ", %" PRIuFAST32
                       ", %" PRIuFAST32 ")!",
                       ix, iy, iz);
          }
        }
      }
    }
  } else if (type == "AMR") {
    // find the smallest number of blocks that fits the requested top level grid
    uint_fast32_t power_of_2_x = get_power_of_two(_ncell.x());
    uint_fast32_t power_of_2_y = get_power_of_two(_ncell.y());
    uint_fast32_t power_of_2_z = get_power_of_two(_ncell.z());
    uint_fast32_t power_of_2 = std::min(power_of_2_x, power_of_2_y);
    power_of_2 = std::min(power_of_2, power_of_2_z);
    CoordinateVector< uint_fast32_t > nblock = _ncell / power_of_2;
    // find out how many cells each block should have at the lowest level
    // this is just the power in power_of_2
    uint_fast8_t level = 0;
    while (power_of_2 > 1) {
      power_of_2 >>= 1;
      ++level;
    }
    // this creates all cells on the highest level
    _amr_grid = new AMRGrid< DensityValues >(_box, nblock);
    _amr_grid->create_all_cells(level);
    // now try to fill the cells
    for (size_t i = 0; i < cell_midpoints.size(); ++i) {
      CoordinateVector<> p = cell_midpoints[i];
      p += _box.get_anchor();
      amrkey_t key = _amr_grid->get_key(p);
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
      for (int_fast32_t j = 0; j < NUMBER_OF_IONNAMES; ++j) {
        IonName ion = static_cast< IonName >(j);
        values.set_ionic_fraction(ion, neutral_fractions[j][i]);
      }
    }
  } else if (type == "Voronoi") {
    _voronoi_generators.resize(cell_midpoints.size());
    _voronoi_densityvalues.resize(cell_midpoints.size());
    for (size_t i = 0; i < cell_midpoints.size(); ++i) {
      _voronoi_generators[i] = cell_midpoints[i] + _box.get_anchor();
      _voronoi_densityvalues[i].set_number_density(cell_densities[i]);
      _voronoi_densityvalues[i].set_temperature(cell_temperatures[i]);
      for (int_fast32_t j = 0; j < NUMBER_OF_IONNAMES; ++j) {
        IonName ion = static_cast< IonName >(j);
        _voronoi_densityvalues[i].set_ionic_fraction(ion,
                                                     neutral_fractions[j][i]);
      }
    }
    _voronoi_pointlocations =
        new PointLocations(_voronoi_generators, 100, _box);
  } else {
    cmac_error("Reconstructing a density field from a %sDensityGrid is not yet "
               "supported!",
               type.c_str());
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the snapshot file to read (required)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
CMacIonizeSnapshotDensityFunction::CMacIonizeSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : CMacIonizeSnapshotDensityFunction(
          params.get_value< std::string >("DensityFunction:filename"), log) {}

/**
 * @brief Destructor.
 *
 * Delete the internal density grid.
 */
CMacIonizeSnapshotDensityFunction::~CMacIonizeSnapshotDensityFunction() {
  if (_cartesian_grid) {
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        delete[] _cartesian_grid[ix][iy];
      }
      delete[] _cartesian_grid[ix];
    }
    delete[] _cartesian_grid;
  } else if (_amr_grid) {
    delete _amr_grid;
  } else if (_voronoi_pointlocations) {
    delete _voronoi_pointlocations;
  }
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues CMacIonizeSnapshotDensityFunction::
operator()(const Cell &cell) const {

  const CoordinateVector<> position = cell.get_cell_midpoint();

  if (_cartesian_grid) {
    // get the indices of the cell containing the position
    const uint_fast32_t ix = _ncell.x() *
                             (position.x() - _box.get_anchor().x()) /
                             _box.get_sides().x();
    const uint_fast32_t iy = _ncell.y() *
                             (position.y() - _box.get_anchor().y()) /
                             _box.get_sides().y();
    const uint_fast32_t iz = _ncell.z() *
                             (position.z() - _box.get_anchor().z()) /
                             _box.get_sides().z();

    cmac_assert_message(ix < _ncell.x(), "%" PRIuFAST32, ix);
    cmac_assert_message(iy < _ncell.y(), "%" PRIuFAST32, iy);
    cmac_assert_message(iz < _ncell.z(), "%" PRIuFAST32, iz);

    return _cartesian_grid[ix][iy][iz];
  } else if (_amr_grid) {
    const amrkey_t key = _amr_grid->get_key(position);
    const AMRGridCell< DensityValues > &amr_cell = (*_amr_grid)[key];
    DensityValues values = amr_cell.value();
    // we make sure cells are refined to the correct level
    const CoordinateVector<> cell_midpoint = amr_cell.get_midpoint();
    if (std::abs(cell_midpoint.x() - position.x()) >
            1.e-9 * std::abs(cell_midpoint.x() + position.x()) &&
        std::abs(cell_midpoint.y() - position.y()) >
            1.e-9 * std::abs(cell_midpoint.y() + position.y()) &&
        std::abs(cell_midpoint.z() - position.z()) >
            1.e-9 * std::abs(cell_midpoint.z() + position.z())) {
      values.set_number_density(-1.);
    }
    return values;
  } else if (_voronoi_pointlocations) {
    const uint_fast32_t index =
        _voronoi_pointlocations->get_closest_neighbour(position);
    return _voronoi_densityvalues[index];
  } else {
    cmac_error("This grid type is not supported (and you should never see this "
               "error)!");
    return DensityValues();
  }
}
