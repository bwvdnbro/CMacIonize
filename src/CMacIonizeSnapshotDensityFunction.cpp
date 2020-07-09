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
 * @param use_density Use the mass density instead of the number density?
 * @param use_pressure Use the pressure instead of the temperature?
 * @param initial_neutral_fraction Initial neutral fraction to use (if
 * neutral fractions are not present in the file).
 * @param log Log to write logging info to.
 */
CMacIonizeSnapshotDensityFunction::CMacIonizeSnapshotDensityFunction(
    std::string filename, const bool use_density, const bool use_pressure,
    const double initial_neutral_fraction, Log *log)
    : _filename(filename), _use_density(use_density),
      _use_pressure(use_pressure),
      _initial_neutral_fraction(initial_neutral_fraction), _log(log),
      _cartesian_grid(nullptr), _amr_grid(nullptr),
      _voronoi_pointlocations(nullptr) {

  std::ifstream file(filename);
  if (!file.is_open()) {
    cmac_error("Could not open file \"%s\"!", filename.c_str());
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the snapshot file to read (required)
 *  - use density: Use the mass density instead of the number density? (default:
 *    no)
 *  - use pressure: Use the pressure instead of the temperature? (default: no)
 *  - initial neutral fraction: Initial value for the neutral fractions if they
 *    are not present in the file (default: 1.e-6)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
CMacIonizeSnapshotDensityFunction::CMacIonizeSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : CMacIonizeSnapshotDensityFunction(
          params.get_filename("DensityFunction:filename"),
          params.get_value< bool >("DensityFunction:use density", false),
          params.get_value< bool >("DensityFunction:use pressure", false),
          params.get_value< double >("DensityFunction:initial neutral fraction",
                                     1.e-6),
          log) {}

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
 * @brief Perform all computationally expensive initialization that needs to
 * be done before operator() will work.
 *
 * This routine actually reads the file.
 */
void CMacIonizeSnapshotDensityFunction::initialize() {

  if (_log) {
    _log->write_info("Opening file ", _filename, "...");
  }

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(_filename, HDF5Tools::HDF5FILEMODE_READ);

  if (_log) {
    _log->write_info("Done opening file.");
    _log->write_info("Reading parameter block...");
  }

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
  std::string type;
  if (parameters.has_value("DensityGrid:type")) {
    type = parameters.get_value< std::string >("DensityGrid:type");
  } else {
    type = "TaskBased";
  }
  CoordinateVector< uint_fast32_t > numsubgrid;
  if (type == "TaskBased") {
    numsubgrid = parameters.get_value< CoordinateVector< uint_fast32_t > >(
        "DensitySubGridCreator:number of subgrids");
  }

  HDF5Tools::close_group(group);

  if (_log) {
    _log->write_info("Done reading parameters.");
    _log->write_info("Reading unit block...");
  }

  // units
  double unit_length_in_SI = 1.;
  double unit_density_in_SI = 1.;
  double unit_temperature_in_SI = 1.;
  double unit_velocity_in_SI = 1.;
  if (HDF5Tools::group_exists(file, "/Units")) {
    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Units");
    const double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit length in cgs (U_L)");
    const double unit_temperature_in_cgs = HDF5Tools::read_attribute< double >(
        units, "Unit temperature in cgs (U_T)");
    const double unit_time_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit time in cgs (U_t)");
    unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    unit_density_in_SI =
        1. / unit_length_in_SI / unit_length_in_SI / unit_length_in_SI;
    // K is K
    unit_temperature_in_SI = unit_temperature_in_cgs;
    // seconds are seconds
    unit_velocity_in_SI = unit_length_in_SI / unit_time_in_cgs;
    HDF5Tools::close_group(units);
  }

  if (_log) {
    _log->write_info("Done reading units.");
    _log->write_info("Reading particle data...");
  }

  // read cell midpoints, number densities and temperatures
  group = HDF5Tools::open_group(file, "/PartType0");

  if (_log) {
    _log->write_info("Coordinates...");
  }

  std::vector< CoordinateVector<> > cell_midpoints;
  if (type != "TaskBased") {
    cell_midpoints =
        HDF5Tools::read_dataset< CoordinateVector<> >(group, "Coordinates");
  }

  if (_log) {
    _log->write_info("Densities...");
  }

  std::vector< double > cell_densities;
  if (HDF5Tools::group_exists(group, "NumberDensity") && !_use_density) {
    cell_densities = HDF5Tools::read_dataset< double >(group, "NumberDensity");
  } else {
    cell_densities = HDF5Tools::read_dataset< double >(group, "Density");
    unit_density_in_SI /=
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
  }

  if (_log) {
    _log->write_info("Ionic fractions...");
  }

  std::vector< std::vector< double > > neutral_fractions(
      NUMBER_OF_IONNAMES,
      std::vector< double >(cell_densities.size(), _initial_neutral_fraction));
  for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    // skip ionic fractions that do not exist
    if (HDF5Tools::group_exists(group, "NeutralFraction" + get_ion_name(i))) {
      neutral_fractions[i] = HDF5Tools::read_dataset< double >(
          group, "NeutralFraction" + get_ion_name(i));
    }
  }

  if (_log) {
    _log->write_info("Temperatures...");
  }

  std::vector< double > cell_temperatures;
  if (HDF5Tools::group_exists(group, "Temperature") && !_use_pressure) {
    cell_temperatures = HDF5Tools::read_dataset< double >(group, "Temperature");
  } else {
    const double kB =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    cell_temperatures = HDF5Tools::read_dataset< double >(group, "Pressure");
    for (size_t i = 0; i < cell_temperatures.size(); ++i) {
      const double mu = 0.5 * (1. + neutral_fractions[ION_H_n][i]);
      cell_temperatures[i] *= mu / (cell_densities[i] * unit_density_in_SI *
                                    kB * unit_temperature_in_SI);
    }
  }

  if (_log) {
    _log->write_info("Velocities...");
  }

  // velocities (if they exist)
  std::vector< CoordinateVector<> > cell_velocities;
  if (HDF5Tools::group_exists(group, "Velocities")) {
    cell_velocities =
        HDF5Tools::read_dataset< CoordinateVector<> >(group, "Velocities");
  }

  HDF5Tools::close_group(group);

  if (_log) {
    _log->write_info("Done reading particle data.");
    _log->write_info("Closing file.");
  }

  HDF5Tools::close_file(file);

  if (_log) {
    _log->write_info("Converting units...");
  }

  // unit conversion
  for (size_t i = 0; i < cell_densities.size(); ++i) {
    if (cell_midpoints.size() > 0) {
      cell_midpoints[i][0] *= unit_length_in_SI;
      cell_midpoints[i][1] *= unit_length_in_SI;
      cell_midpoints[i][2] *= unit_length_in_SI;
    }
    cell_densities[i] *= unit_density_in_SI;
    cell_temperatures[i] *= unit_temperature_in_SI;
    if (cell_velocities.size() > 0) {
      cell_velocities[i][0] *= unit_velocity_in_SI;
      cell_velocities[i][1] *= unit_velocity_in_SI;
      cell_velocities[i][2] *= unit_velocity_in_SI;
    }
  }

  if (_log) {
    _log->write_info("Done converting units.");
    _log->write_status(
        "Constructing a CMacIonizeSnapshotDensityFunction containing ",
        _ncell.x(), " x ", _ncell.y(), " x ", _ncell.z(), " cells.");
  }

  if (_log) {
    _log->write_info("Creating grid structure...");
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
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        _cartesian_grid[ix][iy][iz].set_ionic_fraction(
            ion, neutral_fractions[ion][i]);
      }
      if (cell_velocities.size() > 0) {
        _cartesian_grid[ix][iy][iz].set_velocity(cell_velocities[i]);
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
  } else if (type == "TaskBased") {
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

    const uint_fast32_t numblockx = _ncell.x() / numsubgrid.x();
    const uint_fast32_t numblocky = _ncell.y() / numsubgrid.y();
    const uint_fast32_t numblockz = _ncell.z() / numsubgrid.z();
    const uint_fast32_t numblocktot = numblockx * numblocky * numblockz;

    for (uint_fast32_t six = 0; six < numsubgrid.x(); ++six) {
      for (uint_fast32_t siy = 0; siy < numsubgrid.y(); ++siy) {
        for (uint_fast32_t siz = 0; siz < numsubgrid.z(); ++siz) {
          const uint_fast32_t subgrid_index =
              six * numsubgrid.y() * numsubgrid.z() + siy * numsubgrid.z() +
              siz;
          for (uint_fast32_t cix = 0; cix < numblockx; ++cix) {
            for (uint_fast32_t ciy = 0; ciy < numblocky; ++ciy) {
              for (uint_fast32_t ciz = 0; ciz < numblockz; ++ciz) {
                const uint_fast32_t cell_index = subgrid_index * numblocktot +
                                                 cix * numblocky * numblockz +
                                                 ciy * numblockz + ciz;
                const uint_fast32_t ix = six * numblockx + cix;
                const uint_fast32_t iy = siy * numblocky + ciy;
                const uint_fast32_t iz = siz * numblockz + ciz;
                _cartesian_grid[ix][iy][iz].set_number_density(
                    cell_densities[cell_index]);
                _cartesian_grid[ix][iy][iz].set_temperature(
                    cell_temperatures[cell_index]);
                for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                  _cartesian_grid[ix][iy][iz].set_ionic_fraction(
                      ion, neutral_fractions[ion][cell_index]);
                }
                if (cell_velocities.size() > 0) {
                  _cartesian_grid[ix][iy][iz].set_velocity(
                      cell_velocities[cell_index]);
                }
              }
            }
          }
        }
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
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        values.set_ionic_fraction(ion, neutral_fractions[ion][i]);
      }
      if (cell_velocities.size() > 0) {
        values.set_velocity(cell_velocities[i]);
      }
    }
  } else if (type == "Voronoi") {
    _voronoi_generators.resize(cell_midpoints.size());
    _voronoi_densityvalues.resize(cell_midpoints.size());
    for (size_t i = 0; i < cell_midpoints.size(); ++i) {
      _voronoi_generators[i] = cell_midpoints[i] + _box.get_anchor();
      _voronoi_densityvalues[i].set_number_density(cell_densities[i]);
      _voronoi_densityvalues[i].set_temperature(cell_temperatures[i]);
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        _voronoi_densityvalues[i].set_ionic_fraction(ion,
                                                     neutral_fractions[ion][i]);
      }
      if (cell_velocities.size() > 0) {
        _voronoi_densityvalues[i].set_velocity(cell_velocities[i]);
      }
    }
    _voronoi_pointlocations =
        new PointLocations(_voronoi_generators, 100, _box);
  } else {
    cmac_error("Reconstructing a density field from a %sDensityGrid is not yet "
               "supported!",
               type.c_str());
  }

  if (_log) {
    _log->write_info("Done creating grid structure.");
  }
}

/**
 * @brief Free up the memory used by the density function. After this,
 * operator() will no longer work.
 */
void CMacIonizeSnapshotDensityFunction::free() {
  if (_cartesian_grid) {
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        delete[] _cartesian_grid[ix][iy];
      }
      delete[] _cartesian_grid[ix];
    }
    delete[] _cartesian_grid;
    _cartesian_grid = nullptr;
  } else if (_amr_grid) {
    delete _amr_grid;
    _amr_grid = nullptr;
  } else if (_voronoi_pointlocations) {
    delete _voronoi_pointlocations;
    _voronoi_pointlocations = nullptr;
  }
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues CMacIonizeSnapshotDensityFunction::operator()(const Cell &cell) {

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
