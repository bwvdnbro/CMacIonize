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
 * @file AsciiFileDensityFunction.cpp
 *
 * @brief AsciiFileDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AsciiFileDensityFunction.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include <cinttypes>
#include <fstream>
#include <sstream>

/**
 * @brief Constructor.
 *
 * @param filename Name of the ASCII text file to read.
 * @param ncell Dimensions of the grid.
 * @param box Box containing the grid (in m).
 * @param temperature Initial temperature of the ISM (in K).
 * @param length_unit_in_SI Length unit used in the ASCII file (in m).
 * @param density_unit_in_SI Density unit used in the ASCII file (in m^-3).
 * @param log Log to write logging info to.
 */
AsciiFileDensityFunction::AsciiFileDensityFunction(
    std::string filename, CoordinateVector< uint_fast32_t > ncell, Box<> box,
    double temperature, double length_unit_in_SI, double density_unit_in_SI,
    Log *log)
    : _ncell(ncell), _box(box), _temperature(temperature), _log(log) {

  _grid = new double **[_ncell.x()];
  for (uint_fast32_t i = 0; i < _ncell.x(); ++i) {
    _grid[i] = new double *[_ncell.y()];
    for (uint_fast32_t j = 0; j < _ncell.y(); ++j) {
      _grid[i][j] = new double[_ncell.z()];
      // initialize all cells with negative values
      // this way, we can check after reading if all cells got a file value
      for (uint_fast32_t k = 0; k < _ncell.z(); ++k) {
        _grid[i][j][k] = -1.;
      }
    }
  }

  std::ifstream file(filename);
  if (!file.is_open()) {
    cmac_error("Could not open file \"%s\"!", filename.c_str());
  }

  std::string line;
  while (getline(file, line)) {
    if (line[0] != '#') {
      double x, y, z, rho;
      std::stringstream linestream(line);
      linestream >> x >> y >> z >> rho;
      // unit conversion
      x *= length_unit_in_SI;
      y *= length_unit_in_SI;
      z *= length_unit_in_SI;
      rho *= density_unit_in_SI;
      // get the cell indices
      uint_fast32_t ix, iy, iz;
      ix = (x - _box.get_anchor().x()) / _box.get_sides().x() * _ncell.x();
      iy = (y - _box.get_anchor().y()) / _box.get_sides().y() * _ncell.y();
      iz = (z - _box.get_anchor().z()) / _box.get_sides().z() * _ncell.z();
      _grid[ix][iy][iz] = rho;
    }
  }

  // check that all cells received a value
  for (uint_fast32_t i = 0; i < _ncell.x(); ++i) {
    for (uint_fast32_t j = 0; j < _ncell.y(); ++j) {
      for (uint_fast32_t k = 0; k < _ncell.z(); ++k) {
        if (_grid[i][j][k] < 0.) {
          cmac_error("No value found for cell [%" PRIuFAST32 ", %" PRIuFAST32
                     ", %" PRIuFAST32 "]!",
                     i, j, k);
        }
      }
    }
  }

  if (_log) {
    _log->write_status("Successfully read in density grid from \"", filename,
                       "\".");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the ASCII file (required)
 *  - number of cells: Number of cells in the ASCII file (default: [64, 64, 64])
 *  - box anchor: Anchor of the box containing the cells (default: [-5. pc, -5.
 *    pc, -5. pc])
 *  - box sides: Side lengths of the box containing the cells (default: [10. pc,
 *    10. pc, 10. pc])
 *  - temperature: Initial temperature of the ISM (default: 8000. K)
 *  - length unit: Length unit (default: 1. m)
 *  - density unit: Density unit (default: 1. m^-3)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging information to.
 */
AsciiFileDensityFunction::AsciiFileDensityFunction(ParameterFile &params,
                                                   Log *log)
    : AsciiFileDensityFunction(
          params.get_value< std::string >("DensityFunction:filename"),
          params.get_value< CoordinateVector< uint_fast32_t > >(
              "DensityFunction:number of cells",
              CoordinateVector< uint_fast32_t >(64)),
          Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                    "DensityFunction:box anchor", "[-5. pc, -5. pc, -5. pc]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "DensityFunction:box sides", "[10. pc, 10. pc, 10. pc]")),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "DensityFunction:temperature", "8000. K"),
          params.get_physical_value< QUANTITY_LENGTH >(
              "DensityFunction:length unit", "1. m"),
          params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
              "DensityFunction:density unit", "1. m^-3"),
          log) {}

/**
 * @brief Destructor.
 *
 * Free memory used by the internal density grid.
 */
AsciiFileDensityFunction::~AsciiFileDensityFunction() {
  for (uint_fast32_t i = 0; i < _ncell.x(); ++i) {
    for (uint_fast32_t j = 0; j < _ncell.y(); ++j) {
      delete[] _grid[i][j];
    }
    delete[] _grid[i];
  }
  delete[] _grid;
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * We calculate in which cell of the grid the given coordinate sits, and return
 * the values in that cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues AsciiFileDensityFunction::operator()(const Cell &cell) const {
  DensityValues values;

  const CoordinateVector<> position = cell.get_cell_midpoint();

  uint_fast32_t ix, iy, iz;
  ix = (position.x() - _box.get_anchor().x()) / _box.get_sides().x() *
       _ncell.x();
  iy = (position.y() - _box.get_anchor().y()) / _box.get_sides().y() *
       _ncell.y();
  iz = (position.z() - _box.get_anchor().z()) / _box.get_sides().z() *
       _ncell.z();

  values.set_number_density(_grid[ix][iy][iz]);
  values.set_temperature(_temperature);
  values.set_ionic_fraction(ION_H_n, 1.e-6);
  values.set_ionic_fraction(ION_He_n, 1.e-6);
  return values;
}

/**
 * @brief Get the total number of hydrogen atoms in the file.
 *
 * @return Total number of hydrogen atoms.
 */
double AsciiFileDensityFunction::get_total_hydrogen_number() const {

  double side_x, side_y, side_z, cellvolume;
  side_x = _box.get_sides().x() / _ncell.x();
  side_y = _box.get_sides().y() / _ncell.y();
  side_z = _box.get_sides().z() / _ncell.z();
  cellvolume = side_x * side_y * side_z;
  double mtot = 0.;
  for (uint_fast32_t i = 0; i < _ncell.x(); ++i) {
    for (uint_fast32_t j = 0; j < _ncell.y(); ++j) {
      for (uint_fast32_t k = 0; k < _ncell.z(); ++k) {
        mtot += _grid[i][j][k] * cellvolume;
      }
    }
  }
  return mtot;
}
