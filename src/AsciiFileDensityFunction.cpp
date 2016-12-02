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
    std::string filename, CoordinateVector< int > ncell, Box box,
    double temperature, double length_unit_in_SI, double density_unit_in_SI,
    Log *log)
    : _ncell(ncell), _box(box), _temperature(temperature), _log(log) {
  _grid = new double **[_ncell.x()];
  for (int i = 0; i < _ncell.x(); ++i) {
    _grid[i] = new double *[_ncell.y()];
    for (int j = 0; j < _ncell.y(); ++j) {
      _grid[i][j] = new double[_ncell.z()];
      // initialize all cells with negative values
      // this way, we can check after reading if all cells got a file value
      for (int k = 0; k < _ncell.z(); ++k) {
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
      int ix, iy, iz;
      ix = (x - _box.get_anchor().x()) / _box.get_sides().x() * _ncell.x();
      iy = (y - _box.get_anchor().y()) / _box.get_sides().y() * _ncell.y();
      iz = (z - _box.get_anchor().z()) / _box.get_sides().z() * _ncell.z();
      _grid[ix][iy][iz] = rho;
    }
  }

  // check that all cells received a value
  for (int i = 0; i < _ncell.x(); ++i) {
    for (int j = 0; j < _ncell.y(); ++j) {
      for (int k = 0; k < _ncell.z(); ++k) {
        if (_grid[i][j][k] < 0.) {
          cmac_error("No value found for cell [%i, %i, %i]!", i, j, k);
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
 * @param params ParameterFile to read from.
 * @param log Log to write logging information to.
 */
AsciiFileDensityFunction::AsciiFileDensityFunction(ParameterFile &params,
                                                   Log *log)
    : AsciiFileDensityFunction(
          params.get_value< std::string >("densityfunction.filename"),
          params.get_value< CoordinateVector< int > >(
              "densityfunction.ncell", CoordinateVector< int >(64)),
          Box(params.get_physical_vector< QUANTITY_LENGTH >(
                  "densityfunction.box_anchor", "[0. m, 0. m, 0. m]"),
              params.get_physical_vector< QUANTITY_LENGTH >(
                  "densityfunction.box_sides", "[1. m, 1. m, 1. m]")),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "densityfunction.temperature", "8000. K"),
          params.get_physical_value< QUANTITY_LENGTH >(
              "densityfunction.length_unit", "1. m"),
          params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
              "densityfunction.density_unit", "1. m^-3"),
          log) {}

/**
 * @brief Destructor.
 *
 * Free memory used by internal density grid.
 */
AsciiFileDensityFunction::~AsciiFileDensityFunction() {
  for (int i = 0; i < _ncell.x(); ++i) {
    for (int j = 0; j < _ncell.y(); ++j) {
      delete[] _grid[i][j];
    }
    delete[] _grid[i];
  }
  delete[] _grid;
}

/**
 * @brief Function that gives the density for a given coordinate.
 *
 * We calculate in which cell of the grid the given coordinate sits, and return
 * the density in that cell.
 *
 * @param position CoordinateVector specifying a coordinate position (in m).
 * @return Density at the given coordinate (in m^-3).
 */
DensityValues AsciiFileDensityFunction::
operator()(CoordinateVector<> position) {
  DensityValues cell;

  int ix, iy, iz;
  ix = (position.x() - _box.get_anchor().x()) / _box.get_sides().x() *
       _ncell.x();
  iy = (position.y() - _box.get_anchor().y()) / _box.get_sides().y() *
       _ncell.y();
  iz = (position.z() - _box.get_anchor().z()) / _box.get_sides().z() *
       _ncell.z();

  cell.set_total_density(_grid[ix][iy][iz]);
  cell.set_temperature(_temperature);
  cell.set_ionic_fraction(ION_H_n, 1.e-6);
  cell.set_ionic_fraction(ION_He_n, 1.e-6);
  return cell;
}

/**
 * @brief Get the total number of hydrogen atoms in the file.
 *
 * @return Total number of hydrogen atoms.
 */
double AsciiFileDensityFunction::get_total_hydrogen_number() {
  double side_x, side_y, side_z, cellvolume;
  side_x = _box.get_sides().x() / _ncell.x();
  side_y = _box.get_sides().y() / _ncell.y();
  side_z = _box.get_sides().z() / _ncell.z();
  cellvolume = side_x * side_y * side_z;
  double mtot = 0.;
  for (int i = 0; i < _ncell.x(); ++i) {
    for (int j = 0; j < _ncell.y(); ++j) {
      for (int k = 0; k < _ncell.z(); ++k) {
        mtot += _grid[i][j][k] * cellvolume;
      }
    }
  }
  return mtot;
}
