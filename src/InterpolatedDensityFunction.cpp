/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file InterpolatedDensityFunction.cpp
 *
 * @brief InterpolatedDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "InterpolatedDensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Utilities.hpp"
#include "YAMLDictionary.hpp"
#include <fstream>
#include <iostream>

/**
 * @brief Constructor.
 *
 * Open the file and parse its contents. We always construct ranges for the
 * three coordinate directions; if no values are present for a specific range in
 * the file, we still create a two value range, spanning from the minimal
 * coordinate value in that direction to the maximal coordinate value. The
 * linear interpolation will then automatically keep the number density constant
 * in that coordinate direction.
 *
 * If the coordinate ranges that are provided do not include the minimal and
 * maximal value in that direction, the linear interpolation will automatically
 * linearly extrapolate for points within the valid range, but outside of the
 * provided range. We do not add extra points to the range to achieve this.
 *
 * @param filename Name of the file containing the data.
 * @param temperature Initial temperature (in K).
 * @param log Log to write logging info to.
 */
InterpolatedDensityFunction::InterpolatedDensityFunction(std::string filename,
                                                         double temperature,
                                                         Log *log)
    : _temperature(temperature) {

  std::ifstream file(filename);
  if (!file) {
    cmac_error("Error while opening file \"%s\"!", filename.c_str());
  }

  if (log) {
    log->write_info("Reading density function from file ", filename, "...");
  }

  // extract the YAML block from the head of the file
  std::string line;
  while (std::getline(file, line) && line != "---") {
  }
  if (line != "---") {
    cmac_error("No YAML block found in file \"%s\"!", filename.c_str());
  }
  // ok, we found the starting line
  // add subsequent lines to a std::string, and look for the ending line
  std::string yaml_block;
  while (std::getline(file, line) && line != "---") {
    yaml_block += line + "\n";
  }
  if (line != "---") {
    cmac_error("Reached end of file \"%s\" while parsing YAML block!",
               filename.c_str());
  }
  std::istringstream yaml_stream(yaml_block);
  YAMLDictionary yaml_dictionary(yaml_stream);

  // get the number of x, y and z points
  unsigned int num_x = yaml_dictionary.get_value< unsigned int >("num_x");
  unsigned int num_y = yaml_dictionary.get_value< unsigned int >("num_y");
  unsigned int num_z = yaml_dictionary.get_value< unsigned int >("num_z");
  // get the extents of the box in x, y and z
  _x_bounds.first =
      yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("xmin");
  _x_bounds.second =
      yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("xmax");
  _y_bounds.first =
      yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("ymin");
  _y_bounds.second =
      yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("ymax");
  _z_bounds.first =
      yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("zmin");
  _z_bounds.second =
      yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("zmax");
  // get the number of columns
  unsigned int num_column =
      yaml_dictionary.get_value< unsigned int >("num_column");
  // get the column names and units
  std::map< std::string, unsigned int > name_to_column;
  std::vector< std::string > units(num_column);
  for (unsigned int i = 0; i < num_column; ++i) {
    std::stringstream column_name_stream;
    column_name_stream << "column_" << i << "_";
    std::string column_name = column_name_stream.str();
    std::string name =
        yaml_dictionary.get_value< std::string >(column_name + "variable");
    std::string unit =
        yaml_dictionary.get_value< std::string >(column_name + "unit");
    name_to_column[name] = i;
    units[i] = unit;
  }

  // now check if everything works out and set up the coordinate vectors
  if (num_x == 0 && num_y == 0 && num_z == 0) {
    cmac_error("No coordinate values provided! We need at least one "
               "non-trivial coordinate axis.");
  }
  if (_x_bounds.first > _x_bounds.second) {
    cmac_error("Minimal X value larger than maximal X value!");
  }
  if (_y_bounds.first > _y_bounds.second) {
    cmac_error("Minimal Y value larger than maximal Y value!");
  }
  if (_z_bounds.first > _z_bounds.second) {
    cmac_error("Minimal Z value larger than maximal Z value!");
  }
  unsigned int x_column = 0;
  unsigned int y_column = 0;
  unsigned int z_column = 0;
  unsigned int num_variable = 1;
  if (num_x != 0) {
    if (name_to_column.count("x") == 0) {
      cmac_error("No column found containing x values!");
    }
    x_column = name_to_column["x"];
  }
  if (num_x > 1) {
    _x_coords.resize(num_x);
    _number_densities.resize(num_x);
    num_variable *= num_x;
  } else {
    _x_coords.resize(2);
    _x_coords[0] = _x_bounds.first;
    _x_coords[1] = _x_bounds.second;
    _number_densities.resize(2);
  }
  if (num_y != 0) {
    if (name_to_column.count("y") == 0) {
      cmac_error("No column found containing y values!");
    }
    y_column = name_to_column["y"];
  }
  if (num_y > 1) {
    _y_coords.resize(num_y);
    for (unsigned int ix = 0; ix < _number_densities.size(); ++ix) {
      _number_densities[ix].resize(num_y);
    }
    num_variable *= num_y;
  } else {
    _y_coords.resize(2);
    _y_coords[0] = _y_bounds.first;
    _y_coords[1] = _y_bounds.second;
    for (unsigned int ix = 0; ix < _number_densities.size(); ++ix) {
      _number_densities[ix].resize(2);
    }
  }
  if (num_z != 0) {
    if (name_to_column.count("z") == 0) {
      cmac_error("No column found containing z values!");
    }
    z_column = name_to_column["z"];
  }
  if (num_z > 1) {
    _z_coords.resize(num_z);
    for (unsigned int ix = 0; ix < _number_densities.size(); ++ix) {
      for (unsigned int iy = 0; iy < _number_densities[ix].size(); ++iy) {
        _number_densities[ix][iy].resize(num_z);
      }
    }
    num_variable *= num_z;
  } else {
    _z_coords.resize(2);
    _z_coords[0] = _z_bounds.first;
    _z_coords[1] = _z_bounds.second;
    for (unsigned int ix = 0; ix < _number_densities.size(); ++ix) {
      for (unsigned int iy = 0; iy < _number_densities[ix].size(); ++iy) {
        _number_densities[ix][iy].resize(2);
      }
    }
  }
  if (name_to_column.count("number density") == 0) {
    cmac_error("No column found containing number density values!");
  }
  unsigned int number_density_column = name_to_column["number density"];

  unsigned int ix = 0;
  unsigned int iy = 0;
  unsigned int iz = 0;
  unsigned int i = 0;
  while (std::getline(file, line)) {
    std::stringstream lstream(line);
    std::vector< double > row(num_column);
    for (unsigned int j = 0; j < num_column; ++j) {
      lstream >> row[j];
    }

    if (num_x != 0) {
      double next_x = UnitConverter::to_SI< QUANTITY_LENGTH >(row[x_column],
                                                              units[x_column]);
      cmac_assert(next_x >= _x_bounds.first && next_x <= _x_bounds.second);
      if (i > 0 && next_x != _x_coords[ix]) {
        ++ix;
        cmac_assert(ix < _x_coords.size());
      }
      _x_coords[ix] = next_x;
    }
    if (num_y != 0) {
      double next_y = UnitConverter::to_SI< QUANTITY_LENGTH >(row[y_column],
                                                              units[y_column]);
      cmac_assert(next_y >= _y_bounds.first && next_y <= _y_bounds.second);
      if (i > 0 && next_y != _y_coords[iy]) {
        ++iy;
        cmac_assert(iy < _y_coords.size());
      }
      _y_coords[iy] = next_y;
    }
    if (num_z != 0) {
      double next_z = UnitConverter::to_SI< QUANTITY_LENGTH >(row[z_column],
                                                              units[z_column]);
      cmac_assert(next_z >= _z_bounds.first && next_z <= _z_bounds.second);
      if (i > 0 && next_z != _z_coords[iz]) {
        ++iz;
        cmac_assert(iz < _z_coords.size());
      }
      _z_coords[iz] = next_z;
    }
    _number_densities[ix][iy][iz] =
        UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(
            row[number_density_column], units[number_density_column]);

    ++i;
    cmac_assert(i <= num_variable);
  }

  cmac_assert(i == num_variable);

  // now copy the arrays to complete the missing data
  if (num_x < 2) {
    cmac_assert(ix == 0);
    for (unsigned int iy = 0; iy < _y_coords.size(); ++iy) {
      for (unsigned int iz = 0; iz < _z_coords.size(); ++iz) {
        _number_densities[1][iy][iz] = _number_densities[0][iy][iz];
      }
    }
  } else {
    cmac_assert(ix == _x_coords.size() - 1);
  }
  if (num_y < 2) {
    cmac_assert(iy == 0);
    for (unsigned int ix = 0; ix < _x_coords.size(); ++ix) {
      for (unsigned int iz = 0; iz < _z_coords.size(); ++iz) {
        _number_densities[ix][1][iz] = _number_densities[ix][0][iz];
      }
    }
  } else {
    cmac_assert(iy == _y_coords.size() - 1);
  }
  if (num_z < 2) {
    cmac_assert(iz == 0);
    for (unsigned int ix = 0; ix < _x_coords.size(); ++ix) {
      for (unsigned int iy = 0; iy < _y_coords.size(); ++iy) {
        _number_densities[ix][iy][1] = _number_densities[ix][iy][0];
      }
    }
  } else {
    cmac_assert(iz == _z_coords.size() - 1);
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
InterpolatedDensityFunction::InterpolatedDensityFunction(ParameterFile &params,
                                                         Log *log)
    : InterpolatedDensityFunction(
          params.get_value< std::string >("densityfunction:filename"),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "densityfunction:temperature", "8000. K"),
          log) {}

/**
 * @brief Function that gives the density for a given coordinate.
 *
 * @param position CoordinateVector specifying a coordinate position (in m).
 * @return Density at the given coordinate (in m^-3).
 */
DensityValues InterpolatedDensityFunction::
operator()(CoordinateVector<> position) const {
  cmac_assert(position.x() >= _x_bounds.first &&
              position.x() <= _x_bounds.second);
  cmac_assert(position.y() >= _y_bounds.first &&
              position.y() <= _y_bounds.second);
  cmac_assert(position.z() >= _z_bounds.first &&
              position.z() <= _z_bounds.second);

  unsigned int ix =
      Utilities::locate(position.x(), &_x_coords[0], _x_coords.size());
  unsigned int iy =
      Utilities::locate(position.y(), &_y_coords[0], _y_coords.size());
  unsigned int iz =
      Utilities::locate(position.z(), &_z_coords[0], _z_coords.size());

  // equations from https://en.wikipedia.org/wiki/Trilinear_interpolation

  double xd =
      (position.x() - _x_coords[ix]) / (_x_coords[ix + 1] - _x_coords[ix]);
  double omxd = 1. - xd;
  double yd =
      (position.y() - _y_coords[iy]) / (_y_coords[iy + 1] - _y_coords[iy]);
  double omyd = 1. - yd;
  double zd =
      (position.z() - _z_coords[iz]) / (_z_coords[iz + 1] - _z_coords[iz]);
  double omzd = 1. - zd;

  double c00 = _number_densities[ix][iy][iz] * omxd +
               _number_densities[ix + 1][iy][iz] * xd;
  double c01 = _number_densities[ix][iy][iz + 1] * omxd +
               _number_densities[ix + 1][iy][iz + 1] * xd;
  double c10 = _number_densities[ix][iy + 1][iz] * omxd +
               _number_densities[ix + 1][iy + 1][iz] * xd;
  double c11 = _number_densities[ix][iy + 1][iz + 1] * omxd +
               _number_densities[ix + 1][iy + 1][iz + 1] * xd;

  double c0 = c00 * omyd + c10 * yd;
  double c1 = c01 * omyd + c11 * yd;

  double number_density = c0 * omzd + c1 * zd;

  DensityValues values;
  values.set_number_density(number_density);
  values.set_temperature(_temperature);
  values.set_ionic_fraction(ION_H_n, 1.e-6);
  values.set_ionic_fraction(ION_He_n, 1.e-6);
  return values;
}
