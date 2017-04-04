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
 * @param filename Name of the file containing the data.
 * @param log Log to write logging info to.
 */
InterpolatedDensityFunction::InterpolatedDensityFunction(std::string filename,
                                                         Log *log) {

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
  double xmin = yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("xmin");
  double xmax = yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("xmax");
  double ymin = yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("ymin");
  double ymax = yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("ymax");
  double zmin = yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("zmin");
  double zmax = yaml_dictionary.get_physical_value< QUANTITY_LENGTH >("zmax");
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
  unsigned int x_column = 0;
  unsigned int y_column = 0;
  unsigned int z_column = 0;
  unsigned int num_variable = 1;
  if (num_x != 0) {
    if (name_to_column.count("x") == 0) {
      cmac_error("No column found containing x values!");
    }
    x_column = name_to_column["x"];
    _x_coords.resize(num_x);
    _number_densities.resize(num_x);
    num_variable *= num_x;
  } else {
    _x_coords.resize(2);
    _x_coords[0] = xmin;
    _x_coords[1] = xmax;
    _number_densities.resize(2);
  }
  if (num_y != 0) {
    if (name_to_column.count("y") == 0) {
      cmac_error("No column found containing y values!");
    }
    y_column = name_to_column["y"];
    _y_coords.resize(num_y);
    for (unsigned int ix = 0; ix < _number_densities.size(); ++ix) {
      _number_densities[ix].resize(num_y);
    }
    num_variable *= num_y;
  } else {
    _y_coords.resize(2);
    _y_coords[0] = ymin;
    _y_coords[1] = ymax;
    for (unsigned int ix = 0; ix < _number_densities.size(); ++ix) {
      _number_densities[ix].resize(2);
    }
  }
  if (num_z != 0) {
    if (name_to_column.count("z") == 0) {
      cmac_error("No column found containing z values!");
    }
    z_column = name_to_column["z"];
    _z_coords.resize(num_z);
    for (unsigned int ix = 0; ix < _number_densities.size(); ++ix) {
      for (unsigned int iy = 0; iy < _number_densities[ix].size(); ++iy) {
        _number_densities[ix][iy].resize(num_z);
      }
    }
    num_variable *= num_z;
  } else {
    _z_coords.resize(2);
    _z_coords[0] = zmin;
    _z_coords[1] = zmax;
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
      cmac_assert(next_x >= xmin && next_x <= xmax);
      if (i > 0 && next_x != _x_coords[ix]) {
        ++ix;
        cmac_assert(ix < _x_coords.size());
      }
      _x_coords[ix] = next_x;
    }
    if (num_y != 0) {
      double next_y = UnitConverter::to_SI< QUANTITY_LENGTH >(row[y_column],
                                                              units[y_column]);
      cmac_assert(next_y >= ymin && next_y <= ymax);
      if (i > 0 && next_y != _y_coords[iy]) {
        ++iy;
        cmac_assert(iy < _y_coords.size());
      }
      _y_coords[iy] = next_y;
    }
    if (num_z != 0) {
      double next_z = UnitConverter::to_SI< QUANTITY_LENGTH >(row[z_column],
                                                              units[z_column]);
      cmac_assert(next_z >= zmin && next_z <= zmax);
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
  if (num_x == 0) {
    cmac_assert(ix == 0);
    for (unsigned int iy = 0; iy < _y_coords.size(); ++iy) {
      for (unsigned int iz = 0; iz < _z_coords.size(); ++iz) {
        _number_densities[1][iy][iz] = _number_densities[0][iy][iz];
      }
    }
  } else {
    cmac_assert(ix == _x_coords.size() - 1);
  }
  if (num_y == 0) {
    cmac_assert(iy == 0);
    for (unsigned int ix = 0; ix < _x_coords.size(); ++ix) {
      for (unsigned int iz = 0; iz < _z_coords.size(); ++iz) {
        _number_densities[ix][1][iz] = _number_densities[ix][0][iz];
      }
    }
  } else {
    cmac_assert(iy == _y_coords.size() - 1);
  }
  if (num_z == 0) {
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
          params.get_value< std::string >("densityfunction:filename"), log) {}

/**
 * @brief Function that gives the density for a given coordinate.
 *
 * @param position CoordinateVector specifying a coordinate position (in m).
 * @return Density at the given coordinate (in m^-3).
 */
DensityValues InterpolatedDensityFunction::
operator()(CoordinateVector<> position) const {
  cmac_assert(position.x() >= _x_coords[0] && position.x() <= _x_coords.back());
  cmac_assert(position.y() >= _y_coords[0] && position.y() <= _y_coords.back());
  cmac_assert(position.z() >= _z_coords[0] && position.z() <= _z_coords.back());

  unsigned int ix =
      Utilities::locate(position.x(), &_x_coords[0], _x_coords.size());
  unsigned int iy =
      Utilities::locate(position.y(), &_y_coords[0], _y_coords.size());
  unsigned int iz =
      Utilities::locate(position.z(), &_z_coords[0], _z_coords.size());

  double number_density = _number_densities[ix][iy][iz];

  DensityValues values;
  values.set_number_density(number_density);
  return values;
}
