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
#include <fstream>

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
