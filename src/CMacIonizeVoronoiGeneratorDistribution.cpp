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
 * @file CMacIonizeVoronoiGeneratorDistribution.cpp
 *
 * @brief CMacIonizeVoronoiGeneratorDistribution implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CMacIonizeVoronoiGeneratorDistribution.hpp"
#include "HDF5Tools.hpp"
#include "ParameterFile.hpp"

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read from.
 */
CMacIonizeVoronoiGeneratorDistribution::CMacIonizeVoronoiGeneratorDistribution(
    ParameterFile &params)
    : _next_index(0) {
  std::string filename = params.get_value< std::string >(
      "densitygrid:voronoi_generator_distribution:filename");

  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

  // get the simulation box anchor
  CoordinateVector<> box_anchor =
      params.get_physical_vector< QUANTITY_LENGTH >("densitygrid:box_anchor");

  // units
  double unit_length_in_SI = 1.;
  if (HDF5Tools::group_exists(file, "/Units")) {
    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Units");
    double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit length in cgs (U_L)");
    unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    HDF5Tools::close_group(units);
  }

  // open the gas particle group
  HDF5Tools::HDF5Group particles = HDF5Tools::open_group(file, "/PartType0");

  // read the positions
  _positions =
      HDF5Tools::read_dataset< CoordinateVector<> >(particles, "Coordinates");

  HDF5Tools::close_group(particles);

  HDF5Tools::close_file(file);

  // unit conversion
  for (unsigned int i = 0; i < _positions.size(); ++i) {
    _positions[i][0] *= unit_length_in_SI;
    _positions[i][1] *= unit_length_in_SI;
    _positions[i][2] *= unit_length_in_SI;
    _positions[i] += box_anchor;
  }
}

/**
 * @brief Virtual destructor.
 */
CMacIonizeVoronoiGeneratorDistribution::
    ~CMacIonizeVoronoiGeneratorDistribution() {}

/**
 * @brief Get the number of unique positions returned by this distribution.
 *
 * @return Number of unique positions returned by this distribution.
 */
unsigned int
CMacIonizeVoronoiGeneratorDistribution::get_number_of_positions() const {
  return _positions.size();
}

/**
 * @brief Get the next position.
 *
 * @return Generator position (in m).
 */
CoordinateVector<> CMacIonizeVoronoiGeneratorDistribution::get_position() {
  cmac_assert(_next_index < _positions.size());

  const unsigned int next_index = _next_index;
  ++_next_index;
  return _positions[next_index];
}
