/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SKIRTAsciiFileDensityFunction.cpp
 *
 * @brief SKIRTAsciiFileDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 */
#include "SKIRTAsciiFileDensityFunction.hpp"
#include "Log.hpp"
#include "Octree.hpp"
#include "ParameterFile.hpp"
#include "SKIRTAsciiFile.hpp"

/**
 * @brief Constructor.
 *
 * @param filename Name of the ASCII text file to read.
 * @param xname Name of the x positions column.
 * @param yname Name of the y positions column.
 * @param zname Name of the z positions column.
 * @param rhoname Name of the density column.
 * @param log Log to write logging info to.
 */
SKIRTAsciiFileDensityFunction::SKIRTAsciiFileDensityFunction(
    const std::string filename, const std::string xname,
    const std::string yname, const std::string zname, const std::string rhoname,
    Log *log)
    : _octree(nullptr) {

  SKIRTAsciiFile file(filename);

  if (!file.has_column(xname) || !file.is_quantity(xname, QUANTITY_LENGTH)) {
    cmac_error("Could not initiate x positions from column \"%s\"!",
               xname.c_str());
  }
  if (!file.has_column(yname) || !file.is_quantity(yname, QUANTITY_LENGTH)) {
    cmac_error("Could not initiate y positions from column \"%s\"!",
               yname.c_str());
  }
  if (!file.has_column(zname) || !file.is_quantity(zname, QUANTITY_LENGTH)) {
    cmac_error("Could not initiate z positions from column \"%s\"!",
               zname.c_str());
  }
  if (!file.has_column(rhoname) ||
      !file.is_quantity(rhoname, QUANTITY_DENSITY)) {
    cmac_error("Could not initiate density from column \"%s\"!",
               rhoname.c_str());
  }

  const std::vector< double > &x = file.get_column(xname);
  const std::vector< double > &y = file.get_column(yname);
  const std::vector< double > &z = file.get_column(zname);
  const std::vector< double > &rho = file.get_column(rhoname);

  CoordinateVector<> minpos(x[0], x[1], x[2]), maxpos(x[0], x[1], x[2]);
  _positions.resize(file.number_of_rows());
  _number_densities.resize(file.number_of_rows());
  for (size_t i = 0; i < file.number_of_rows(); ++i) {
    _positions[i][0] = x[i];
    _positions[i][1] = y[i];
    _positions[i][2] = z[i];
    minpos = CoordinateVector<>::min(minpos, _positions[i]);
    maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    _number_densities[i] = rho[i] / PhysicalConstants::get_physical_constant(
                                        PHYSICALCONSTANT_PROTON_MASS);
  }

  // add some margins to the box to avoid problems with the Octree construction
  maxpos -= minpos;
  minpos = minpos - 0.005 * maxpos;
  maxpos = 1.01 * maxpos;
  Box<> box(minpos, maxpos);
  _octree = new Octree(_positions, box, false);

  if (log) {
    log->write_status("Successfully read in positions and densities from \"",
                      filename, "\".");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the ASCII file (required)
 *  - xname: Name of the x positions column (default: x-coordinate)
 *  - yname: Name of the y positions column (default: y-coordinate)
 *  - zname: Name of the z positions column (default: z-coordinate)
 *  - rhoname: Name of the density column (default: density)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging information to.
 */
SKIRTAsciiFileDensityFunction::SKIRTAsciiFileDensityFunction(
    ParameterFile &params, Log *log)
    : SKIRTAsciiFileDensityFunction(
          params.get_filename("DensityFunction:filename"),
          params.get_value< std::string >("DensityFuntion:xname",
                                          "x-coordinate"),
          params.get_value< std::string >("DensityFuntion:yname",
                                          "y-coordinate"),
          params.get_value< std::string >("DensityFuntion:zname",
                                          "z-coordinate"),
          params.get_value< std::string >("DensityFuntion:rhoname", "density"),
          log) {}

/**
 * @brief Destructor.
 *
 * Free memory used by the Octree.
 */
SKIRTAsciiFileDensityFunction::~SKIRTAsciiFileDensityFunction() {
  delete _octree;
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * We find the closests position in the input file and return those value.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues SKIRTAsciiFileDensityFunction::operator()(const Cell &cell) {

  DensityValues values;

  const uint_fast32_t closest =
      _octree->get_closest_ngb(cell.get_cell_midpoint());
  values.set_number_density(_number_densities[closest]);
  values.set_temperature(8000.);
  for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
    values.set_ionic_fraction(ion, 1.e-6);
  }

  return values;
}