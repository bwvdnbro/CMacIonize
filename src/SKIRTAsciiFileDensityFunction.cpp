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
 * @param x_name Name of the x positions column.
 * @param y_name Name of the y positions column.
 * @param z_name Name of the z positions column.
 * @param rho_name Name of the density column.
 * @param He_name Name of the helium abundance column.
 * @param C_name Name of the carbon abundance column.
 * @param N_name Name of the nitrogen abundance column.
 * @param O_name Name of the oxygen abundance column.
 * @param Ne_name Name of the neon abundance column.
 * @param S_name Name of the sulphur abundance column.
 * @param log Log to write logging info to.
 */
SKIRTAsciiFileDensityFunction::SKIRTAsciiFileDensityFunction(
    const std::string filename, const std::string x_name,
    const std::string y_name, const std::string z_name,
    const std::string rho_name, const std::string He_name,
    const std::string C_name, const std::string N_name,
    const std::string O_name, const std::string Ne_name,
    const std::string S_name, Log *log)
    : _octree(nullptr) {

  if (log) {
    log->write_info("Initialising SKIRTAsciiFileDensityFunction from file \"",
                    filename, "\", reading fields \"", x_name, "\", \"", y_name,
                    "\", \"", z_name, "\" and \"", rho_name, "\".");
  }

  SKIRTAsciiFile file(filename);

  if (!file.has_column(x_name) || !file.is_quantity(x_name, QUANTITY_LENGTH)) {
    cmac_error("Could not initiate x positions from column \"%s\"!",
               x_name.c_str());
  }
  if (!file.has_column(y_name) || !file.is_quantity(y_name, QUANTITY_LENGTH)) {
    cmac_error("Could not initiate y positions from column \"%s\"!",
               y_name.c_str());
  }
  if (!file.has_column(z_name) || !file.is_quantity(z_name, QUANTITY_LENGTH)) {
    cmac_error("Could not initiate z positions from column \"%s\"!",
               z_name.c_str());
  }
  if (!file.has_column(rho_name) ||
      !file.is_quantity(rho_name, QUANTITY_DENSITY)) {
    cmac_error("Could not initiate density from column \"%s\"!",
               rho_name.c_str());
  }
  if (!file.has_column(He_name)) {
    cmac_error("Could not initiate helium abundance from column \"%s\"!",
               He_name.c_str());
  }
  if (!file.has_column(C_name)) {
    cmac_error("Could not initiate carbon abundance from column \"%s\"!",
               C_name.c_str());
  }
  if (!file.has_column(N_name)) {
    cmac_error("Could not initiate nitrogen abundance from column \"%s\"!",
               N_name.c_str());
  }
  if (!file.has_column(O_name)) {
    cmac_error("Could not initiate oxygen abundance from column \"%s\"!",
               O_name.c_str());
  }
  if (!file.has_column(Ne_name)) {
    cmac_error("Could not initiate neon abundance from column \"%s\"!",
               Ne_name.c_str());
  }
  if (!file.has_column(S_name)) {
    cmac_error("Could not initiate sulphur abundance from column \"%s\"!",
               S_name.c_str());
  }

#ifndef VARIABLE_ABUNDANCES
  if (log) {
    log->write_warning(
        "Code was not configured with ACTIVATE_VARIABLE_ABUNDANCES, so "
        "abundance values in SKIRTAsciiFileDensityFunction are ignored!");
  }
#endif

  const std::vector< double > &x = file.get_column(x_name);
  const std::vector< double > &y = file.get_column(y_name);
  const std::vector< double > &z = file.get_column(z_name);
  const std::vector< double > &rho = file.get_column(rho_name);
  const std::vector< double > &aHe = file.get_column(He_name);
  const std::vector< double > &aC = file.get_column(C_name);
  const std::vector< double > &aN = file.get_column(N_name);
  const std::vector< double > &aO = file.get_column(O_name);
  const std::vector< double > &aNe = file.get_column(Ne_name);
  const std::vector< double > &aS = file.get_column(S_name);

  CoordinateVector<> minpos(x[0], x[1], x[2]), maxpos(x[0], x[1], x[2]);
  _positions.resize(file.number_of_rows());
  _number_densities.resize(file.number_of_rows());
  _abundances.resize(file.number_of_rows());
  for (size_t i = 0; i < file.number_of_rows(); ++i) {
    _positions[i][0] = x[i];
    _positions[i][1] = y[i];
    _positions[i][2] = z[i];
    minpos = CoordinateVector<>::min(minpos, _positions[i]);
    maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    _number_densities[i] = rho[i] / PhysicalConstants::get_physical_constant(
                                        PHYSICALCONSTANT_PROTON_MASS);
#ifdef HAS_HELIUM
    _abundances[i].set_abundance(ELEMENT_He, aHe[i]);
#endif
#ifdef HAS_CARBON
    _abundances[i].set_abundance(ELEMENT_C, aC[i]);
#endif
#ifdef HAS_NITROGEN
    _abundances[i].set_abundance(ELEMENT_N, aN[i]);
#endif
#ifdef HAS_OXYGEN
    _abundances[i].set_abundance(ELEMENT_O, aO[i]);
#endif
#ifdef HAS_NEON
    _abundances[i].set_abundance(ELEMENT_Ne, aNe[i]);
#endif
#ifdef HAS_SULPHUR
    _abundances[i].set_abundance(ELEMENT_S, aS[i]);
#endif
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
 *  - x_name: Name of the x positions column (default: x-coordinate)
 *  - y_name: Name of the y positions column (default: y-coordinate)
 *  - z_name: Name of the z positions column (default: z-coordinate)
 *  - rho_name: Name of the density column (default: density)
 *  - He_name: Name of the helium abundance column (default: helium-abundance)
 *  - C_name: Name of the carbon abundance column (default: carbon-abundance)
 *  - N_name: Name of the nitrogen abundance column (default:
 * nitrogen-abundance)
 *  - O_name: Name of the oxygen abundance column (default: oxygen-abundance)
 *  - Ne_name: Name of the neon abundance column (default: neon-abundance)
 *  - S_name: Name of the sulphur abundance column (default: sulphur-abundance)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging information to.
 */
SKIRTAsciiFileDensityFunction::SKIRTAsciiFileDensityFunction(
    ParameterFile &params, Log *log)
    : SKIRTAsciiFileDensityFunction(
          params.get_filename("DensityFunction:filename"),
          params.get_value< std::string >("DensityFunction:x_name",
                                          "x-coordinate"),
          params.get_value< std::string >("DensityFunction:y_name",
                                          "y-coordinate"),
          params.get_value< std::string >("DensityFunction:z_name",
                                          "z-coordinate"),
          params.get_value< std::string >("DensityFunction:rho_name",
                                          "density"),
          params.get_value< std::string >("DensityFunction:He_name",
                                          "helium-abundance"),
          params.get_value< std::string >("DensityFunction:C_name",
                                          "carbon-abundance"),
          params.get_value< std::string >("DensityFunction:N_name",
                                          "nitrogen-abundance"),
          params.get_value< std::string >("DensityFunction:O_name",
                                          "oxygen-abundance"),
          params.get_value< std::string >("DensityFunction:Ne_name",
                                          "neon-abundance"),
          params.get_value< std::string >("DensityFunction:S_name",
                                          "sulphur-abundance"),
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
#ifdef VARIABLE_ABUNDANCES
  values.set_abundances(_abundances[closest]);
#endif

  return values;
}
