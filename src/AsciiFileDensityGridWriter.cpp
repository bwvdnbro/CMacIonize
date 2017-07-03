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
 * @file AsciiFileDensityGridWriter.cpp
 *
 * @brief AsciiFileDensityGridWriter implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "AsciiFileDensityGridWriter.hpp"
#include "DensityGrid.hpp"
#include "Utilities.hpp"
#include <fstream>

/**
 * @brief Constructor.
 *
 * @param prefix Prefix of snapshot file names.
 * @param grid DensityGrid to write out.
 * @param output_folder Name of the folder where output files should be placed.
 * @param log Log to write logging information to.
 */
AsciiFileDensityGridWriter::AsciiFileDensityGridWriter(
    std::string prefix, DensityGrid &grid, std::string output_folder, Log *log)
    : DensityGridWriter(grid, output_folder, log), _prefix(prefix) {}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read.
 * @param grid DensityGrid to write out.
 * @param log Log to write logging information to.
 */
AsciiFileDensityGridWriter::AsciiFileDensityGridWriter(ParameterFile &params,
                                                       DensityGrid &grid,
                                                       Log *log)
    : AsciiFileDensityGridWriter(params.get_value< std::string >(
                                     "densitygridwriter:prefix", "snapshot"),
                                 grid, params.get_value< std::string >(
                                           "densitygridwriter:folder", "."),
                                 log) {}

/**
 * @brief Write a snapshot.
 *
 * @param iteration Iteration number to use in the snapshot file name(s).
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void AsciiFileDensityGridWriter::write(unsigned int iteration,
                                       ParameterFile &params, double time) {
  std::string filename =
      Utilities::compose_filename(_output_folder, _prefix, "txt", iteration, 3);
  std::ofstream file(filename);

  file << "#x (m)\ty (m)\tz (m)\tn (m^-3)\t neutral H fraction\n";   // Maya

  for (auto it = _grid.begin(); it != _grid.end(); ++it) {
    CoordinateVector<> x = it.get_cell_midpoint();
    double n = it.get_ionization_variables().get_number_density();        
    //IonName ion = static_cast< IonName >(0);     // Maya
    IonName ion = ION_H_n;
    double frac = it.get_ionization_variables().get_ionic_fraction(ion);    // Maya
    file << x.x() << "\t" << x.y() << "\t" << x.z() << "\t" << n << "\t" << frac <<"\n";   // Maya

  }
}
