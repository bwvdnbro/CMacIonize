/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2017 Maya Petkova (map32@st-andrews.ac.uk)
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
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */

#include "AsciiFileDensityGridWriter.hpp"
#include "DensityGrid.hpp"
#include "DensitySubGridCreator.hpp"
#include "Utilities.hpp"
#include <fstream>

/**
 * @brief Constructor.
 *
 * @param prefix Prefix of snapshot file names.
 * @param output_folder Name of the folder where output files should be placed.
 * @param log Log to write logging information to.
 */
AsciiFileDensityGridWriter::AsciiFileDensityGridWriter(
    std::string prefix, std::string output_folder, Log *log)
    : DensityGridWriter(output_folder, false, DensityGridWriterFields(false),
                        log),
      _prefix(prefix) {}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - prefix: Prefix that will be prepended to all snapshot file names (default:
 *    snapshot)
 *
 * @param output_folder Name of the folder where output files should be placed.
 * @param params ParameterFile to read.
 * @param log Log to write logging information to.
 */
AsciiFileDensityGridWriter::AsciiFileDensityGridWriter(
    std::string output_folder, ParameterFile &params, Log *log)
    : AsciiFileDensityGridWriter(params.get_value< std::string >(
                                     "DensityGridWriter:prefix", "snapshot"),
                                 output_folder, log) {}

/**
 * @brief Write a snapshot.
 *
 * @param grid DensityGrid to write out.
 * @param iteration Iteration number to use in the snapshot file name(s).
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 * @param hydro_units Internal unit system for the hydro.
 */
void AsciiFileDensityGridWriter::write(DensityGrid &grid,
                                       uint_fast32_t iteration,
                                       ParameterFile &params, double time,
                                       const InternalHydroUnits *hydro_units) {

  std::string filename =
      Utilities::compose_filename(_output_folder, _prefix, "txt", iteration, 3);
  std::ofstream file(filename);

  file << "#x (m)\ty (m)\tz (m)\tn (m^-3)\tvolume (m^3)\tneutral H fraction\n";

  for (auto it = grid.begin(); it != grid.end(); ++it) {
    CoordinateVector<> x = it.get_cell_midpoint();
    double n = it.get_ionization_variables().get_number_density();
    IonName ion = ION_H_n;
    double frac = it.get_ionization_variables().get_ionic_fraction(ion);
    double volume = it.get_volume();
    file << x.x() << "\t" << x.y() << "\t" << x.z() << "\t" << n << "\t"
         << volume << "\t" << frac << "\n";
  }
}

/**
 * @brief Write a snapshot for a split grid.
 *
 * @param grid_creator Grid.
 * @param counter Counter value to add to the snapshot file name.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 * @param hydro_units Internal unit system for the hydrodynamic quantities.
 */
void AsciiFileDensityGridWriter::write(DensitySubGridCreator &grid_creator,
                                       const uint_fast32_t counter,
                                       ParameterFile &params, double time,
                                       const InternalHydroUnits *hydro_units) {

  std::string filename =
      Utilities::compose_filename(_output_folder, _prefix, "txt", counter, 3);
  std::ofstream file(filename);

  file << "#x (m)\ty (m)\tz (m)\tn (m^-3)\tvolume (m^3)\tneutral H fraction\n";

  for (auto gridit = grid_creator.begin();
       gridit != grid_creator.original_end(); ++gridit) {
    for (auto cellit = (*gridit).begin(); cellit != (*gridit).end(); ++cellit) {
      const CoordinateVector<> x = cellit.get_cell_midpoint();
      const double n = cellit.get_number_density();
      const double xH = cellit.get_neutral_fraction();
      const double volume = cellit.get_volume();
      file << x.x() << "\t" << x.y() << "\t" << x.z() << "\t" << n << "\t"
           << volume << "\t" << xH << "\n";
    }
  }
}
