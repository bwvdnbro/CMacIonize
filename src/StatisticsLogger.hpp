/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file StatisticsLogger.hpp
 *
 * @brief Runtime log output for statistical properties of the simulation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef STATISTICSLOGGER_HPP
#define STATISTICSLOGGER_HPP

#include "DensityGrid.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"

#include <fstream>
#include <unistd.h>

/**
 * @brief Runtime log output for statistical properties of the simulation.
 */
class StatisticsLogger {
private:
  /*! @brief Log file to write to. */
  std::ofstream _logfile;

public:
  /**
   * @brief Constructor.
   */
  inline StatisticsLogger() : _logfile("StatisticsLogger.txt") {
    _logfile << "time (s)\ttotal mass (kg)\ttotal momentum x (kg m s^-1)\t"
                "total momentum y (kg m s^-1)\ttotal momentum z (kg m s^-1)\t"
                "total energy (kg m^2 s^-2)\n";
  }

  /**
   * @brief Write statistics for the given grid at the given time.
   *
   * @param time Current simulation time (s).
   * @param grid DensityGrid to analyse.
   */
  inline void write_statistics(const double time, DensityGrid &grid) {

    double mass = 0.;
    CoordinateVector<> momentum(0., 0., 0.);
    double energy = 0.;

    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const HydroVariables hydro = it.get_hydro_variables();
      mass += hydro.get_conserved_mass();
      momentum += hydro.get_conserved_momentum();
      energy += hydro.get_conserved_total_energy();
    }

    _logfile << time << "\t" << mass << "\t" << momentum.x() << "\t"
             << momentum.y() << "\t" << momentum.z() << "\t" << energy << "\n";
    _logfile.flush();
  }

  /**
   * @brief Write the logger to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  inline void write_restart_file(RestartWriter &restart_writer) {

    const auto filepos = _logfile.tellp();
    restart_writer.write(filepos);
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline StatisticsLogger(RestartReader &restart_reader) {

    const std::streampos filepos = restart_reader.read< std::streampos >();
    // truncate the original file to the size we were at
    if (truncate("StatisticsLogger.txt", filepos) != 0) {
      cmac_error("Error while truncating output file!");
    }
    // now open the file in append mode
    _logfile = std::ofstream("StatisticsLogger.txt", std::ios_base::app);
  }
};

#endif // STATISTICSLOGGER_HPP
