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
 * @file IonizationSimulation.hpp
 *
 * @brief Ionization radiative transfer simulation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONSIMULATION_HPP
#define IONIZATIONSIMULATION_HPP

#include <string>

class CommandLineParser;
class Log;
class MPICommunicator;
class Timer;

/**
 * @brief Ionization radiative transfer simulation.
 */
class IonizationSimulation {
private:
  /*! @brief Log to write logging info to. */
  const Log *_log;

public:
  IonizationSimulation(const bool write_output,
                       const bool every_iteration_output, const int num_threads,
                       const std::string parameterfile, const bool dry_run,
                       MPICommunicator &comm, Log *log = nullptr);
};

#endif // IONIZATIONSIMULATION_HPP
