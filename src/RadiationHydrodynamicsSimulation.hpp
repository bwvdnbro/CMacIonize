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
 * @file RadiationHydrodynamicsSimulation.hpp
 *
 * @brief Radiation hydrodynamics (RHD) simulation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RADIATIONHYDRODYNAMICSSIMULATION_HPP
#define RADIATIONHYDRODYNAMICSSIMULATION_HPP

class CommandLineParser;
class Log;
class Timer;

/**
 * @brief Radiation hydrodynamics (RHD) simulation.
 */
class RadiationHydrodynamicsSimulation {
public:
  static int do_simulation(CommandLineParser &parser, bool write_output,
                           Timer &programtimer, Log *log = nullptr);
};

#endif // RADIATIONHYDRODYNAMICSSIMULATION_HPP
