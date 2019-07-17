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
 * @file EmissivityCalculationSimulation.hpp
 *
 * @brief Program to compute the emissivities based on a CMacIonize snapshot.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EMISSIVITYCALCULATIONSIMULATION_HPP
#define EMISSIVITYCALCULATIONSIMULATION_HPP

class CommandLineParser;
class Log;
class Timer;

/**
 * @brief Program to compute the emissivities based on a CMacIonize snapshot.
 */
class EmissivityCalculationSimulation {
public:
  static void add_command_line_parameters(CommandLineParser &parser);

  static int do_simulation(CommandLineParser &parser, bool write_output,
                           Timer &programtimer, Log *log = nullptr);
};

#endif // EMISSIVITYCALCULATIONSIMULATION_HPP
