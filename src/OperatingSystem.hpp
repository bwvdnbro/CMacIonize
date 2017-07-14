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
 * @file OperatingSystem.hpp
 *
 * @brief Operating system specific functionality.
 *
 * This file only defines the namespace functions, the explicit implementation
 * is done in separate, operating system dependent files.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef OPERATINGSYSTEM_HPP
#define OPERATINGSYSTEM_HPP

#include "Configuration.hpp"

#include <string>

// first define the namespace...
/**
 * @brief Operating system specific functionality.
 */
namespace OperatingSystem {

/**
 * @brief Class used to store time information.
 */
class TimeValue;

/**
 * @brief Reset the given TimeValue object.
 *
 * @param time_value TimeValue object to reset.
 */
static void clear_time_value(TimeValue &time_value);

/**
 * @brief Set the value of the given TimeValue object to the current system
 * time value.
 *
 * @param time_value TimeValue object to set.
 */
static void get_time_value(TimeValue &time_value);

/**
 * @brief Subtract the given two TimeValue objects and store the result in the
 * third given TimeValue object.
 *
 * @param first_term First term.
 * @param second_term Second term.
 * @param result TimeValue object to store the result in.
 */
static void subtract_time_values(const TimeValue &first_term,
                                 const TimeValue &second_term,
                                 TimeValue &result);

/**
 * @brief Add the given two TimeValue objects and store the result in the third
 * given TimeValue object.
 *
 * @param first_term First term.
 * @param second_term Second term.
 * @param result TimeValue object to store the result in.
 */
static void add_time_values(const TimeValue &first_term,
                            const TimeValue &second_term, TimeValue &result);

/**
 * @brief Convert the given TimeValue object to a time in seconds.
 *
 * @param time_value TimeValue object to convert.
 * @return Time in seconds (with microsecond precision).
 */
static double convert_to_seconds(const TimeValue &time_value);

/**
 * @brief Get the absolute path corresponding to the given relative or absolute
 * path.
 *
 * @param path Path, can be either absolute or relative.
 * @return Absolute path.
 */
static std::string absolute_path(std::string path);
}

// ...then include the correct implementation
#ifdef HAVE_WINDOWS
#include "Windows.hpp"
#else
#include "Unix.hpp"
#endif

#endif // OPERATINGSYSTEM_HPP
