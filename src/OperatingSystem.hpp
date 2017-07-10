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
