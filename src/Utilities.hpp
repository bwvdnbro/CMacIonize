/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Utilities.hpp
 *
 * @brief General functions that are used throughout the program and are not
 * part of the standard library.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cstdlib>

namespace Utilities {

/**
 * @brief Get a random double precision floating point value in between 0 and 1.
 *
 * @return Random uniform double precision floating point value.
 */
inline double random_double() { return ((double)rand()) / ((double)RAND_MAX); }
}

#endif // UTILITIES_HPP
