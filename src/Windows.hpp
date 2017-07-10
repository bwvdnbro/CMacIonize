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
 * @file Windows.hpp
 *
 * @brief Windows specific implementation of OperatingSystem functionality.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef WINDOWS_HPP
#define WINDOWS_HPP

#include "OperatingSystem.hpp"

inline std::string OperatingSystem::absolute_path(std::string path) {
  // a maximum length of 1000 should be more than enough...
  char *absolute_path_ptr = _fullpath(nullptr, path.c_str(), 1000);
  if (absolute_path_ptr == nullptr) {
    cmac_error("Unable to resolve path \"%s\"!", path.c_str());
  }
  std::string absolute_path(absolute_path_ptr);
  free(absolute_path_ptr);
  return absolute_path;
}

#endif // WINDOWS_HPP
