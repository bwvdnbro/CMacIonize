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
 * @file HDF5Tools.hpp
 *
 * @brief Custom wrappers around some HDF5 library functions that feel more like
 * C++.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HDF5TOOLS_HPP
#define HDF5TOOLS_HPP

#include <hdf5.h>
#include <string>

/**
 * @brief Custom wrappers around some HDF5 library functions that feel more like
 * C++.
 */
namespace HDF5Tools {

/*! @brief More convenient name for a HDF5 file handle. */
typedef hid_t HDF5FileHandle;

/*! @brief Modes with which an HDF5 file can be opened. */
enum HDF5FileMode {
  /*! @brief Read mode (actually: read only). */
  HDF5FILEMODE_READ = 0,
  /*! @brief Write mode. */
  HDF5FILEMODE_WRITE
};

inline HDF5FileHandle open(std::string name, int mode) {
  hid_t file_mode;
  if (mode == HDF5FILEMODE_READ) {
    file_mode = H5F_ACC_RDONLY;
  } else if (mode == HDF5FILEMODE_WRITE) {
    file_mode = H5F_ACC_TRUNC;
  } else {
    //    error("Unknown file mode: %i", mode);
  }
  hid_t file = H5Fopen(name.c_str(), file_mode, H5P_DEFAULT);
  if (file > 0) {
    return file;
  } else {
    //    error("Unable to open file: %s", name.c_str());
    return -1;
  }
}
}

#endif // HDF5TOOLS_HPP
