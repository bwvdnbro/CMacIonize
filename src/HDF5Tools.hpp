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

#include "Error.hpp"
#include <hdf5.h>
#include <string>

/**
 * @brief Custom wrappers around some HDF5 library functions that feel more like
 * C++.
 */
namespace HDF5Tools {

/*! @brief More convenient name for a HDF5 file handle. */
typedef hid_t HDF5File;

/*! @brief More convenient name for a HDF5 group handle. */
typedef hid_t HDF5Group;

/*! @brief Modes with which an HDF5 file can be opened. */
enum HDF5FileMode {
  /*! @brief Read mode (actually: read only). */
  HDF5FILEMODE_READ = 0,
  /*! @brief Write mode. */
  HDF5FILEMODE_WRITE
};

/**
 * @brief Turn off default HDF5 error handling.
 */
inline void initialize() {
#ifdef HDF5_OLD_API
  herr_t status = H5Eset_auto(NULL, NULL);
#else
  herr_t status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
#endif
  if (status < 0) {
    error("Unable to turn off default HDF5 error handling!");
  }
}

/**
 * @brief Open the file with the given name for reading and/or writing.
 *
 * @param name Name of the file to open.
 * @param mode HDF5FileMode specifying the way the file will be used. Possible
 * values are HDF5FILEMODE_READ for reading and HDF5FILEMODE_WRITE for writing.
 * @return HDF5File handle to the open file that can be used by other methods.
 */
inline HDF5File open(std::string name, int mode) {
  hid_t file_mode;
  if (mode == HDF5FILEMODE_READ) {
    file_mode = H5F_ACC_RDONLY;
  } else if (mode == HDF5FILEMODE_WRITE) {
    file_mode = H5F_ACC_TRUNC;
  } else {
    error("Unknown file mode: %i", mode);
  }
  hid_t file = H5Fopen(name.c_str(), file_mode, H5P_DEFAULT);
  if (file < 0) {
    error("Unable to open file: %s", name.c_str());
  }

  return file;
}

/**
 * @brief Open the HDF5 group with the given name from the given file.
 *
 * @param file HDF5File handle to an open HDF5 file.
 * @param name Name of the group to open.
 * @return HDF5Group handle to the open group that can be used by other methods.
 */
inline HDF5Group open_group(hid_t file, std::string name) {
#ifdef HDF5_OLD_API
  hid_t group = H5Gopen(file, name.c_str());
#else
  hid_t group = H5Gopen(file, name.c_str(), H5P_DEFAULT);
#endif
  if (group < 0) {
    error("Unable to open group: %s", name.c_str());
  }

  return group;
}

/**
 * @brief Get the HDF5 data type corresponding to the template typename.
 *
 * This template function needs to be specialized for every typename that is
 * used.
 *
 * @return hid_t handle for the corresponding HDF5 data type.
 */
template <typename T> inline hid_t get_datatype_name();

/**
 * @brief get_datatype_name specialization for a double precision floating point
 * value.
 *
 * @return H5T_NATIVE_DOUBLE.
 */
template <> inline hid_t get_datatype_name<double>() {
  return H5T_NATIVE_DOUBLE;
}

/**
 * @brief Read the attribute with the given name of the given group.
 *
 * @param group HDF5Group handle to an open HDF5 group.
 * @param name Name of the attribute to read.
 * @return Value of the attribute.
 */
template <typename T> inline T read_attribute(hid_t group, std::string name) {
  hid_t datatype = get_datatype_name<double>();
  // open attribute
  hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    error("Failed to open attribute %s!", name.c_str());
  }

  // read attribute
  T value;
  herr_t status = H5Aread(attr, datatype, &value);
  if (status < 0) {
    error("Failed to read attribute %s!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute %s!", name.c_str());
  }

  return value;
};

/**
 * @brief read_attribute specialization for std::string.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to read.
 * @return std::string containing the value of the attribute.
 */
template <>
inline std::string read_attribute<std::string>(hid_t group, std::string name) {
  // open attribute
  hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    error("Failed to open attribute %s!", name.c_str());
  }

  // retrieve the length of the string
  H5A_info_t info;
  herr_t status = H5Aget_info(attr, &info);
  if (status < 0) {
    error("Failed to retrieve info for attribute %s!", name.c_str());
  }
  hsize_t length = info.data_size;

  // C-string buffer to store the result in.
  char *data = (char *)malloc(length);

  // create C-string datatype
  hid_t strtype = H5Tcopy(H5T_C_S1);
  if (strtype < 0) {
    error("Failed to create C-string datatype for attribute %s!", name.c_str());
  }

  // set datatype length to variable
  status = H5Tset_size(strtype, length);
  if (status < 0) {
    error("Failed to set size of C-string datatype for attribute %s!",
          name.c_str());
  }

  // read attribute
  status = H5Aread(attr, strtype, data);
  if (status < 0) {
    error("Failed to read string attribute %s!", name.c_str());
  }

  // close string type
  status = H5Tclose(strtype);
  if (status < 0) {
    error("Failed to close C-string datatype for attribute %s!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute %s!", name.c_str());
  }

  std::string value(data);
  free(data);

  return value;
};
}

#endif // HDF5TOOLS_HPP
