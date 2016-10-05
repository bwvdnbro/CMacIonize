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

#include "CoordinateVector.hpp"
#include "Error.hpp"
#include <hdf5.h>
#include <string>
#include <vector>

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
  /*! @brief Write mode: file is created. */
  HDF5FILEMODE_WRITE,
  /*! @brief Append mode: file is opened for reading and writing. */
  HDF5FILEMODE_APPEND
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
inline HDF5File open_file(std::string name, int mode) {
  hid_t file;
  if (mode == HDF5FILEMODE_READ) {
    file = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      error("Unable to open file \"%s\"", name.c_str());
    }
  } else if (mode == HDF5FILEMODE_WRITE) {
    file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
      error("Unable to create file \"%s\"", name.c_str());
    }
  } else if (mode == HDF5FILEMODE_APPEND) {
    file = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
      error("Unable to open file \"%s\"", name.c_str());
    }
  } else {
    error("Unknown file mode: %i", mode);
  }

  return file;
}

/**
 * @brief Close the given open file.
 *
 * @param file HDF5File handle to an open file.
 */
inline void close_file(hid_t file) {
  herr_t status = H5Fclose(file);
  if (status < 0) {
    error("Failed to close file!");
  }
}

/**
 * @brief Create the HDF5 group with the given name in the given file.
 *
 * @param file HDF5File handle to an HDF5 file that is open in write mode.
 * @param name Name of the group to create.
 */
inline HDF5Group create_group(hid_t file, std::string name) {
#ifdef HDF5_OLD_API
  hid_t group = H5Gcreate(file, name.c_str(), -1);
#else
  hid_t group =
      H5Gcreate(file, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

  if (group < 0) {
    error("Unable to open group \"%s\"", name.c_str());
  }

  return group;
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
    error("Unable to open group \"%s\"", name.c_str());
  }

  return group;
}

/**
 * @brief Close the given open group.
 *
 * @param group HDF5Group handle to an open group.
 */
inline void close_group(hid_t group) {
  herr_t status = H5Gclose(group);
  if (status < 0) {
    error("Failed to close group!");
  }
}

/**
 * @brief Get the HDF5 data type corresponding to the template typename.
 *
 * This template function needs to be specialized for every typename that is
 * used.
 *
 * @return hid_t handle for the corresponding HDF5 data type.
 */
template < typename T > inline hid_t get_datatype_name();

/**
 * @brief get_datatype_name specialization for a double precision floating point
 * value.
 *
 * @return H5T_NATIVE_DOUBLE.
 */
template <> inline hid_t get_datatype_name< double >() {
  return H5T_NATIVE_DOUBLE;
}

/**
 * @brief get_datatype_name specialization for a single precision floating point
 * value.
 *
 * @return H5T_NATIVE_FLOAT.
 */
template <> inline hid_t get_datatype_name< float >() {
  return H5T_NATIVE_FLOAT;
}

/**
 * @brief get_datatype_name specialization for a 32 bit unsigned integer.
 *
 * @return H5T_NATIVE_UINT32.
 */
template <> inline hid_t get_datatype_name< unsigned int >() {
  return H5T_NATIVE_UINT32;
}

/**
 * @brief get_datatype_name specialization for a 32 bit signed integer.
 *
 * @return H5T_NATIVE_INT32.
 */
template <> inline hid_t get_datatype_name< int >() { return H5T_NATIVE_INT32; }

/**
 * @brief get_datatype_name specialization for a 64 bit unsigned integer.
 *
 * @return H5T_NATIVE_UINT64.
 */
template <> inline hid_t get_datatype_name< unsigned long long >() {
  return H5T_NATIVE_UINT64;
}

/**
 * @brief Read the attribute with the given name of the given group.
 *
 * @param group HDF5Group handle to an open HDF5 group.
 * @param name Name of the attribute to read.
 * @return Value of the attribute.
 */
template < typename T > inline T read_attribute(hid_t group, std::string name) {
  hid_t datatype = get_datatype_name< T >();
  // open attribute
  hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // read attribute
  T value;
  herr_t status = H5Aread(attr, datatype, &value);
  if (status < 0) {
    error("Failed to read attribute \"%s\"!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute \"%s\"!", name.c_str());
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
inline std::string read_attribute< std::string >(hid_t group,
                                                 std::string name) {
  // open attribute
  hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // retrieve the length of the string
  H5A_info_t info;
  herr_t status = H5Aget_info(attr, &info);
  if (status < 0) {
    error("Failed to retrieve info for attribute \"%s\"!", name.c_str());
  }
  hsize_t length = info.data_size;

  // C-string buffer to store the result in.
  char *data = (char *)malloc(length);

  // create C-string datatype
  hid_t strtype = H5Tcopy(H5T_C_S1);
  if (strtype < 0) {
    error("Failed to create C-string datatype for attribute \"%s\"!",
          name.c_str());
  }

  // set datatype length to variable
  status = H5Tset_size(strtype, length);
  if (status < 0) {
    error("Failed to set size of C-string datatype for attribute \"%s\"!",
          name.c_str());
  }

  // read attribute
  status = H5Aread(attr, strtype, data);
  if (status < 0) {
    error("Failed to read string attribute \"%s\"!", name.c_str());
  }

  // close string type
  status = H5Tclose(strtype);
  if (status < 0) {
    error("Failed to close C-string datatype for attribute \"%s\"!",
          name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute \"%s\"!", name.c_str());
  }

  std::string value(data);
  free(data);

  return value;
};

/**
 * @brief read_attribute specialization for CoordinateVector.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to read.
 * @return CoordinateVector containing the values of the attribute.
 */
template <>
inline CoordinateVector<>
read_attribute< CoordinateVector<> >(hid_t group, std::string name) {
  hid_t datatype = get_datatype_name< double >();
  // open attribute
  hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // read attribute
  CoordinateVector<> value;
  herr_t status = H5Aread(attr, datatype, &value);
  if (status < 0) {
    error("Failed to read attribute \"%s\"!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute \"%s\"!", name.c_str());
  }

  return value;
};

/**
 * @brief Write the attribute with the given name to the given group.
 *
 * @param group HDF5Group handle to an open HDF5 group.
 * @param name Name of the attribute to write.
 * @param value Value of the attribute.
 */
template < typename T >
inline void write_attribute(hid_t group, std::string name, T &value) {
  hid_t datatype = get_datatype_name< T >();
  // create dataspace
  hid_t attspace = H5Screate(H5S_SCALAR);
  if (attspace < 0) {
    error("Failed to create dataspace for attribute \"%s\"!", name.c_str());
  }

// create attribute
#ifdef HDF5_OLD_API
  hid_t attr = H5Acreate(group, name.c_str(), datatype, attspace, H5P_DEFAULT);
#else
  hid_t attr = H5Acreate(group, name.c_str(), datatype, attspace, H5P_DEFAULT,
                         H5P_DEFAULT);
#endif
  if (attr < 0) {
    error("Failed to create attribute \"%s\"!", name.c_str());
  }

  // write attribute
  herr_t status = H5Awrite(attr, datatype, &value);
  if (status < 0) {
    error("Failed to write attribute \"%s\"!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute \"%s\"!", name.c_str());
  }

  // close dataspace
  status = H5Sclose(attspace);
  if (status < 0) {
    error("Failed to close dataspace for attribute \"%s\"!", name.c_str());
  }
};

/**
 * @brief write_attribute specialization for std::string.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to write.
 * @param value std::string containing the value to write.
 */
template <>
inline void write_attribute< std::string >(hid_t group, std::string name,
                                           std::string &value) {
  // create C-string datatype
  hid_t strtype = H5Tcopy(H5T_C_S1);
  if (strtype < 0) {
    error("Failed to copy C-string datatype for attribute \"%s\"!",
          name.c_str());
  }

  // set datatype length to length of string
  // note that we need to add an extra character for the string termination
  // character
  herr_t status = H5Tset_size(strtype, value.size() + 1);
  if (status < 0) {
    error("Failed to set size of C-string datatype for attribute \"%s\"!",
          name.c_str());
  }

  // create dataspace
  hid_t attspace = H5Screate(H5S_SCALAR);
  if (attspace < 0) {
    error("Failed to create dataspace for attribute \"%s\"!", name.c_str());
  }

// create attribute
#ifdef HDF5_OLD_API
  hid_t attr = H5Acreate(group, name.c_str(), strtype, attspace, H5P_DEFAULT);
#else
  hid_t attr = H5Acreate(group, name.c_str(), strtype, attspace, H5P_DEFAULT,
                         H5P_DEFAULT);
#endif
  if (attr < 0) {
    error("Failed to create attribute \"%s\"!", name.c_str());
  }

  // write attribute
  status = H5Awrite(attr, strtype, value.c_str());
  if (status < 0) {
    error("Failed to write string attribute \"%s\"!", name.c_str());
  }

  // close string type
  status = H5Tclose(strtype);
  if (status < 0) {
    error("Failed to close C-string datatype for attribute \"%s\"!",
          name.c_str());
  }

  // close dataspace
  status = H5Sclose(attspace);
  if (status < 0) {
    error("Failed to close dataspace for attribute \"%s\"!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute \"%s\"!", name.c_str());
  }
};

/**
 * @brief write_attribute specialization for CoordinateVector.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to write.
 * @param value CoordinateVector containing the values to write.
 */
template < typename T >
inline void write_attribute(hid_t group, std::string name,
                            CoordinateVector<> &value) {
  hid_t datatype = get_datatype_name< double >();
  // create dataspace
  hsize_t dims[1] = {3};
  hid_t attspace = H5Screate_simple(1, dims, NULL);
  if (attspace < 0) {
    error("Failed to create dataspace for attribute \"%s\"!", name.c_str());
  }

// create attribute
#ifdef HDF5_OLD_API
  hid_t attr = H5Acreate(group, name.c_str(), datatype, attspace, H5P_DEFAULT);
#else
  hid_t attr = H5Acreate(group, name.c_str(), datatype, attspace, H5P_DEFAULT,
                         H5P_DEFAULT);
#endif
  if (attr < 0) {
    error("Failed to create attribute \"%s\"!", name.c_str());
  }

  // write attribute
  herr_t status = H5Awrite(attr, datatype, &value);
  if (status < 0) {
    error("Failed to write attribute \"%s\"!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute \"%s\"!", name.c_str());
  }

  // close dataspace
  status = H5Sclose(attspace);
  if (status < 0) {
    error("Failed to close dataspace for attribute \"%s\"!", name.c_str());
  }
};

/**
 * @brief Read the dataset with the given name from the given group.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to read.
 * @return std::vector containing the contents of the dataset.
 */
template < typename T >
inline std::vector< T > read_dataset(hid_t group, std::string name) {
  hid_t datatype = get_datatype_name< T >();

// open dataset
#ifdef HDF5_OLD_API
  hid_t dataset = H5Dopen(group, name.c_str());
#else
  hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    error("Failed to open dataset \"%s\"", name.c_str());
  }

  // open dataspace
  hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    error("Failed to open dataspace of dataset \"%s\"", name.c_str());
  }

  // query dataspace extents
  hsize_t size[1];
  hsize_t maxsize[1];
  int ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
  if (ndim < 0) {
    error("Unable to query extent of dataset \"%s\"", name.c_str());
  }

  // read dataset
  std::vector< T > data(size[0]);
  herr_t status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
  if (status < 0) {
    error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close dataspace
  status = H5Sclose(filespace);
  if (status < 0) {
    error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  status = H5Dclose(dataset);
  if (status < 0) {
    error("Failed to close dataset \"%s\"", name.c_str());
  }

  return data;
}

/**
 * @brief read_dataset specialization for a CoordinateVector dataset.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to read.
 * @return std::vector containing the contents of the dataset.
 */
template <>
inline std::vector< CoordinateVector<> >
read_dataset< CoordinateVector<> >(hid_t group, std::string name) {
  hid_t datatype = get_datatype_name< double >();

// open dataset
#ifdef HDF5_OLD_API
  hid_t dataset = H5Dopen(group, name.c_str());
#else
  hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    error("Failed to open dataset \"%s\"", name.c_str());
  }

  // open dataspace
  hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    error("Failed to open dataspace of dataset \"%s\"", name.c_str());
  }

  // query dataspace extents
  hsize_t size[2];
  hsize_t maxsize[2];
  int ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
  if (ndim < 0) {
    error("Unable to query extent of dataset \"%s\"", name.c_str());
  }

  // read dataset
  std::vector< CoordinateVector<> > data(size[0]);
  herr_t status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
  if (status < 0) {
    error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close dataspace
  status = H5Sclose(filespace);
  if (status < 0) {
    error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  status = H5Dclose(dataset);
  if (status < 0) {
    error("Failed to close dataset \"%s\"", name.c_str());
  }

  return data;
}

/**
 * @brief Write the dataset with the given name to the given group.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to write.
 * @param values std::vector containing the contents of the dataset.
 */
template < typename T >
inline void write_dataset(hid_t group, std::string name,
                          std::vector< T > &values) {
  hid_t datatype = get_datatype_name< T >();

  hsize_t dims[1] = {values.size()};
  // create dataspace
  hid_t filespace = H5Screate_simple(1, dims, NULL);
  if (filespace < 0) {
    error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

// create dataset
#ifdef HDF5_OLD_API
  hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, H5P_DEFAULT);
#else
  hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    error("Failed to create dataset \"%s\"", name.c_str());
  }

  // write dataset
  herr_t status =
      H5Dwrite(dataset, datatype, H5S_ALL, filespace, H5P_DEFAULT, &values[0]);
  if (status < 0) {
    error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close dataspace
  status = H5Sclose(filespace);
  if (status < 0) {
    error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  status = H5Dclose(dataset);
  if (status < 0) {
    error("Failed to close dataset \"%s\"", name.c_str());
  }
}

/**
 * @brief write_dataset specialization for a CoordinateVector dataset.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to write.
 * @param values std::vector containing the contents of the dataset.
 */
template <>
inline void write_dataset(hid_t group, std::string name,
                          std::vector< CoordinateVector<> > &values) {
  hid_t datatype = get_datatype_name< double >();

  hsize_t dims[2] = {values.size(), 3};
  // create dataspace
  hid_t filespace = H5Screate_simple(2, dims, NULL);
  if (filespace < 0) {
    error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

// create dataset
#ifdef HDF5_OLD_API
  hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, H5P_DEFAULT);
#else
  hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    error("Failed to create dataset \"%s\"", name.c_str());
  }

  // write dataset
  herr_t status =
      H5Dwrite(dataset, datatype, H5S_ALL, filespace, H5P_DEFAULT, &values[0]);
  if (status < 0) {
    error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close dataspace
  status = H5Sclose(filespace);
  if (status < 0) {
    error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  status = H5Dclose(dataset);
  if (status < 0) {
    error("Failed to close dataset \"%s\"", name.c_str());
  }
}
}

#endif // HDF5TOOLS_HPP
