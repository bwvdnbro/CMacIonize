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

#include <array>
#include <hdf5.h>
#include <map>
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
  herr_t status = H5Eset_auto(nullptr, nullptr);
#else
  herr_t status = H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
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
}

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
}

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
}

/**
 * @brief Template read_attribute() version for std::vector.
 *
 * This function exists because we cannot specialize read_attribute() with a
 * template std::vector. Which means that we have to specialize read_attribute()
 * for every std::vector we want to use separately. These specializations all
 * call this function.
 *
 * If that is not clear: we cannot do this:
 * \code{.cpp}
 *  template<typename T>
 *  inline std::vector<T> read_attribute< std::vector<T> >(group, name){}
 * \endcode
 * because this type of template usage is not allowed (no idea why...).
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to read.
 * @return std::vector containing the values of the attribute.
 */
template < typename T >
inline std::vector< T > read_vector_attribute(hid_t group, std::string name) {
  hid_t datatype = get_datatype_name< T >();
  // open attribute
  hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // open attribute dataspace
  hid_t space = H5Aget_space(attr);
  if (space < 0) {
    error("Failed to open dataspace of attributes \"%s\"!", name.c_str());
  }

  // query dataspace size
  hsize_t size[1];
  hsize_t maxsize[1];
  int ndim = H5Sget_simple_extent_dims(space, size, maxsize);
  if (ndim < 0) {
    error("Unable to query extent of attribute \"%s\"!", name.c_str());
  }
  if (!ndim) {
    size[0] = 1;
  }

  // read attribute
  std::vector< T > value(size[0]);
  herr_t status = H5Aread(attr, datatype, &value[0]);
  if (status < 0) {
    error("Failed to read attribute \"%s\"!", name.c_str());
  }

  // close dataspace
  status = H5Sclose(space);
  if (status < 0) {
    error("Failed to close dataspace of attribute \"%s\"!", name.c_str());
  }

  // close attribute
  status = H5Aclose(attr);
  if (status < 0) {
    error("Failed to close attribute \"%s\"!", name.c_str());
  }

  return value;
}

/**
 * @brief read_attribute() specialization for std::vector<unsigned int>.
 *
 * This method calls read_vector_attribute().
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to read.
 * @return std::vector<unsigned int> containing the values of the attribute.
 */
template <>
inline std::vector< unsigned int >
read_attribute< std::vector< unsigned int > >(hid_t group, std::string name) {
  return read_vector_attribute< unsigned int >(group, name);
}

/**
 * @brief read_attribute() specialization for std::vector<double>.
 *
 * This method calls read_vector_attribute().
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to read.
 * @return std::vector<double> containing the values of the attribute.
 */
template <>
inline std::vector< double >
read_attribute< std::vector< double > >(hid_t group, std::string name) {
  return read_vector_attribute< double >(group, name);
}

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
}

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
}

/**
 * @brief write_attribute specialization for CoordinateVector.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to write.
 * @param value CoordinateVector containing the values to write.
 */
template <>
inline void write_attribute< CoordinateVector<> >(hid_t group, std::string name,
                                                  CoordinateVector<> &value) {
  hid_t datatype = get_datatype_name< double >();
  // create dataspace
  hsize_t dims[1] = {3};
  hid_t attspace = H5Screate_simple(1, dims, nullptr);
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
}

/**
 * @brief write_attribute specialization for a general std::vector.
 *
 * Just as for reading, we cannot use partial template specialization, and have
 * to do a workaround. Note that we could omit the <> and make it work without
 * the workaround, but then the read and write syntax would be different, which
 * could be confusing.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to write.
 * @param value std::vector containing the values to write.
 */
template < typename T >
inline void write_vector_attribute(hid_t group, std::string name,
                                   std::vector< T > &value) {
  hid_t datatype = get_datatype_name< T >();
  // create dataspace
  hsize_t dims[1] = {value.size()};
  hid_t attspace = H5Screate_simple(1, dims, nullptr);
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
  herr_t status = H5Awrite(attr, datatype, &value[0]);
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
}

/**
 * @brief write_attribute specialization for std::vector<double>.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to write.
 * @param value std::vector<double> containing the values to write.
 */
template <>
inline void
write_attribute< std::vector< double > >(hid_t group, std::string name,
                                         std::vector< double > &value) {
  write_vector_attribute(group, name, value);
}

/**
 * @brief write_attribute specialization for std::vector<unsigned int>.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the attribute to write.
 * @param value std::vector<double> containing the values to write.
 */
template <>
inline void write_attribute< std::vector< unsigned int > >(
    hid_t group, std::string name, std::vector< unsigned int > &value) {
  write_vector_attribute(group, name, value);
}

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
 * @brief Multidimensional data block.
 */
template < typename T, unsigned char Tsize > class HDF5DataBlock {
private:
  /*! @brief Size of the multidimensional array. */
  std::array< unsigned int, Tsize > _size;

  /*! @brief Data. */
  T *_data;

public:
  /**
   * @brief Constructor.
   *
   * @param dimensions Size of the data block in each dimension.
   * @param data Data for the block, should be an array with total size equal to
   * the product of the given dimensions.
   */
  HDF5DataBlock(std::array< unsigned int, Tsize > dimensions, T *data)
      : _size(dimensions) {
    unsigned int datasize = 1;
    for (unsigned char i = 0; i < Tsize; ++i) {
      datasize *= _size[i];
    }
    _data = new T[datasize];
    for (unsigned int i = 0; i < datasize; ++i) {
      _data[i] = data[i];
    }
  }

  ~HDF5DataBlock() { delete[] _data; }

  /**
   * @brief Access the element at the given position.
   *
   * @param index Multidimensional index.
   * @return Element at that position.
   */
  inline T &operator[](std::array< unsigned int, Tsize > index) {
    unsigned int dataindex = 0;
    unsigned int product = 1;
    for (unsigned char i = 0; i < Tsize; ++i) {
      dataindex += index[Tsize - 1 - i] * product;
      product *= _size[Tsize - 1 - i];
    }
    return _data[dataindex];
  }

  /**
   * @brief Get the size of the multidimensional array.
   *
   * @return Size of the array.
   */
  inline std::array< unsigned int, Tsize > size() { return _size; }
};

/**
 * @brief read_dataset specialization for a HDF5DataBlock, a multidimensional
 * data array.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to read.
 * @return HDF5DataBlock containing the contents of the dataset.
 */
template < typename T, unsigned char Tsize >
HDF5DataBlock< T, Tsize > read_dataset(hid_t group, std::string name) {
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
  hsize_t size[Tsize];
  hsize_t maxsize[Tsize];
  int ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
  if (ndim < 0) {
    error("Unable to query extent of dataset \"%s\"", name.c_str());
  }

  // read dataset
  std::array< unsigned int, Tsize > dimensions;
  unsigned int dprod = 1;
  for (unsigned char i = 0; i < Tsize; ++i) {
    dimensions[i] = size[i];
    dprod *= size[i];
  }
  T *data = new T[dprod];
  herr_t status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
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

  HDF5DataBlock< T, Tsize > block(dimensions, data);

  delete[] data;

  return block;
}

/**
 * @brief Struct used to read in compound datasets consisting of a key and a
 * value, like in FLASH snapshots.
 */
template < typename T > struct HDF5CompoundKeyValueType {
  /*! @brief Key name. */
  char _name[20];
  /*! @brief Value. */
  T _value;
};

/**
 * @brief Wrapper for std::map that checks if accessed elements exist.
 */
template < typename T > class HDF5Dictionary {
private:
  /*! @brief std::map for which this class is a wrapper. */
  std::map< std::string, T > _map;

public:
  /**
   * @brief Constructor.
   *
   * @param map std::map for which this class is a wrapper.
   */
  inline HDF5Dictionary(std::map< std::string, T > &map) : _map(map) {}

  /**
   * @brief Access operator.
   *
   * Contrary to a real std::map, this operator checks if the requested element
   * exists.
   *
   * @param key Key in the dictionary that we want to access.
   * @return Element belonging to that key.
   */
  inline T &operator[](std::string key) {
    auto it = _map.find(key);
    if (it == _map.end()) {
      error("Element \"%s\" not found in dictionary!", key.c_str());
    }
    return it->second;
  }
};

/**
 * @brief Read in a compound dataset as a dictionary.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the compound dataset to read.
 * @return Contents of the compound dataset, as a dictionary.
 */
template < typename T >
inline HDF5Dictionary< T > read_dictionary(hid_t group, std::string name) {
  hid_t valuetype = get_datatype_name< T >();

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

  // create compound data type
  hid_t datatype =
      H5Tcreate(H5T_COMPOUND, sizeof(HDF5CompoundKeyValueType< T >));
  if (datatype < 0) {
    error("Failed to create datatype for dataset \"%s\"", name.c_str());
  }

  // set the contents of the compound data type
  hid_t string20 = H5Tcopy(H5T_C_S1);
  herr_t status = H5Tset_size(string20, 20);
  if (status < 0) {
    error("Failed to initialize string type for dataset \"%s\"", name.c_str());
  }

  status = H5Tinsert(datatype, "name",
                     HOFFSET(HDF5CompoundKeyValueType< T >, _name), string20);
  if (status < 0) {
    error("Failed to insert name type for dataset \"%s\"", name.c_str());
  }
  status = H5Tinsert(datatype, "value",
                     HOFFSET(HDF5CompoundKeyValueType< T >, _value), valuetype);
  if (status < 0) {
    error("Failed to insert value type for dataset \"%s\"", name.c_str());
  }

  // read the data
  HDF5CompoundKeyValueType< T > *data =
      new HDF5CompoundKeyValueType< T >[ size[0] ];

  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (status < 0) {
    error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close the datatype
  status = H5Tclose(string20);
  if (status < 0) {
    error("Failed to close string type for dataset \"%s\"", name.c_str());
  }

  status = H5Tclose(datatype);
  if (status < 0) {
    error("Failed to close data type for dataset \"%s\"", name.c_str());
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

  // construct the dictionary
  std::map< std::string, T > dictionary;
  for (int i = 0; i < size[0]; ++i) {
    // strip spaces at the end of the string
    unsigned int j = 18;
    while (data[i]._name[j] == ' ') {
      data[i]._name[j] = '\0';
      --j;
    }
    std::string key(data[i]._name);
    dictionary[key] = data[i]._value;
  }

  // free memory
  delete[] data;

  return HDF5Dictionary< T >(dictionary);
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
  hid_t filespace = H5Screate_simple(1, dims, nullptr);
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
  hid_t filespace = H5Screate_simple(2, dims, nullptr);
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
