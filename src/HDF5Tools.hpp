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
#include <cinttypes>
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
  const herr_t hdf5status = H5Eset_auto(nullptr, nullptr);
#else
  const herr_t hdf5status = H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
#endif
  if (hdf5status < 0) {
    cmac_error("Unable to turn off default HDF5 error handling!");
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
inline HDF5File open_file(std::string name, int_fast32_t mode) {

  hid_t file;
  if (mode == HDF5FILEMODE_READ) {
    file = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
      cmac_error("Unable to open file \"%s\"", name.c_str());
    }
  } else if (mode == HDF5FILEMODE_WRITE) {
    file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
      cmac_error("Unable to create file \"%s\"", name.c_str());
    }
  } else if (mode == HDF5FILEMODE_APPEND) {
    file = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
      cmac_error("Unable to open file \"%s\"", name.c_str());
    }
  } else {
    cmac_error("Unknown file mode: %" PRIuFAST32, mode);
  }

  return file;
}

/**
 * @brief Close the given open file.
 *
 * @param file HDF5File handle to an open file.
 */
inline void close_file(hid_t file) {
  const herr_t hdf5status = H5Fclose(file);
  if (hdf5status < 0) {
    cmac_error("Failed to close file!");
  }
}

/**
 * @brief Check whether a group with the given name exists in the given file.
 *
 * @param file HDF5File handle to an open file.
 * @param name Name of a group that might or might not exist.
 * @return True if the group exists.
 */
inline bool group_exists(hid_t file, std::string name) {
  const htri_t check = H5Lexists(file, name.c_str(), H5P_DEFAULT);
  return check > 0;
}

/**
 * @brief Create the HDF5 group with the given name in the given file.
 *
 * @param file HDF5File handle to an HDF5 file that is open in write mode.
 * @param name Name of the group to create.
 * @return HDF5Group handle to the newly created group.
 */
inline HDF5Group create_group(hid_t file, std::string name) {
#ifdef HDF5_OLD_API
  const hid_t group = H5Gcreate(file, name.c_str(), -1);
#else
  const hid_t group =
      H5Gcreate(file, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

  if (group < 0) {
    cmac_error("Unable to open group \"%s\"", name.c_str());
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
  const hid_t group = H5Gopen(file, name.c_str());
#else
  const hid_t group = H5Gopen(file, name.c_str(), H5P_DEFAULT);
#endif
  if (group < 0) {
    cmac_error("Unable to open group \"%s\"", name.c_str());
  }

  return group;
}

/**
 * @brief Close the given open group.
 *
 * @param group HDF5Group handle to an open group.
 */
inline void close_group(hid_t group) {
  const herr_t hdf5status = H5Gclose(group);
  if (hdf5status < 0) {
    cmac_error("Failed to close group!");
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
template < typename _datatype_ > inline hid_t get_datatype_name();

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
template <> inline hid_t get_datatype_name< uint32_t >() {
  return H5T_NATIVE_UINT32;
}

/**
 * @brief get_datatype_name specialization for a 64 bit unsigned integer.
 *
 * @return H5T_NATIVE_UINT64.
 */
template <> inline hid_t get_datatype_name< uint64_t >() {
  return H5T_NATIVE_UINT64;
}

/**
 * @brief get_datatype_name specialization for a 32 bit signed integer.
 *
 * @return H5T_NATIVE_INT32.
 */
template <> inline hid_t get_datatype_name< int32_t >() {
  return H5T_NATIVE_INT32;
}

/**
 * @brief Read the attribute with the given name of the given group.
 *
 * @param group HDF5Group handle to an open HDF5 group.
 * @param name Name of the attribute to read.
 * @return Value of the attribute.
 */
template < typename _datatype_ >
inline _datatype_ read_attribute(hid_t group, std::string name) {

  const hid_t datatype = get_datatype_name< _datatype_ >();
  // open attribute
  const hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    cmac_error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // read attribute
  _datatype_ value;
  herr_t hdf5status = H5Aread(attr, datatype, &value);
  if (hdf5status < 0) {
    cmac_error("Failed to read attribute \"%s\"!", name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
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
  const hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    cmac_error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // retrieve the length of the string
  H5A_info_t info;
  herr_t hdf5status = H5Aget_info(attr, &info);
  if (hdf5status < 0) {
    cmac_error("Failed to retrieve info for attribute \"%s\"!", name.c_str());
  }
  const hsize_t length = info.data_size;

  // C-string buffer to store the result in.
  char *data = new char[length + 1];

  // create C-string datatype
  const hid_t strtype = H5Tcopy(H5T_C_S1);
  if (strtype < 0) {
    cmac_error("Failed to create C-string datatype for attribute \"%s\"!",
               name.c_str());
  }

  // set datatype length to variable
  hdf5status = H5Tset_size(strtype, length + 1);
  if (hdf5status < 0) {
    cmac_error("Failed to set size of C-string datatype for attribute \"%s\"!",
               name.c_str());
  }

  // read attribute
  hdf5status = H5Aread(attr, strtype, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read string attribute \"%s\"!", name.c_str());
  }
  data[length] = '\0';

  // close string type
  hdf5status = H5Tclose(strtype);
  if (hdf5status < 0) {
    cmac_error("Failed to close C-string datatype for attribute \"%s\"!",
               name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
  }

  const std::string value(data);

  delete[] data;

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

  const hid_t datatype = get_datatype_name< double >();
  // open attribute
  const hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    cmac_error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // read attribute
  double data[3];
  herr_t hdf5status = H5Aread(attr, datatype, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read attribute \"%s\"!", name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
  }

  return CoordinateVector<>(data[0], data[1], data[2]);
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
template < typename _datatype_ >
inline std::vector< _datatype_ > read_vector_attribute(hid_t group,
                                                       std::string name) {

  const hid_t datatype = get_datatype_name< _datatype_ >();
  // open attribute
  const hid_t attr = H5Aopen(group, name.c_str(), H5P_DEFAULT);
  if (attr < 0) {
    cmac_error("Failed to open attribute \"%s\"!", name.c_str());
  }

  // open attribute dataspace
  const hid_t space = H5Aget_space(attr);
  if (space < 0) {
    cmac_error("Failed to open dataspace of attributes \"%s\"!", name.c_str());
  }

  // query dataspace size
  hsize_t size[1];
  hsize_t maxsize[1];
  const int_fast32_t ndim = H5Sget_simple_extent_dims(space, size, maxsize);
  if (ndim < 0) {
    cmac_error("Unable to query extent of attribute \"%s\"!", name.c_str());
  }
  if (!ndim) {
    size[0] = 1;
  }

  // read attribute
  _datatype_ *data = new _datatype_[size[0]];
  herr_t hdf5status = H5Aread(attr, datatype, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read attribute \"%s\"!", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(space);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of attribute \"%s\"!", name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
  }

  std::vector< _datatype_ > data_vector(size[0]);
  for (hsize_t i = 0; i < size[0]; ++i) {
    data_vector[i] = data[i];
  }

  delete[] data;

  return data_vector;
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
inline std::vector< uint32_t >
read_attribute< std::vector< uint32_t > >(hid_t group, std::string name) {
  return read_vector_attribute< uint32_t >(group, name);
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
 * @brief Add an attribute name to the list. Function used by H5Aiterate.
 *
 * @param group HDF5Group handle to the group we are reading.
 * @param c_name C-string name of the attribute.
 * @param info Attribute info.
 * @param data Extra data passed on to the function.
 * @return Error code used by H5Aiterate; just zero in our case.
 */
inline herr_t add_attribute_name(hid_t group, const char *c_name,
                                 const H5A_info_t *info, void *data) {

  std::vector< std::string > &names = *((std::vector< std::string > *)data);
  const std::string name(c_name);
  names.push_back(name);
  return 0;
}

/**
 * @brief Get the names of all the attributes of the given group.
 *
 * @param group HDF5Group handle to an open group.
 * @return std::vector containing the names of all attributes in the group.
 */
inline std::vector< std::string > get_attribute_names(hid_t group) {

  std::vector< std::string > names;
  const herr_t hdf5status = H5Aiterate(group, H5_INDEX_NAME, H5_ITER_INC,
                                       nullptr, add_attribute_name, &names);
  if (hdf5status < 0) {
    cmac_error("Failed to read attribute names for group!");
  }

  return names;
}

/**
 * @brief Write the attribute with the given name to the given group.
 *
 * @param group HDF5Group handle to an open HDF5 group.
 * @param name Name of the attribute to write.
 * @param value Value of the attribute.
 */
template < typename _datatype_ >
inline void write_attribute(hid_t group, std::string name, _datatype_ &value) {

  const hid_t datatype = get_datatype_name< _datatype_ >();
  // create dataspace
  const hid_t attspace = H5Screate(H5S_SCALAR);
  if (attspace < 0) {
    cmac_error("Failed to create dataspace for attribute \"%s\"!",
               name.c_str());
  }

// create attribute
#ifdef HDF5_OLD_API
  const hid_t attr =
      H5Acreate(group, name.c_str(), datatype, attspace, H5P_DEFAULT);
#else
  const hid_t attr = H5Acreate(group, name.c_str(), datatype, attspace,
                               H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0) {
    cmac_error("Failed to create attribute \"%s\"!", name.c_str());
  }

  // write attribute
  herr_t hdf5status = H5Awrite(attr, datatype, &value);
  if (hdf5status < 0) {
    cmac_error("Failed to write attribute \"%s\"!", name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(attspace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace for attribute \"%s\"!", name.c_str());
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
  const hid_t strtype = H5Tcopy(H5T_C_S1);
  if (strtype < 0) {
    cmac_error("Failed to copy C-string datatype for attribute \"%s\"!",
               name.c_str());
  }

  // set datatype length to length of string
  // note that we need to add an extra character for the string termination
  // character
  herr_t hdf5status = H5Tset_size(strtype, value.size() + 1);
  if (hdf5status < 0) {
    cmac_error("Failed to set size of C-string datatype for attribute \"%s\"!",
               name.c_str());
  }

  // create dataspace
  const hid_t attspace = H5Screate(H5S_SCALAR);
  if (attspace < 0) {
    cmac_error("Failed to create dataspace for attribute \"%s\"!",
               name.c_str());
  }

// create attribute
#ifdef HDF5_OLD_API
  const hid_t attr =
      H5Acreate(group, name.c_str(), strtype, attspace, H5P_DEFAULT);
#else
  const hid_t attr = H5Acreate(group, name.c_str(), strtype, attspace,
                               H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0) {
    cmac_error("Failed to create attribute \"%s\"!", name.c_str());
  }

  // write attribute
  hdf5status = H5Awrite(attr, strtype, value.c_str());
  if (hdf5status < 0) {
    cmac_error("Failed to write string attribute \"%s\"!", name.c_str());
  }

  // close string type
  hdf5status = H5Tclose(strtype);
  if (hdf5status < 0) {
    cmac_error("Failed to close C-string datatype for attribute \"%s\"!",
               name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(attspace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace for attribute \"%s\"!", name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
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

  const hid_t datatype = get_datatype_name< double >();
  // create dataspace
  const hsize_t dims[1] = {3};
  const hid_t attspace = H5Screate_simple(1, dims, nullptr);
  if (attspace < 0) {
    cmac_error("Failed to create dataspace for attribute \"%s\"!",
               name.c_str());
  }

// create attribute
#ifdef HDF5_OLD_API
  const hid_t attr =
      H5Acreate(group, name.c_str(), datatype, attspace, H5P_DEFAULT);
#else
  const hid_t attr = H5Acreate(group, name.c_str(), datatype, attspace,
                               H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0) {
    cmac_error("Failed to create attribute \"%s\"!", name.c_str());
  }

  // write attribute
  double data[3] = {value.x(), value.y(), value.z()};
  herr_t hdf5status = H5Awrite(attr, datatype, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write attribute \"%s\"!", name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(attspace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace for attribute \"%s\"!", name.c_str());
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
template < typename _datatype_ >
inline void write_vector_attribute(hid_t group, std::string name,
                                   std::vector< _datatype_ > &value) {

  const hid_t datatype = get_datatype_name< _datatype_ >();
  // create dataspace
  const hsize_t dims[1] = {value.size()};
  const hid_t attspace = H5Screate_simple(1, dims, nullptr);
  if (attspace < 0) {
    cmac_error("Failed to create dataspace for attribute \"%s\"!",
               name.c_str());
  }

// create attribute
#ifdef HDF5_OLD_API
  const hid_t attr =
      H5Acreate(group, name.c_str(), datatype, attspace, H5P_DEFAULT);
#else
  const hid_t attr = H5Acreate(group, name.c_str(), datatype, attspace,
                               H5P_DEFAULT, H5P_DEFAULT);
#endif
  if (attr < 0) {
    cmac_error("Failed to create attribute \"%s\"!", name.c_str());
  }

  // write attribute
  _datatype_ *data = new _datatype_[value.size()];
  for (size_t i = 0; i < value.size(); ++i) {
    data[i] = value[i];
  }
  herr_t hdf5status = H5Awrite(attr, datatype, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write attribute \"%s\"!", name.c_str());
  }

  // close attribute
  hdf5status = H5Aclose(attr);
  if (hdf5status < 0) {
    cmac_error("Failed to close attribute \"%s\"!", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(attspace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace for attribute \"%s\"!", name.c_str());
  }

  delete[] data;
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
inline void
write_attribute< std::vector< uint32_t > >(hid_t group, std::string name,
                                           std::vector< uint32_t > &value) {
  write_vector_attribute(group, name, value);
}

/**
 * @brief Read the dataset with the given name from the given group.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to read.
 * @return std::vector containing the contents of the dataset.
 */
template < typename _datatype_ >
inline std::vector< _datatype_ > read_dataset(hid_t group, std::string name) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  // open dataspace
  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to open dataspace of dataset \"%s\"", name.c_str());
  }

  // query dataspace extents
  hsize_t size[1];
  hsize_t maxsize[1];
  const int_fast32_t ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
  if (ndim < 0) {
    cmac_error("Unable to query extent of dataset \"%s\"", name.c_str());
  }

  // read dataset
  _datatype_ *data = new _datatype_[size[0]];
  herr_t hdf5status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  std::vector< _datatype_ > datavector(size[0]);
  for (hsize_t i = 0; i < size[0]; ++i) {
    datavector[i] = data[i];
  }

  delete[] data;

  return datavector;
}

/**
 * @brief Read part of the dataset with the given name from the given group.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to read.
 * @param part_offset Offset of the part that needs to be read.
 * @param part_size Size of the part that needs to be read.
 * @return std::vector containing the contents of the dataset.
 */
template < typename _datatype_ >
inline std::vector< _datatype_ >
read_dataset_part(const hid_t group, const std::string name,
                  const hsize_t part_offset, const hsize_t part_size) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  // open dataspace
  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to open dataspace of dataset \"%s\"", name.c_str());
  }

  // select the hyperslab in filespace we want to write to
  const hsize_t dims[1] = {part_size};
  const hsize_t offs[1] = {part_offset};
  herr_t hdf5status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offs,
                                          nullptr, dims, nullptr);
  if (hdf5status < 0) {
    cmac_error("Failed to select hyperslab in file space of dataset \"%s\"!",
               name.c_str());
  }

  // create memory space
  const hid_t memspace = H5Screate_simple(1, dims, nullptr);
  if (memspace < 0) {
    cmac_error("Failed to create memory space to write dataset \"%s\"!",
               name.c_str());
  }

  // read dataset
  _datatype_ *data = new _datatype_[part_size];
  hdf5status =
      H5Dread(dataset, datatype, memspace, filespace, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  std::vector< _datatype_ > datavector(part_size);
  for (hsize_t i = 0; i < part_size; ++i) {
    datavector[i] = data[i];
  }

  delete[] data;

  return datavector;
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

  const hid_t datatype = get_datatype_name< double >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  // open dataspace
  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to open dataspace of dataset \"%s\"", name.c_str());
  }

  // query dataspace extents
  hsize_t size[2];
  hsize_t maxsize[2];
  const int_fast32_t ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
  if (ndim < 0) {
    cmac_error("Unable to query extent of dataset \"%s\"", name.c_str());
  }

  // read dataset
  double *data = new double[size[0] * 3];
  herr_t hdf5status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  std::vector< CoordinateVector<> > datavector(size[0]);
  for (hsize_t i = 0; i < size[0]; ++i) {
    datavector[i][0] = data[3 * i];
    datavector[i][1] = data[3 * i + 1];
    datavector[i][2] = data[3 * i + 2];
  }

  delete[] data;

  return datavector;
}

/**
 * @brief Multidimensional data block.
 */
template < typename _datatype_, uint_fast8_t _size_ > class HDF5DataBlock {
private:
  /*! @brief Size of the multidimensional array. */
  std::array< size_t, _size_ > _size;

  /*! @brief Data. */
  _datatype_ *_data;

public:
  /**
   * @brief Constructor.
   *
   * @param dimensions Size of the data block in each dimension.
   * @param data Data for the block, should be an array with total size equal to
   * the product of the given dimensions.
   */
  HDF5DataBlock(std::array< size_t, _size_ > dimensions, _datatype_ *data)
      : _size(dimensions) {

    size_t datasize = 1;
    for (uint_fast8_t i = 0; i < _size_; ++i) {
      datasize *= _size[i];
    }
    _data = new _datatype_[datasize];
    for (size_t i = 0; i < datasize; ++i) {
      _data[i] = data[i];
    }
  }

  /**
   * @brief Assignment operator.
   *
   * We have to implement this ourselves, since the default assignment operator
   * does not know how to do the memory allocations.
   *
   * @param block HDF5DataBlock that is copied into this one.
   */
  void operator=(const HDF5DataBlock &block) {

    _size = block._size;
    size_t datasize = 1;
    for (uint_fast8_t i = 0; i < _size_; ++i) {
      datasize *= _size[i];
    }
    _data = new _datatype_[datasize];
    for (size_t i = 0; i < datasize; ++i) {
      _data[i] = block._data[i];
    }
  }

  /**
   * @brief Destructor.
   *
   * Deletes the internal data array.
   */
  ~HDF5DataBlock() { delete[] _data; }

  /**
   * @brief Access the element at the given position.
   *
   * @param index Multidimensional index.
   * @return Element at that position.
   */
  inline _datatype_ &operator[](std::array< size_t, _size_ > index) {

    size_t dataindex = 0;
    size_t product = 1;
    for (uint_fast8_t i = 0; i < _size_; ++i) {
      dataindex += index[_size_ - 1 - i] * product;
      product *= _size[_size_ - 1 - i];
    }
    return _data[dataindex];
  }

  /**
   * @brief Access the element at the given position.
   *
   * This version only allows read access.
   *
   * @param index Multidimensional index.
   * @return Element at that position.
   */
  inline const _datatype_ &
  operator[](std::array< size_t, _size_ > index) const {

    size_t dataindex = 0;
    size_t product = 1;
    for (uint_fast8_t i = 0; i < _size_; ++i) {
      dataindex += index[_size_ - 1 - i] * product;
      product *= _size[_size_ - 1 - i];
    }
    return _data[dataindex];
  }

  /**
   * @brief Get the size of the multidimensional array.
   *
   * @return Size of the array.
   */
  inline std::array< size_t, _size_ > size() const { return _size; }
};

/**
 * @brief read_dataset specialization for a HDF5DataBlock, a multidimensional
 * data array.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to read.
 * @return HDF5DataBlock containing the contents of the dataset.
 */
template < typename _datatype_, uint_fast8_t _size_ >
HDF5DataBlock< _datatype_, _size_ > read_dataset(hid_t group,
                                                 std::string name) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  // open dataspace
  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to open dataspace of dataset \"%s\"", name.c_str());
  }

  // query dataspace extents
  hsize_t size[_size_];
  hsize_t maxsize[_size_];
  const int_fast32_t ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
  if (ndim < 0) {
    cmac_error("Unable to query extent of dataset \"%s\"", name.c_str());
  }

  // read dataset
  std::array< size_t, _size_ > dimensions;
  size_t dprod = 1;
  for (uint_fast8_t i = 0; i < _size_; ++i) {
    dimensions[i] = size[i];
    dprod *= size[i];
  }
  _datatype_ *data = new _datatype_[dprod];
  herr_t hdf5status =
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  const HDF5DataBlock< _datatype_, _size_ > block(dimensions, data);

  delete[] data;

  return block;
}

/**
 * @brief Struct used to read in compound datasets consisting of a key and a
 * value, like in FLASH snapshots.
 */
template < typename _datatype_ > struct HDF5CompoundKeyValueType {
  /*! @brief Key name. */
  char _name[20];
  /*! @brief Value. */
  _datatype_ _value;
};

/**
 * @brief Wrapper for std::map that checks if accessed elements exist.
 */
template < typename _datatype_ > class HDF5Dictionary {
private:
  /*! @brief std::map for which this class is a wrapper. */
  std::map< std::string, _datatype_ > _map;

public:
  /**
   * @brief Constructor.
   *
   * @param map std::map for which this class is a wrapper.
   */
  inline HDF5Dictionary(std::map< std::string, _datatype_ > &map) : _map(map) {}

  /**
   * @brief Access operator.
   *
   * Contrary to a real std::map, this operator checks if the requested element
   * exists.
   *
   * @param key Key in the dictionary that we want to access.
   * @return Element belonging to that key.
   */
  inline _datatype_ &operator[](std::string key) {
    const auto it = _map.find(key);
    if (it == _map.end()) {
      cmac_error("Element \"%s\" not found in dictionary!", key.c_str());
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
template < typename _datatype_ >
inline HDF5Dictionary< _datatype_ > read_dictionary(hid_t group,
                                                    std::string name) {

  const hid_t valuetype = get_datatype_name< _datatype_ >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  // open dataspace
  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to open dataspace of dataset \"%s\"", name.c_str());
  }

  // query dataspace extents
  hsize_t size[1];
  hsize_t maxsize[1];
  const int_fast32_t ndim = H5Sget_simple_extent_dims(filespace, size, maxsize);
  if (ndim < 0) {
    cmac_error("Unable to query extent of dataset \"%s\"", name.c_str());
  }

  // create compound data type
  const hid_t datatype =
      H5Tcreate(H5T_COMPOUND, sizeof(HDF5CompoundKeyValueType< _datatype_ >));
  if (datatype < 0) {
    cmac_error("Failed to create datatype for dataset \"%s\"", name.c_str());
  }

  // set the contents of the compound data type
  const hid_t string20 = H5Tcopy(H5T_C_S1);
  herr_t hdf5status = H5Tset_size(string20, 20);
  if (hdf5status < 0) {
    cmac_error("Failed to initialize string type for dataset \"%s\"",
               name.c_str());
  }

  hdf5status = H5Tinsert(datatype, "name",
                         HOFFSET(HDF5CompoundKeyValueType< _datatype_ >, _name),
                         string20);
  if (hdf5status < 0) {
    cmac_error("Failed to insert name type for dataset \"%s\"", name.c_str());
  }
  hdf5status = H5Tinsert(
      datatype, "value",
      HOFFSET(HDF5CompoundKeyValueType< _datatype_ >, _value), valuetype);
  if (hdf5status < 0) {
    cmac_error("Failed to insert value type for dataset \"%s\"", name.c_str());
  }

  // read the data
  HDF5CompoundKeyValueType< _datatype_ > *data =
      new HDF5CompoundKeyValueType< _datatype_ >[size[0]];

  hdf5status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to read dataset \"%s\"", name.c_str());
  }

  // close the datatype
  hdf5status = H5Tclose(string20);
  if (hdf5status < 0) {
    cmac_error("Failed to close string type for dataset \"%s\"", name.c_str());
  }

  hdf5status = H5Tclose(datatype);
  if (hdf5status < 0) {
    cmac_error("Failed to close data type for dataset \"%s\"", name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  // construct the dictionary
  std::map< std::string, _datatype_ > dictionary;
  for (hsize_t i = 0; i < size[0]; ++i) {
    // strip spaces at the end of the string
    uint_fast8_t j = 18;
    while (data[i]._name[j] == ' ') {
      data[i]._name[j] = '\0';
      --j;
    }
    const std::string key(data[i]._name);
    dictionary[key] = data[i]._value;
  }

  // free memory
  delete[] data;

  return HDF5Dictionary< _datatype_ >(dictionary);
}

/**
 * @brief Write the dataset with the given name to the given group.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to write.
 * @param values std::vector containing the contents of the dataset.
 * @param compress Apply compression to the dataset?
 */
template < typename _datatype_ >
inline void write_dataset(hid_t group, std::string name,
                          std::vector< _datatype_ > &values,
                          const bool compress = false) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

  // create dataspace
  const uint_fast32_t vsize = values.size();
  const uint_fast32_t limit = 1 << 10;
  const hsize_t dims[1] = {vsize};
  const hsize_t chunk[1] = {std::min(vsize, limit)};
  const hid_t filespace = H5Screate_simple(1, dims, nullptr);
  if (filespace < 0) {
    cmac_error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

  // enable data compression
  const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  herr_t hdf5status = H5Pset_chunk(prop, 1, chunk);
  if (hdf5status < 0) {
    cmac_error("Failed to set chunk size for dataset \"%s\"", name.c_str());
  }

  if (compress) {
    hdf5status = H5Pset_fletcher32(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set Fletcher32 filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_shuffle(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set shuffle filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_deflate(prop, 9);
    if (hdf5status < 0) {
      cmac_error("Failed to set compression for dataset \"%s\"", name.c_str());
    }
  }

  // create dataset
#ifdef HDF5_OLD_API
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, prop);
#else
  const hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to create dataset \"%s\"", name.c_str());
  }

  // write dataset
  _datatype_ *data = new _datatype_[vsize];
  for (size_t i = 0; i < vsize; ++i) {
    data[i] = values[i];
  }
  hdf5status =
      H5Dwrite(dataset, datatype, H5S_ALL, filespace, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close creation properties
  hdf5status = H5Pclose(prop);
  if (hdf5status < 0) {
    cmac_error("Failed to close creation properties for dataset \"%s\"",
               name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  delete[] data;
}

/**
 * @brief write_dataset specialization for a CoordinateVector dataset.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to write.
 * @param values std::vector containing the contents of the dataset.
 * @param compress Apply compression to the dataset?
 */
template <>
inline void write_dataset(hid_t group, std::string name,
                          std::vector< CoordinateVector<> > &values,
                          const bool compress) {

  const hid_t datatype = get_datatype_name< double >();

  // create dataspace
  const uint_fast32_t vsize = values.size();
  const uint_fast32_t limit = 1 << 10;
  const hsize_t dims[2] = {vsize, 3};
  const hsize_t chunk[2] = {std::min(vsize, limit), 3};
  const hid_t filespace = H5Screate_simple(2, dims, nullptr);
  if (filespace < 0) {
    cmac_error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

  // enable data compression
  const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  herr_t hdf5status = H5Pset_chunk(prop, 2, chunk);
  if (hdf5status < 0) {
    cmac_error("Failed to set chunk size for dataset \"%s\"", name.c_str());
  }

  if (compress) {
    hdf5status = H5Pset_fletcher32(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set Fletcher32 filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_shuffle(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set shuffle filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_deflate(prop, 9);
    if (hdf5status < 0) {
      cmac_error("Failed to set compression for dataset \"%s\"", name.c_str());
    }
  }

// create dataset
#ifdef HDF5_OLD_API
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, prop);
#else
  const hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to create dataset \"%s\"", name.c_str());
  }

  // write dataset
  double *data = new double[3 * vsize];
  for (size_t i = 0; i < vsize; ++i) {
    data[3 * i] = values[i].x();
    data[3 * i + 1] = values[i].y();
    data[3 * i + 2] = values[i].z();
  }
  hdf5status =
      H5Dwrite(dataset, datatype, H5S_ALL, filespace, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close creation properties
  hdf5status = H5Pclose(prop);
  if (hdf5status < 0) {
    cmac_error("Failed to close creation properties for dataset \"%s\"",
               name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  delete[] data;
}

/**
 * @brief write_dataset specialization for a std::string dataset.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to write.
 * @param values std::vector containing the contents of the dataset.
 * @param compress Apply compression to the dataset?
 */
template <>
inline void write_dataset(hid_t group, std::string name,
                          std::vector< std::string > &values,
                          const bool compress) {

  const hid_t datatype = H5Tcopy(H5T_C_S1);
  if (datatype < 0) {
    cmac_error("Failed to create C-string datatype for attribute \"%s\"!",
               name.c_str());
  }

  // determine the maximum length of all strings
  const uint_fast32_t vsize = values.size();
  size_t length = 0;
  for (uint_fast32_t i = 0; i < vsize; ++i) {
    length = std::max(length, values[i].length());
  }

  // set datatype length to variable
  herr_t hdf5status = H5Tset_size(datatype, length + 1);
  if (hdf5status < 0) {
    cmac_error("Failed to set size of C-string datatype for dataset \"%s\"!",
               name.c_str());
  }

  // create dataspace
  const uint_fast32_t limit = 1 << 10;
  const hsize_t dims[1] = {vsize};
  const hsize_t chunk[1] = {std::min(vsize, limit)};
  const hid_t filespace = H5Screate_simple(1, dims, nullptr);
  if (filespace < 0) {
    cmac_error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

  // enable data compression
  const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  hdf5status = H5Pset_chunk(prop, 1, chunk);
  if (hdf5status < 0) {
    cmac_error("Failed to set chunk size for dataset \"%s\"", name.c_str());
  }

  if (compress) {
    hdf5status = H5Pset_fletcher32(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set Fletcher32 filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_shuffle(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set shuffle filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_deflate(prop, 9);
    if (hdf5status < 0) {
      cmac_error("Failed to set compression for dataset \"%s\"", name.c_str());
    }
  }

// create dataset
#ifdef HDF5_OLD_API
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, prop);
#else
  const hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to create dataset \"%s\"", name.c_str());
  }

  // write dataset
  char *data = new char[(length + 1) * vsize];
  for (size_t i = 0; i < vsize; ++i) {
    size_t j;
    for (j = 0; j < values[i].size(); ++j) {
      data[(length + 1) * i + j] = values[i][j];
    }
    data[(length + 1) * i + j] = '\0';
  }
  hdf5status =
      H5Dwrite(dataset, datatype, H5S_ALL, filespace, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close creation properties
  hdf5status = H5Pclose(prop);
  if (hdf5status < 0) {
    cmac_error("Failed to close creation properties for dataset \"%s\"",
               name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  delete[] data;
}

/**
 * @brief Create a new dataset with the given name and size in the given group.
 *
 * Once created, the dataset can be filled using HDF5Tools::append_dataset.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to create.
 * @param size Size of the dataset.
 * @param compress Apply compression to the dataset?
 */
template < typename _datatype_ >
inline void create_dataset(hid_t group, std::string name, hsize_t size,
                           const bool compress = false) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

  // create dataspace
  const hsize_t limit = 1 << 10;
  const hsize_t dims[1] = {size};
  const hsize_t chunk[1] = {std::min(size, limit)};
  const hid_t filespace = H5Screate_simple(1, dims, nullptr);
  if (filespace < 0) {
    cmac_error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

  // enable data compression
  const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  herr_t hdf5status = H5Pset_chunk(prop, 1, chunk);
  if (hdf5status < 0) {
    cmac_error("Failed to set chunk size for dataset \"%s\"", name.c_str());
  }

  if (compress) {
    hdf5status = H5Pset_fletcher32(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set Fletcher32 filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_shuffle(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set shuffle filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_deflate(prop, 9);
    if (hdf5status < 0) {
      cmac_error("Failed to set compression for dataset \"%s\"", name.c_str());
    }
  }

// create dataset
#ifdef HDF5_OLD_API
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, prop);
#else
  const hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to create dataset \"%s\"", name.c_str());
  }

  // close creation properties
  hdf5status = H5Pclose(prop);
  if (hdf5status < 0) {
    cmac_error("Failed to close creation properties for dataset \"%s\"",
               name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }
}

/**
 * @brief Create a new dataset with the given name and size in the given group.
 *
 * Template specialization for a dataset containing CoordinateVector<>s.
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to create.
 * @param size Size of the dataset.
 * @param compress Apply compression to the dataset?
 */
template <>
inline void create_dataset< CoordinateVector<> >(hid_t group, std::string name,
                                                 hsize_t size,
                                                 const bool compress) {

  const hid_t datatype = get_datatype_name< double >();

  // create dataspace
  const hsize_t limit = 1 << 10;
  const hsize_t dims[2] = {size, 3};
  const hsize_t chunk[2] = {std::min(size, limit), 3};
  const hid_t filespace = H5Screate_simple(2, dims, nullptr);
  if (filespace < 0) {
    cmac_error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

  // enable data compression
  const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  herr_t hdf5status = H5Pset_chunk(prop, 2, chunk);
  if (hdf5status < 0) {
    cmac_error("Failed to set chunk size for dataset \"%s\"", name.c_str());
  }

  if (compress) {
    hdf5status = H5Pset_fletcher32(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set Fletcher32 filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_shuffle(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set shuffle filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_deflate(prop, 9);
    if (hdf5status < 0) {
      cmac_error("Failed to set compression for dataset \"%s\"", name.c_str());
    }
  }

// create dataset
#ifdef HDF5_OLD_API
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, prop);
#else
  const hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to create dataset \"%s\"", name.c_str());
  }

  // close creation properties
  hdf5status = H5Pclose(prop);
  if (hdf5status < 0) {
    cmac_error("Failed to close creation properties for dataset \"%s\"",
               name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }
}

/**
 * @brief Append the given data to the dataset with the given name.
 *
 * @param group HDF5Group handle to an open group that contains the given
 * dataset.
 * @param name Name of the dataset to append to.
 * @param offset Offset of the new data within the dataset.
 * @param values std::vector containing the data to append.
 */
template < typename _datatype_ >
inline void append_dataset(hid_t group, std::string name, hsize_t offset,
                           std::vector< _datatype_ > &values) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to obtain file space of dataset \"%s\"!", name.c_str());
  }

  // select the hyperslab in filespace we want to write to
  const hsize_t dims[1] = {values.size()};
  const hsize_t offs[1] = {offset};
  herr_t hdf5status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offs,
                                          nullptr, dims, nullptr);
  if (hdf5status < 0) {
    cmac_error("Failed to select hyperslab in file space of dataset \"%s\"!",
               name.c_str());
  }

  // create memory space
  const hid_t memspace = H5Screate_simple(1, dims, nullptr);
  if (memspace < 0) {
    cmac_error("Failed to create memory space to write dataset \"%s\"!",
               name.c_str());
  }

  // write dataset
  _datatype_ *data = new _datatype_[values.size()];
  for (size_t i = 0; i < values.size(); ++i) {
    data[i] = values[i];
  }
  hdf5status =
      H5Dwrite(dataset, datatype, memspace, filespace, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close memory space
  hdf5status = H5Sclose(memspace);
  if (hdf5status < 0) {
    cmac_error("Failed to close memory space for dataset \"%s\"!",
               name.c_str());
  }

  // close file space
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close file space for dataset \"%s\"!", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  delete[] data;
}

/**
 * @brief Append the given data to the dataset with the given name.
 *
 * Template specialization for a CoordinateVector<> dataset.
 *
 * @param group HDF5Group handle to an open group that contains the given
 * dataset.
 * @param name Name of the dataset to append to.
 * @param offset Offset of the new data within the dataset.
 * @param values std::vector containing the data to append.
 */
template <>
inline void append_dataset(hid_t group, std::string name, hsize_t offset,
                           std::vector< CoordinateVector<> > &values) {

  const hid_t datatype = get_datatype_name< double >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to obtain file space of dataset \"%s\"!", name.c_str());
  }

  // select the hyperslab in filespace we want to write to
  const hsize_t dims[2] = {values.size(), 3};
  const hsize_t offs[2] = {offset, 0};
  herr_t hdf5status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offs,
                                          nullptr, dims, nullptr);
  if (hdf5status < 0) {
    cmac_error("Failed to select hyperslab in file space of dataset \"%s\"!",
               name.c_str());
  }

  // create memory space
  const hid_t memspace = H5Screate_simple(2, dims, nullptr);
  if (memspace < 0) {
    cmac_error("Failed to create memory space to write dataset \"%s\"!",
               name.c_str());
  }

  // write dataset
  double *data = new double[3 * values.size()];
  for (size_t i = 0; i < values.size(); ++i) {
    data[3 * i] = values[i].x();
    data[3 * i + 1] = values[i].y();
    data[3 * i + 2] = values[i].z();
  }
  hdf5status =
      H5Dwrite(dataset, datatype, memspace, filespace, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close memory space
  hdf5status = H5Sclose(memspace);
  if (hdf5status < 0) {
    cmac_error("Failed to close memory space for dataset \"%s\"!",
               name.c_str());
  }

  // close file space
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close file space for dataset \"%s\"!", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  delete[] data;
}

/**
 * @brief Create a new data table with the given name, number of rows and number
 * of columns in the given group.
 *
 * Once created, the data table can be filled using HDF5Tools::fill_row().
 *
 * @param group HDF5Group handle to an open group.
 * @param name Name of the dataset to create.
 * @param number_of_rows Number of rows in the dataset.
 * @param number_of_columns Number of columns in the dataset.
 * @param compress Apply compression to the dataset?
 */
template < typename _datatype_ >
inline void create_datatable(hid_t group, std::string name,
                             hsize_t number_of_rows, hsize_t number_of_columns,
                             const bool compress = false) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

  // create dataspace
  const hsize_t limit = 1 << 10;
  const hsize_t dims[2] = {number_of_rows, number_of_columns};
  const hsize_t chunk[2] = {std::min(number_of_rows, limit),
                            std::min(number_of_columns, limit)};
  const hid_t filespace = H5Screate_simple(2, dims, nullptr);
  if (filespace < 0) {
    cmac_error("Failed to create dataspace for dataset \"%s\"!", name.c_str());
  }

  // enable data compression
  const hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  herr_t hdf5status = H5Pset_chunk(prop, 2, chunk);
  if (hdf5status < 0) {
    cmac_error("Failed to set chunk size for dataset \"%s\"", name.c_str());
  }

  if (compress) {
    hdf5status = H5Pset_fletcher32(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set Fletcher32 filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_shuffle(prop);
    if (hdf5status < 0) {
      cmac_error("Failed to set shuffle filter for dataset \"%s\"",
                 name.c_str());
    }
    hdf5status = H5Pset_deflate(prop, 9);
    if (hdf5status < 0) {
      cmac_error("Failed to set compression for dataset \"%s\"", name.c_str());
    }
  }

// create dataset
#ifdef HDF5_OLD_API
  const hid_t dataset =
      H5Dcreate(group, name.c_str(), datatype, filespace, prop);
#else
  const hid_t dataset = H5Dcreate(group, name.c_str(), datatype, filespace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to create dataset \"%s\"", name.c_str());
  }

  // close creation properties
  hdf5status = H5Pclose(prop);
  if (hdf5status < 0) {
    cmac_error("Failed to close creation properties for dataset \"%s\"",
               name.c_str());
  }

  // close dataspace
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataspace of dataset \"%s\"", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }
}

/**
 * @brief Fill a row in the data table with the given name.
 *
 * @param group HDF5Group handle to an open group that contains the given
 * dataset.
 * @param name Name of the dataset to append to.
 * @param row_index Index of the row to fill.
 * @param values std::vector containing the data to append.
 */
template < typename _datatype_ >
inline void fill_row(hid_t group, std::string name, hsize_t row_index,
                     std::vector< _datatype_ > &values) {

  const hid_t datatype = get_datatype_name< _datatype_ >();

// open dataset
#ifdef HDF5_OLD_API
  const hid_t dataset = H5Dopen(group, name.c_str());
#else
  const hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
#endif
  if (dataset < 0) {
    cmac_error("Failed to open dataset \"%s\"", name.c_str());
  }

  const hid_t filespace = H5Dget_space(dataset);
  if (filespace < 0) {
    cmac_error("Failed to obtain file space of dataset \"%s\"!", name.c_str());
  }

  // select the hyperslab in filespace we want to write to
  const hsize_t dims[2] = {1, values.size()};
  const hsize_t offs[2] = {row_index, 0};
  herr_t hdf5status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offs,
                                          nullptr, dims, nullptr);
  if (hdf5status < 0) {
    cmac_error("Failed to select hyperslab in file space of dataset \"%s\"!",
               name.c_str());
  }

  // create memory space
  const hid_t memspace = H5Screate_simple(2, dims, nullptr);
  if (memspace < 0) {
    cmac_error("Failed to create memory space to write dataset \"%s\"!",
               name.c_str());
  }

  // write dataset
  _datatype_ *data = new _datatype_[values.size()];
  for (size_t i = 0; i < values.size(); ++i) {
    data[i] = values[i];
  }
  hdf5status =
      H5Dwrite(dataset, datatype, memspace, filespace, H5P_DEFAULT, data);
  if (hdf5status < 0) {
    cmac_error("Failed to write dataset \"%s\"", name.c_str());
  }

  // close memory space
  hdf5status = H5Sclose(memspace);
  if (hdf5status < 0) {
    cmac_error("Failed to close memory space for dataset \"%s\"!",
               name.c_str());
  }

  // close file space
  hdf5status = H5Sclose(filespace);
  if (hdf5status < 0) {
    cmac_error("Failed to close file space for dataset \"%s\"!", name.c_str());
  }

  // close dataset
  hdf5status = H5Dclose(dataset);
  if (hdf5status < 0) {
    cmac_error("Failed to close dataset \"%s\"", name.c_str());
  }

  delete[] data;
}

} // namespace HDF5Tools

#endif // HDF5TOOLS_HPP
