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

#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "OperatingSystem.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Utility functions that are not really related to a single class.
 */
namespace Utilities {

/**
 * @brief Get a random double precision floating point value in between 0 and 1.
 *
 * @return Random uniform double precision floating point value.
 */
inline double random_double() { return ((double)rand()) / ((double)RAND_MAX); }

/**
 * @brief Get a random double precision position in a box with origin (0, 0, 0)
 * and a unit side length.
 *
 * @return Random uniform double precision position.
 */
inline CoordinateVector<> random_position() {
  return CoordinateVector<>(random_double(), random_double(), random_double());
}

/**
 * @brief Get a random integer value within the given range.
 *
 * @param min_value Minimum value (inclusive).
 * @param max_value Maximum value (exclusive).
 * @return Uniform random value in the range [min_value, max_value[.
 */
inline unsigned int random_int(int min_value, int max_value) {
  return min_value + random_double() * (max_value - min_value);
}

/**
 * @brief Split a string of the form [str1, str2, str3] into its parts.
 *
 * @param value std::string having the form mentioned above.
 * @param str1 Variable to store the first part in.
 * @param str2 Variable to store the second part in.
 * @param str3 Variable to store the third part in.
 */
inline void split_string(const std::string &value, std::string &str1,
                         std::string &str2, std::string &str3) {
  size_t pos1 = value.find('[') + 1;
  size_t pos2 = value.find(',', pos1);
  str1 = value.substr(pos1, pos2 - pos1);
  pos1 = pos2 + 1;
  pos2 = value.find(',', pos1);
  str2 = value.substr(pos1, pos2 - pos1);
  pos1 = pos2 + 1;
  pos2 = value.find(']', pos1);
  str3 = value.substr(pos1, pos2 - pos1);
}

/**
 * @brief Convert the given std::string to a given template integer type.
 *
 * We use this version and not strtol, since our version also supports exponent
 * notation (e.g. '1e6').
 *
 * @param value std::string to convert.
 * @return Integer containing the value stored in the std::string.
 */
template < typename _integer_type_ >
_integer_type_ string_to_integer(const std::string &value) {

  if (value.size() == 0) {
    cmac_error("Cannot extract an integer from an empty string!");
  }

  unsigned int idx = 0;
  // strip trailing whitespace
  while (idx < value.size() && value[idx] == ' ') {
    ++idx;
  }
  if (idx == value.size()) {
    cmac_error("String does not contain an integer: \"%s\"!", value.c_str());
  }

  // check if a minus sign is present
  char sign = 1;
  if (value[idx] == '-') {
    sign = -1;
    ++idx;
  }

  // now parse the integer
  std::vector< char > digits;
  _integer_type_ ivalue = 0;
  if (value[idx] == '0' && (value[idx + 1] == 'x' || value[idx + 1] == 'X')) {
    // hexadecimal notation
    idx += 2;
    while (idx < value.size() && isxdigit(value[idx])) {
      if (isdigit(value[idx])) {
        // we here exploit the fact that ASCII representations of numbers have
        // successive binary codes
        digits.push_back(value[idx] - '0');
      } else {
        // a-f (A-F) also have successive codes
        const char cval = tolower(value[idx]);
        digits.push_back(10 + cval - 'a');
      }
      ++idx;
    }
    if (digits.size() == 0) {
      cmac_error("Cannot convert string to integer: \"%s\"!", value.c_str());
    }
    if (digits.size() * 4 >
        sizeof(_integer_type_) * 8 - std::is_signed< _integer_type_ >::value) {
      cmac_error("Integer value too large: \"%s\"!", value.c_str());
    }
    unsigned int base = 1;
    // we use a reverse iterator, since the least significant digit comes last
    for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
      ivalue += base * (*it);
      base *= 16;
    }
  } else {
    // normal notation
    while (idx < value.size() && isdigit(value[idx])) {
      // we here exploit the fact that ASCII representations of numbers have
      // successive binary codes
      digits.push_back(value[idx] - '0');
      ++idx;
    }
    unsigned int base = 1;
    // we use a reverse iterator, since the least significant digit comes last
    for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
      ivalue += base * (*it);
      base *= 10;
    }

    // check for an exponent
    if (idx < value.size() && (value[idx] == 'e' || value[idx] == 'E')) {
      ++idx;
      digits.clear();
      while (idx < value.size() && isdigit(value[idx])) {
        // we here exploit the fact that ASCII representations of numbers have
        // successive binary codes
        digits.push_back(value[idx] - '0');
        ++idx;
      }
      unsigned int exponent = 0;
      unsigned int base = 1;
      // we use a reverse iterator, since the least significant digit comes last
      for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
        exponent += base * (*it);
        base *= 10;
      }
      _integer_type_ extra_part = 1;
      for (unsigned int i = 0; i < exponent; ++i) {
        extra_part *= 10;
      }
      ivalue *= extra_part;
    }
  }

  ivalue *= sign;

  return ivalue;
}

/**
 * @brief Convert the given string to a variable of the given template type.
 *
 * @param value std::string value.
 * @return Variable of the given template type containing the parsed contents of
 * the given std::string.
 */
template < typename _datatype_ > _datatype_ convert(const std::string &value);

/**
 * @brief Convert the given string to a double precision floating point value.
 *
 * @param value std::string value.
 * @return Double precision floating point stored in the string.
 */
template <> inline double convert< double >(const std::string &value) {
  char *str_end;
  const double dvalue = strtod(value.c_str(), &str_end);
  if (str_end == value.c_str()) {
    cmac_error("Error converting \"%s\" to a floating point value!",
               value.c_str());
  }
  return dvalue;
}

/**
 * @brief Convert the given string to a floating point CoordinateVector.
 *
 * @param value std::string to convert.
 * @return CoordinateVector containing the components found.
 */
template <>
inline CoordinateVector<>
convert< CoordinateVector<> >(const std::string &value) {

  CoordinateVector<> vvalue;
  const int num_found = sscanf(value.c_str(), "[%lf,%lf,%lf]", &vvalue[0],
                               &vvalue[1], &vvalue[2]);
  if (num_found != 3) {
    cmac_error("Error converting \"%s\" to a floating point CoordinateVector!",
               value.c_str());
  }
  return vvalue;
}

/**
 * @brief Convert the given string to an integer value.
 *
 * @param value std::string value.
 * @return Integer stored in the string.
 */
template <> inline long int convert< long int >(const std::string &value) {
  const long int ivalue = string_to_integer< long int >(value);
  return ivalue;
}

/**
 * @brief Convert the given string to an unsigned long integer value.
 *
 * @param value std::string value.
 * @return Unsigned long integer stored in the string.
 */
template <>
inline unsigned long int
convert< unsigned long int >(const std::string &value) {
  const unsigned long int ivalue =
      string_to_integer< unsigned long int >(value);
  return ivalue;
}

/**
 * @brief Convert the given string to an unsigned char value.
 *
 * @param value std::string value.
 * @return Unsigned char stored in the string.
 */
template <>
inline unsigned char convert< unsigned char >(const std::string &value) {
  const unsigned char ivalue = string_to_integer< unsigned char >(value);
  return ivalue;
}

/**
 * @brief Convert the given string to an integer CoordinateVector.
 *
 * @param value std::string to convert.
 * @return CoordinateVector containing the components found.
 */
template <>
inline CoordinateVector< long int >
convert< CoordinateVector< long int > >(const std::string &value) {
  CoordinateVector< long int > vvalue;
  std::string x, y, z;
  split_string(value, x, y, z);
  vvalue[0] = convert< long int >(x);
  vvalue[1] = convert< long int >(y);
  vvalue[2] = convert< long int >(z);
  return vvalue;
}

/**
 * @brief Convert the given string to an unsigned integer CoordinateVector.
 *
 * @param value std::string to convert.
 * @return CoordinateVector containing the components found.
 */
template <>
inline CoordinateVector< unsigned long int >
convert< CoordinateVector< unsigned long int > >(const std::string &value) {
  CoordinateVector< unsigned long int > vvalue;
  std::string x, y, z;
  split_string(value, x, y, z);
  vvalue[0] = convert< unsigned long int >(x);
  vvalue[1] = convert< unsigned long int >(y);
  vvalue[2] = convert< unsigned long int >(z);
  return vvalue;
}

/**
 * @brief Convert the given string to a boolean value.
 *
 * The following string literals map to true: "true", "yes", "on", "y".
 * The following string literals map to false: "false", "no", "off", "n".
 * The string is converted to lowercase before it is parsed, so upper case or
 * mixed case versions, e.g. "True", "FALSE", "oFf" will also be correctly
 * parsed. All other string literals will result in an error.
 *
 * @param value std::string value.
 * @return True or false.
 */
template <> inline bool convert< bool >(const std::string &value) {

  std::string value_copy(value);
  // convert to lowercase
  std::transform(value_copy.begin(), value_copy.end(), value_copy.begin(),
                 ::tolower);
  // strip trailing whitespace
  unsigned int i = 0;
  while (value_copy[i] == ' ') {
    ++i;
  }
  value_copy = value_copy.substr(i);
  i = value_copy.size() - 1;
  while (value_copy[i] == ' ') {
    --i;
  }
  value_copy = value_copy.substr(0, i + 1);
  if (value_copy == "true" || value_copy == "yes" || value_copy == "on" ||
      value_copy == "y") {
    return true;
  } else if (value_copy == "false" || value_copy == "no" ||
             value_copy == "off" || value_copy == "n") {
    return false;
  } else {
    cmac_error("Error converting \"%s\" to a boolean value!",
               value_copy.c_str());
    return false;
  }
}

/**
 * @brief Convert the given string to a boolean CoordinateVector.
 *
 * @param value std::string to convert.
 * @return CoordinateVector containing the components found.
 */
template <>
inline CoordinateVector< bool >
convert< CoordinateVector< bool > >(const std::string &value) {
  CoordinateVector< bool > vvalue;
  std::string x, y, z;
  split_string(value, x, y, z);
  vvalue[0] = convert< bool >(x);
  vvalue[1] = convert< bool >(y);
  vvalue[2] = convert< bool >(z);
  return vvalue;
}

/**
 * @brief Convert the given value to a std::string.
 *
 * @param value Value to convert.
 * @return std::string.
 */
template < typename _datatype_ > std::string to_string(_datatype_ value) {
  std::stringstream sstream;
  sstream << value;
  return sstream.str();
}

/**
 * @brief to_string specialization for unsigned char values.
 *
 * We have to convert the unsigned char to an integer before outputting,
 * otherwise it is outputted as the character it is supposed to represent, which
 * yields garbage.
 *
 * @param value Value to convert.
 * @return std::string.
 */
template <> inline std::string to_string< unsigned char >(unsigned char value) {
  std::stringstream sstream;
  const unsigned int ivalue = value;
  sstream << ivalue;
  return sstream.str();
}

/**
 * @brief to_string specialization for boolean values.
 *
 * @param value Bool value.
 * @return "true" or "false".
 */
template <> inline std::string to_string< bool >(bool value) {
  if (value) {
    return "true";
  } else {
    return "false";
  }
}

/**
 * @brief to_string specialization for a floating point CoordinateVector.
 *
 * @param value Floating point CoordinateVector.
 * @return std::string containing the 3 components of the CoordinateVector.
 */
template <>
inline std::string
to_string< CoordinateVector<> >(const CoordinateVector<> value) {
  std::stringstream sstream;
  sstream << "[" << value.x() << ", " << value.y() << ", " << value.z() << "]";
  return sstream.str();
}

/**
 * @brief to_string specialization for an integer CoordinateVector.
 *
 * @param value Integer CoordinateVector.
 * @return std::string containing the 3 components of the CoordinateVector.
 */
template <>
inline std::string to_string< CoordinateVector< long int > >(
    const CoordinateVector< long int > value) {
  std::stringstream sstream;
  sstream << "[" << value.x() << ", " << value.y() << ", " << value.z() << "]";
  return sstream.str();
}

/**
 * @brief to_string specialization for an unsigned integer CoordinateVector.
 *
 * @param value Unsigned integer CoordinateVector.
 * @return std::string containing the 3 components of the CoordinateVector.
 */
template <>
inline std::string to_string< CoordinateVector< unsigned long int > >(
    const CoordinateVector< unsigned long int > value) {
  std::stringstream sstream;
  sstream << "[" << value.x() << ", " << value.y() << ", " << value.z() << "]";
  return sstream.str();
}

/**
 * @brief to_string specialization for a boolean CoordinateVector.
 *
 * @param value Boolean CoordinateVector.
 * @return std::string containing the 3 components of the CoordinateVector.
 */
template <>
inline std::string
to_string< CoordinateVector< bool > >(const CoordinateVector< bool > value) {
  std::stringstream sstream;
  sstream << "[" << to_string< bool >(value.x()) << ", "
          << to_string< bool >(value.y()) << ", "
          << to_string< bool >(value.z()) << "]";
  return sstream.str();
}

/**
 * @brief Split the given string containing a value and an associated unit into
 * a std::pair.
 *
 * @param svalue std::string containing a value - unit pair.
 * @return std::pair containing the value and unit.
 */
inline std::pair< double, std::string > split_value(const std::string &svalue) {
  size_t idx;
  double value;
  try {
    value = std::stod(svalue, &idx);
  } catch (std::invalid_argument e) {
    cmac_error("Error extracting value from \"%s\" unit-value pair!",
               svalue.c_str());
  }

  while (svalue[idx] == ' ') {
    ++idx;
  }
  const std::string unit = svalue.substr(idx);
  return make_pair(value, unit);
}

/**
 * @brief Get the index of the last element in the given ordered array that is
 * smaller than the given value.
 *
 * This routine uses bisection, and always returns a value in the range
 * [0, length-2], even if the given value is outside the given array.
 *
 * @param x Value to locate.
 * @param xarr Array in which to search.
 * @param length Length of the array.
 * @return Index of the last element in the ordered array that is smaller than
 * the given value, i.e. value is in between xarr[index] and xarr[index+1].
 */
inline unsigned int locate(double x, const double *xarr, unsigned int length) {
  unsigned int jl = 0;
  unsigned int ju = length;
  while (ju - jl > 1) {
    unsigned int jm = (ju + jl) / 2;
    if (x > xarr[jm]) {
      jl = jm;
    } else {
      ju = jm;
    }
  }
  if (jl == length - 1) {
    --jl;
  }
  return jl;
}

/**
 * @brief Compose a filename made up by the given prefix and counter value,
 * appropriately zero padded.
 *
 * @param folder Folder to add to the filename.
 * @param prefix Prefix for the filename.
 * @param extension Extension for the filename.
 * @param counter Value of the counter.
 * @param padding Number of digits the counter should have.
 * @return std::string with format: "<prefix>XX<counter>XX.<extension>", where
 * the number of Xs is equal to padding.
 */
inline std::string compose_filename(const std::string &folder,
                                    const std::string &prefix,
                                    const std::string &extension,
                                    unsigned int counter,
                                    unsigned int padding) {
  std::stringstream namestring;
  if (!folder.empty()) {
    namestring << folder << "/";
  }
  namestring << prefix;
  namestring.fill('0');
  namestring.width(padding);
  namestring << counter;
  namestring << "." << extension;
  return namestring.str();
}

/**
 * @brief Get the absolute path corresponding to the given path.
 *
 * @param path Path, can be relative or absolute.
 * @return Absolute path.
 */
inline std::string get_absolute_path(std::string path) {
  // strip trailing / from path
  if (path[path.size() - 1] == '/') {
    path = path.substr(0, path.size() - 1);
  }

  return OperatingSystem::absolute_path(path);
}

/**
 * @brief Get a time stamp of the form day/month/year, hour:minutes:seconds.
 *
 * @return Time stamp.
 */
inline std::string get_timestamp() {
  const std::time_t timestamp = std::time(nullptr);
  const std::tm *time = std::localtime(&timestamp);
  std::stringstream timestream;
  if (time->tm_mday < 10) {
    timestream << "0";
  }
  timestream << time->tm_mday << "/";
  if (time->tm_mon < 9) {
    timestream << "0";
  }
  // tm_mon counts from 0 to 11
  // tm_year is the number of years since 1900
  timestream << (time->tm_mon + 1) << "/" << (time->tm_year + 1900) << ", ";
  if (time->tm_hour < 10) {
    timestream << "0";
  }
  timestream << time->tm_hour << ":";
  if (time->tm_min < 10) {
    timestream << "0";
  }
  timestream << time->tm_min << ":";
  if (time->tm_sec < 10) {
    timestream << "0";
  }
  timestream << time->tm_sec;
  return timestream.str();
}

/**
 * @brief Convert a floating point time value in seconds to a human readable
 * time string.
 *
 * If the time is larger than a minute, it is displayed as Xm Ys, and so on for
 * larger subdivisions of time, up to years.
 *
 * @param time Double precision floating point time value.
 * @return std::string containing a human readable version of the time value.
 */
inline std::string human_readable_time(double time) {
  std::stringstream timestream;
  // 2^32 minutes is 8166 years. We can safely assume no internal timer will
  // ever reach that value
  unsigned int minutes = time / 60.;
  const double seconds = time - 60. * minutes;
  if (minutes > 0) {
    unsigned int hours = minutes / 60;
    if (hours > 0) {
      minutes -= 60 * hours;
      unsigned int days = hours / 24;
      if (days > 0) {
        hours -= 24 * days;
        const unsigned int years = days / 365;
        if (years > 0) {
          days -= 365 * years;
          timestream << years << "y ";
        }
        timestream << days << "d ";
      }
      timestream << hours << "h ";
    }
    timestream << minutes << "m ";
  }
  timestream << seconds << "s";
  return timestream.str();
}

/**
 * @brief Get the unit of @f$2^{10e}@f$ bytes.
 *
 * @param exponent Exponent @f$e@f$.
 * @return Name for @f$2^{10e}@f$ bytes: (@f$2^{10}@f$ bytes = KB, ...).
 */
inline std::string byte_unit(unsigned char exponent) {
  switch (exponent) {
  case 0:
    return "bytes";
  case 1:
    return "KB";
  case 2:
    return "MB";
  case 3:
    return "GB";
  case 4:
    return "TB";
  default:
    cmac_error("Name for 2^(10*%u) bytes was not implemented!", exponent);
    return "";
  }
}

/**
 * @brief Convert the given number of bytes to a human readable string.
 *
 * @param bytes Number of bytes.
 * @return std::string containing the given number of bytes in "bytes", "KB",
 * "MB", "GB"...
 */
inline std::string human_readable_bytes(unsigned long bytes) {
  unsigned char sizecount = 0;
  double bytefloat = bytes;
  while ((bytes >> 10) > 0) {
    bytes >>= 10;
    ++sizecount;
    bytefloat /= 1024.;
  }
  std::stringstream bytestream;
  bytefloat = std::round(100. * bytefloat) * 0.01;
  bytestream << bytefloat << " " << byte_unit(sizecount);
  return bytestream.str();
}

/**
 * @brief Check if the given std::string ends with the given std::string.
 *
 * @param haystack std::string to search in.
 * @param needle std::string to search.
 * @return True if haystack contains needle at the end.
 */
inline bool string_ends_with(const std::string &haystack,
                             const std::string &needle) {
  if (needle.size() > haystack.size()) {
    return false;
  }
  const size_t check = haystack.rfind(needle);
  // make sure we only flag needle at the end of the string
  return (check == haystack.size() - needle.size());
}

/**
 * @brief Decompose the given number into integer prime factors.
 *
 * @param number Number to decompose.
 * @return std::vector containing the prime integer components.
 */
inline std::vector< int > decompose(int number) {
  std::vector< int > components;
  int divisor = 2;
  while (divisor <= number && number > 1) {
    if (number % divisor == 0) {
      number /= divisor;
      components.push_back(divisor);
    } else {
      ++divisor;
    }
  }
  return components;
}

/**
 * @brief Subdivide the grid with the given resolution in approximately the
 * given number of equal size blocks.
 *
 * If such a decomposition is not possible, a larger number of blocks is used.
 *
 * @param ncell Grid resolution.
 * @param numblock Number of blocks to subdivide in.
 * @return Resolution of a single block.
 */
inline CoordinateVector< int > subdivide(CoordinateVector< int > ncell,
                                         int numblock) {

  std::vector< std::vector< int > > d_ncell(3);
  d_ncell[0] = decompose(ncell.x());
  d_ncell[1] = decompose(ncell.y());
  d_ncell[2] = decompose(ncell.z());

  // remove large factors until we reach at least the requested number of blocks
  int factor = 1;
  unsigned int idim = 0;
  while (factor < numblock) {
    factor *= d_ncell[idim].back();
    d_ncell[idim].pop_back();
    ++idim;
    if (idim == 3) {
      idim = 0;
    }
  }

  // collapse the remaining factors to get the block size
  CoordinateVector< int > block(1);
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < d_ncell[i].size(); ++j) {
      block[i] *= d_ncell[i][j];
    }
  }
  return block;
}

/**
 * @brief Convert the given long value to a binary sequence that can be written
 * out to a terminal or file.
 *
 * @param long_value Long value to convert.
 * @return std::string containing a binary representation of the long value.
 */
inline std::string as_binary_sequence(uint64_t long_value) {

  std::stringstream binary_stream;
  uint64_t mask = 0x8000000000000000;
  for (unsigned char i = 0; i < 64; ++i) {
    if (i > 0 && i % 4 == 0) {
      binary_stream << " ";
    }
    binary_stream << ((long_value & mask) >> (63 - i));
    mask >>= 1;
  }
  return binary_stream.str();
}

/**
 * @brief Return a std::vector that contains the indices that would sort the
 * given vector.
 *
 * @param values std::vector to sort.
 * @return std::vector containing the indices of the elements in the given
 * std::vector in an order that would sort the std::vector.
 */
template < typename _datatype_ >
std::vector< unsigned int > argsort(const std::vector< _datatype_ > &values) {
  std::vector< unsigned int > idx(values.size());
  for (unsigned int i = 0; i < values.size(); ++i) {
    idx[i] = i;
  }
  std::sort(idx.begin(), idx.end(), [&values](size_t i1, size_t i2) {
    return values[i1] < values[i2];
  });
  return idx;
}
}

#endif // UTILITIES_HPP
