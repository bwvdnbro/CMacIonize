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
 * @file SPHNGSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction implementation that reads a density field from an
 * SPHNG snapshot file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHNGSNAPSHOTDENSITYFUNCTION_HPP
#define SPHNGSNAPSHOTDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"

#include <fstream>
#include <map>
#include <sstream>
#include <vector>

class Log;
class Octree;
class ParameterFile;

/**
 * @brief DensityFunction implementation that reads a density field from an
 * SPHNG snapshot file.
 */
class SPHNGSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Positions of the SPH particles in the snapshot (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Masses of the SPH particles in the snapshot (in kg). */
  std::vector< double > _masses;

  /*! @brief Smoothing lengths of the SPH particles in the snapshot (in m). */
  std::vector< double > _smoothing_lengths;

  /*! @brief Box containing all particles (in m). */
  Box<> _partbox;

  /*! @brief Octree used to speed up neighbour finding. */
  Octree *_octree;

  /*! @brief Initial temperature of the gas (in K). */
  double _initial_temperature;

  /*! @brief Number of bins to use when writing particle statistics. */
  unsigned int _stats_numbin;

  /*! @brief Minimum distance to use for particle statistics. */
  double _stats_mindist;

  /*! @brief Maximum distance to use for particle statistics. */
  double _stats_maxdist;

  /*! @brief Name of the file with particle statistics. */
  std::string _stats_filename;

  /*! @brief Log to write logging info to. */
  Log *_log;

  static double kernel(const double q, const double h);

  static double full_integral_new(double phi, double r0, double R_0, double h);
  static double mass_contribution(const Cell &cell, const CoordinateVector<> particle, const double h);

  /**
   * @brief Skip a block from the given Fortran unformatted binary file.
   *
   * @param ifile Reference to an open Fortran unformatted binary file.
   */
  inline static void skip_block(std::ifstream &ifile) {
    unsigned int length1, length2;
    ifile.read(reinterpret_cast< char * >(&length1), sizeof(unsigned int));
    ifile.seekg(length1, std::ios_base::cur);
    ifile.read(reinterpret_cast< char * >(&length2), sizeof(unsigned int));
    if (length1 != length2) {
      cmac_error("Wrong block size!");
    }
  }

  /**
   * @brief Get the size of the given template datatype.
   *
   * @param value Reference to a value of the template datatype.
   * @return Size of the template datatype.
   */
  template < typename _datatype_ >
  inline static unsigned int get_size(_datatype_ &value) {
    return sizeof(_datatype_);
  }

  /**
   * @brief Get the size of the given template datatypes (recursively).
   *
   * @param value Next value in the list.
   * @param args Other values in the list.
   * @return Total size of all values in the list.
   */
  template < typename _datatype_, typename... _arguments_ >
  inline static unsigned int get_size(_datatype_ &value,
                                      _arguments_ &... args) {
    return sizeof(_datatype_) + get_size(args...);
  }

  /**
   * @brief Fill the given referenced parameter by reading from the given
   * Fortran unformatted binary file.
   *
   * @param ifile Reference to an open Fortran unformatted binary file.
   * @param value Next (and last) value to read from the file.
   */
  template < typename _datatype_ >
  inline static void read_value(std::ifstream &ifile, _datatype_ &value) {
    ifile.read(reinterpret_cast< char * >(&value), sizeof(value));
  }

  /**
   * @brief Fill the given referenced template parameters by reading from the
   * given Fortran unformatted binary file.
   *
   * @param ifile Reference to an open Fortran unformatted binary file.
   * @param value Next value to read from the file.
   * @param args Other values to read from the file.
   */
  template < typename _datatype_, typename... _arguments_ >
  inline static void read_value(std::ifstream &ifile, _datatype_ &value,
                                _arguments_ &... args) {
    ifile.read(reinterpret_cast< char * >(&value), sizeof(value));
    read_value(ifile, args...);
  }

  /**
   * @brief Read a block  from a Fortran unformatted binary file and fill the
   * given referenced template parameters with its contents.
   *
   * An error will be thrown if the size (in bytes) of all parameters does not
   * match the size of the block.
   *
   * @param ifile Reference to an open Fortran unformatted binary file.
   * @param args References to variables that should be filled with the contents
   * of the block (in the order they are passed to this routine).
   */
  template < typename... _arguments_ >
  inline static void read_block(std::ifstream &ifile, _arguments_ &... args) {
    unsigned int length1, length2;
    ifile.read(reinterpret_cast< char * >(&length1), sizeof(unsigned int));
    unsigned int blocksize = get_size(args...);
    if (length1 != blocksize) {
      cmac_error("Wrong number of variables passed on to read_block()! Block "
                 "size is %u, but size of variables is %u.",
                 length1, blocksize);
    }
    read_value(ifile, args...);
    ifile.read(reinterpret_cast< char * >(&length2), sizeof(unsigned int));
    if (length1 != length2) {
      cmac_error("Wrong block size!");
    }
  }

  /**
   * @brief Read a dictionary containing tag-value pairs from the given Fortran
   * unformatted binary file.
   *
   * This routine assumes a 3 block structure, whereby the first block contains
   * a single integer giving the number of elements in the second and third
   * block. The second block contains 16-byte tags, while the third block
   * contains a value of the given template type for each tag.
   * If the file is not tagged, the second block is absent. In this case, the
   * tagged flag should be set to false, and the tags will simply be "tag",
   * "tag1"...
   *
   * @param ifile Reference to an open Fortran unformatted binary file to read
   * from.
   * @param tagged Flag indicating whether the file is tagged or not.
   * @return std::map containing the contents of the three blocks as a
   * dictionary.
   */
  template < typename _datatype_ >
  inline std::map< std::string, _datatype_ > read_dict(std::ifstream &ifile,
                                                       bool tagged = true) {
    unsigned int size;
    read_block(ifile, size);
    std::vector< std::string > tags(size);
    if (tagged) {
      read_block(ifile, tags);
    } else {
      for (unsigned int i = 0; i < size; ++i) {
        tags[i] = "tag";
      }
    }
    std::vector< _datatype_ > vals(size);
    read_block(ifile, vals);
    std::map< std::string, _datatype_ > dict;
    for (unsigned int i = 0; i < size; ++i) {
      // check for duplicates and add numbers to duplicate tag names
      if (dict.count(tags[i]) == 1) {
        unsigned int count = 1;
        std::stringstream newtag;
        newtag << tags[i] << count;
        while (dict.count(newtag.str()) == 1) {
          ++count;
          newtag.str("");
          newtag << tags[i] << count;
        }
        tags[i] = newtag.str();
      }
      dict[tags[i]] = vals[i];
    }
    return dict;
  }

public:
  SPHNGSnapshotDensityFunction(std::string filename, double initial_temperature,
                               bool write_stats, unsigned int stats_numbin,
                               double stats_mindist, double stats_maxdist,
                               std::string stats_filename, Log *log = nullptr);

  SPHNGSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  ~SPHNGSnapshotDensityFunction();

  virtual void initialize();

  CoordinateVector<> get_position(unsigned int index);
  double get_mass(unsigned int index);
  double get_smoothing_length(unsigned int index);

  virtual DensityValues operator()(const Cell &cell) const;
};

/**
 * @brief Fill the given referenced parameter by reading from the given
 * Fortran unformatted binary file.
 *
 * Template specialization for a std::vector of integers.
 *
 * @param ifile Reference to an open Fortran unformatted binary file.
 * @param value Next (and last) value to read from the file.
 */
template <>
inline void SPHNGSnapshotDensityFunction::read_value< std::vector< int > >(
    std::ifstream &ifile, std::vector< int > &value) {
  ifile.read(reinterpret_cast< char * >(&value[0]), value.size() * sizeof(int));
}

/**
 * @brief Fill the given referenced parameter by reading from the given
 * Fortran unformatted binary file.
 *
 * Template specialization for a std::vector of unsigned integers.
 *
 * @param ifile Reference to an open Fortran unformatted binary file.
 * @param value Next (and last) value to read from the file.
 */
template <>
inline void
SPHNGSnapshotDensityFunction::read_value< std::vector< unsigned int > >(
    std::ifstream &ifile, std::vector< unsigned int > &value) {
  ifile.read(reinterpret_cast< char * >(&value[0]),
             value.size() * sizeof(unsigned int));
}

/**
 * @brief Fill the given referenced parameter by reading from the given
 * Fortran unformatted binary file.
 *
 * Template specialization for a std::vector of unsigned long integers.
 *
 * @param ifile Reference to an open Fortran unformatted binary file.
 * @param value Next (and last) value to read from the file.
 */
template <>
inline void
SPHNGSnapshotDensityFunction::read_value< std::vector< unsigned long > >(
    std::ifstream &ifile, std::vector< unsigned long > &value) {
  ifile.read(reinterpret_cast< char * >(&value[0]),
             value.size() * sizeof(unsigned long));
}

/**
 * @brief Fill the given referenced parameter by reading from the given
 * Fortran unformatted binary file.
 *
 * Template specialization for a std::vector of chars.
 *
 * @param ifile Reference to an open Fortran unformatted binary file.
 * @param value Next (and last) value to read from the file.
 */
template <>
inline void SPHNGSnapshotDensityFunction::read_value< std::vector< char > >(
    std::ifstream &ifile, std::vector< char > &value) {
  ifile.read(&value[0], value.size());
}

/**
 * @brief Fill the given referenced parameter by reading from the given
 * Fortran unformatted binary file.
 *
 * Template specialization for a std::vector of double precision floating point
 * values.
 *
 * @param ifile Reference to an open Fortran unformatted binary file.
 * @param value Next (and last) value to read from the file.
 */
template <>
inline void SPHNGSnapshotDensityFunction::read_value< std::vector< double > >(
    std::ifstream &ifile, std::vector< double > &value) {
  ifile.read(reinterpret_cast< char * >(&value[0]),
             value.size() * sizeof(double));
}

/**
 * @brief Get the size of the given template datatype.
 *
 * Template specialization for a std::vector containing signed integers.
 *
 * @param value Reference to a value of the template datatype.
 * @return Size of the template datatype.
 */
template <>
inline unsigned int
SPHNGSnapshotDensityFunction::get_size< std::vector< int > >(
    std::vector< int > &value) {
  return value.size() * sizeof(int);
}

/**
 * @brief Get the size of the given template datatype.
 *
 * Template specialization for a std::vector containing unsigned integers.
 *
 * @param value Reference to a value of the template datatype.
 * @return Size of the template datatype.
 */
template <>
inline unsigned int
SPHNGSnapshotDensityFunction::get_size< std::vector< unsigned int > >(
    std::vector< unsigned int > &value) {
  return value.size() * sizeof(unsigned int);
}

/**
 * @brief Get the size of the given template datatype.
 *
 * Template specialization for a std::vector containing unsigned long integers.
 *
 * @param value Reference to a value of the template datatype.
 * @return Size of the template datatype.
 */
template <>
inline unsigned int
SPHNGSnapshotDensityFunction::get_size< std::vector< unsigned long > >(
    std::vector< unsigned long > &value) {
  return value.size() * sizeof(unsigned long);
}

/**
 * @brief Get the size of the given template datatype.
 *
 * Template specialization for a std::vector containing signed chars.
 *
 * @param value Reference to a value of the template datatype.
 * @return Size of the template datatype.
 */
template <>
inline unsigned int
SPHNGSnapshotDensityFunction::get_size< std::vector< char > >(
    std::vector< char > &value) {
  return value.size() * sizeof(char);
}

/**
 * @brief Get the size of the given template datatype.
 *
 * Template specialization for a std::vector containing double precision
 * floating point values.
 *
 * @param value Reference to a value of the template datatype.
 * @return Size of the template datatype.
 */
template <>
inline unsigned int
SPHNGSnapshotDensityFunction::get_size< std::vector< double > >(
    std::vector< double > &value) {
  return value.size() * sizeof(double);
}

/**
 * @brief Read a block  from a Fortran unformatted binary file and fill the
 * given referenced template parameters with its contents.
 *
 * Specialization that reads the entire block into a single string.
 *
 * @param ifile Reference to an open Fortran unformatted binary file.
 * @param value Reference to the std::string parameter that should be filled.
 */
template <>
inline void SPHNGSnapshotDensityFunction::read_block(std::ifstream &ifile,
                                                     std::string &value) {
  unsigned int length1, length2;
  ifile.read(reinterpret_cast< char * >(&length1), sizeof(unsigned int));

  char *cstr = new char[length1 + 1];
  ifile.read(cstr, length1);
  unsigned int i = length1;
  // add string termination character
  cstr[i] = '\0';
  // strip trailing whitespace
  while (i > 0 && cstr[i - 1] == ' ') {
    --i;
    cstr[i] = '\0';
  }
  // fill the given variable
  value = std::string(cstr);
  // free the temporary buffer
  delete[] cstr;

  ifile.read(reinterpret_cast< char * >(&length2), sizeof(unsigned int));
  if (length1 != length2) {
    cmac_error("Wrong block size!");
  }
}

/**
 * @brief Read a block  from a Fortran unformatted binary file and fill the
 * given referenced template parameters with its contents.
 *
 * Specialization that reads the entire block into a std::vector of string,
 * assuming a 16 character tag string for each element.
 *
 * If the total size of the block does not match 16 times the size of the given
 * vector, an error is thrown.
 *
 * @param ifile Reference to an open Fortran unformatted binary file.
 * @param value Reference to the std::string parameter that should be filled.
 */
template <>
inline void
SPHNGSnapshotDensityFunction::read_block(std::ifstream &ifile,
                                         std::vector< std::string > &value) {
  unsigned int length1, length2;
  ifile.read(reinterpret_cast< char * >(&length1), sizeof(unsigned int));

  if (length1 % 16 != 0) {
    cmac_error("Block has the wrong size to contain a list of tags!");
  }
  if (value.size() * 16 != length1) {
    cmac_error("Vector of wrong size given!");
  }

  // create a temporary read buffer
  char *cstr = new char[17];
  // add string termination character
  cstr[16] = '\0';
  for (unsigned int i = 0; i < value.size(); ++i) {
    ifile.read(cstr, 16);
    // strip trailing whitespace
    unsigned int j = 16;
    while (j > 0 && cstr[j - 1] == ' ') {
      --j;
      cstr[j] = '\0';
    }
    // fill the given variable
    value[i] = std::string(cstr);
  }
  // free the temporary buffer
  delete[] cstr;

  ifile.read(reinterpret_cast< char * >(&length2), sizeof(unsigned int));
  if (length1 != length2) {
    cmac_error("Wrong block size!");
  }
}

#endif // SPHNGSNAPSHOTDENSITYFUNCTION_HPP
