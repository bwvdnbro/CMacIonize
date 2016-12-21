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

#include "DensityFunction.hpp"

#include <fstream>

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

  /*! @brief Octree used to speed up neighbour finding. */
  Octree *_octree;

  static double kernel(const double q, const double h);

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
    read_value(ifile, args...);
    ifile.read(reinterpret_cast< char * >(&length2), sizeof(unsigned int));
    if (length1 != length2) {
      cmac_error("Wrong block size!");
    }
  }

public:
  SPHNGSnapshotDensityFunction(std::string filename, Log *log = nullptr);

  SPHNGSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  ~SPHNGSnapshotDensityFunction();

  CoordinateVector<> get_position(unsigned int index);
  double get_mass(unsigned int index);
  double get_smoothing_length(unsigned int index);

  DensityValues operator()(CoordinateVector<> position) const;
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

#endif // SPHNGSNAPSHOTDENSITYFUNCTION_HPP
