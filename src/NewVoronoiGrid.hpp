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
 * @file NewVoronoiGrid.hpp
 *
 * @brief Voronoi grid implementation that uses an incremental Delaunay
 * construction algorithm (that should work in all cases, even in highly
 * degenerate grids).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEWVORONOIGRID_HPP
#define NEWVORONOIGRID_HPP

#include "NewVoronoiCell.hpp"

#include <vector>

/**
 * @brief Voronoi grid implementation that uses an incremental Delaunay
 * construction algorithm (that should work in all cases, even in highly
 * degenerate grids).
 */
class NewVoronoiGrid {
private:
  /*! @brief Simulation box (in m). */
  const Box<> _box;

  /*! @brief Reference to the mesh generating positions (in m). */
  const std::vector< CoordinateVector<> > &_real_generator_positions;

  /*! @brief Real VoronoiBox (in m). */
  const VoronoiBox< double > _real_voronoi_box;

  /*! @brief Integer representation of the mesh generating positions. */
  std::vector< CoordinateVector< unsigned long > > _integer_generator_positions;

  /*! @brief Integer VoronoiBox. */
  VoronoiBox< unsigned long > _integer_voronoi_box;

  /*! @brief Voronoi cells. */
  std::vector< NewVoronoiCell > _cells;

public:
  NewVoronoiGrid(const std::vector< CoordinateVector<> > &positions,
                 const Box<> box);

  void construct();

  /**
   * @brief Auxiliary typedef used to extract the mantissa from a double
   * precision floating point value.
   *
   * The variables in the union occupy the same memory, which allows us to
   * access the bytes used to store the double precision floating point value.
   */
  typedef union {
    /*! @brief Double precision floating point value. */
    double dvalue;
    /**
     * @brief Anonymous struct containing the 3 parts of a general double
     * precision floating point value.
     *
     * A general double precision floating point value has a sign \f$s\f$, a
     * mantissa \f$m\f$, and an exponent \f$e\f$, so that the value \f$v\f$ of
     * the double is given by
     * \f[
     * v = s \times 1.m \times 2^{e-1023}
     * \f]
     * (see
     * https://en.wikipedia.org/wiki/Double-precision_floating-point_format).
     *
     * The IEEE 754 standard specifies that \f$s\f$, \f$m\f$ and \f$e\f$ have
     * respectively 1 bit, 52 bit and 11 bit precision. In memory (low to high
     * bits), they are ordered as follows:
     * \f[
     * m e s
     * \f]
     */
    struct {
      /*! @brief Mantissa \f$m\f$. */
      unsigned long mantissa : 52;
      /*! @brief Exponent \f$e\f$. */
      unsigned long exponent : 11;
      /*! @brief Sign \f$s\f$. */
      unsigned long sign : 1;
    } parts;
  } binary_double;

  /**
   * @brief Get the 52 bit mantissa of the given double precision floating point
   * value.
   *
   * @param value Double precision floating point value.
   * @return 52 bit mantissa of that same value.
   */
  static inline unsigned long get_mantissa(double value) {
    binary_double dvalue;
    dvalue.dvalue = value;
    return dvalue.parts.mantissa;
  }
};

#endif // NEWVORONOIGRID_HPP
