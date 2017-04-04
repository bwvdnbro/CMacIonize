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
 * @file InterpolatedDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density field from a 1D, 2D or 3D data
 * cube in a text file and interpolates on it.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef INTERPOLATEDDENSITYFUNCTION_HPP
#define INTERPOLATEDDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads a density field from a 1D, 2D or 3D data
 * cube in a text file and interpolates on it.
 *
 * The file should consist of two parts: a header that contains information
 * about the data in the file, and a data block containing the actual data.
 * The data block should only contain numerical values, organized in rows with
 * an equal number of columns.
 *
 * The header block uses YAML-syntax (the same syntax used in the ParameterFile)
 * and is contained within a starting and an ending line that contain three
 * dashes.
 * The header block should look like this (capital words denote values of the
 * type indicated in their name):
 * @verbatim
   ---
   num_x: POSITIVE_OR_ZERO_INTEGER
   xmin: FLOAT_WITH_LENGTH_UNIT
   xmax: FLOAT_WITH_LENGTH_UNIT
   num_y: POSITIVE_OR_ZERO_INTEGER
   ymin: FLOAT_WITH_LENGTH_UNIT
   ymax: FLOAT_WITH_LENGTH_UNIT
   num_z: POSITIVE_OR_ZERO_INTEGER
   zmin: FLOAT_WITH_LENGTH_UNIT
   zmax: FLOAT_WITH_LENGTH_UNIT
   num_column: POSITIVE_INTEGER
   column_X_variable: STRING
   column_X_unit: STRING_UNIT_NAME
   ---
   @endverbatim
 * @c X represents a counter, with @c X taking values from 0 to @c num_column-1.
 * All these fields need to be present. @c num_column also denotes the number of
 * columns present in the data block; the fields @c column_X_variable and
 * @c column_X_unit are respectively the name and unit of the variable stored in
 * that specific column. Valid names and units are: @c x (length unit), @c y
 * (length unit), @c z (length unit), and @c number @c density (number density
 * unit). A column containing the @c number @c density is required.
 *
 * @c num_x, @c num_y, and @c num_z denote the number of @e different values for
 * these respective coordinates that are present in the data block. If set to a
 * non-zero value, the corresponding column needs to be present in the data
 * block. If set to zero, we simply assume the density is constant in that
 * coordinate direction, and the corresponding column in the data block is
 * ignored (if present).
 *
 * The total number of rows in the data block needs to be equal to @c num_x @c *
 * @c num_y @c * @c num_z (dimensions for which @c num_D is zero are ignored in
 * this expression), and we expect the different coordinate columns (if present)
 * to traverse all possible permutations of the available coordinate values. The
 * order in which the different coordinate values increase is detected
 * aumatically, but we do expect it to be consistent (and they do need to
 * @e increase). Likewise, the order of the columns in the data block does not
 * matter, as long as it is consistent with the variable names and units given
 * in the header block.
 *
 * The data block should follow immediately after the header, and no extra lines
 * should be present in the file below the data block. However, all lines above
 * the header are ignored, as are YAML-syntax comments inside the header.
 */
class InterpolatedDensityFunction : public DensityFunction {
private:
  /*! @brief X-coordinates to interpolate on (in m). */
  std::vector< double > _x_coords;

  /*! @brief Minimal and maximal X-coordinate allowed (in m). */
  std::pair< double, double > _x_bounds;

  /*! @brief Y-coordinates to interpolate on (in m). */
  std::vector< double > _y_coords;

  /*! @brief Minimal and maximal Y-coordinate allowed (in m). */
  std::pair< double, double > _y_bounds;

  /*! @brief Z-coordinates to interpolate on (in m). */
  std::vector< double > _z_coords;

  /*! @brief Minimal and maximal Z-coordinate allowed (in m). */
  std::pair< double, double > _z_bounds;

  /*! @brief Number density values to interpolate on (in m^-3). */
  std::vector< std::vector< std::vector< double > > > _number_densities;

public:
  InterpolatedDensityFunction(std::string filename, Log *log = nullptr);

  InterpolatedDensityFunction(ParameterFile &params, Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~InterpolatedDensityFunction() {}

  virtual DensityValues operator()(CoordinateVector<> position) const;
};

#endif // INTERPOLATEDDENSITYFUNCTION_HPP
