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
 */
class InterpolatedDensityFunction : public DensityFunction {
private:
  /*! @brief X-coordinates to interpolate on (in m). */
  std::vector< double > _x_coords;

  /*! @brief Y-coordinates to interpolate on (in m). */
  std::vector< double > _y_coords;

  /*! @brief Z-coordinates to interpolate on (in m). */
  std::vector< double > _z_coords;

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
