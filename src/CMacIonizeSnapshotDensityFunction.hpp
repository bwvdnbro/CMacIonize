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
 * @file CMacIonizeSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density grid from a CMacIonize snapshot.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CMACIONIZESNAPSHOTDENSITYFUNCTION_HPP
#define CMACIONIZESNAPSHOTDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads a density grid from a CMacIonize snapshot.
 */
class CMacIonizeSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Number of cells in each dimension. */
  CoordinateVector< int > _ncell;

  /*! @brief Density grid. */
  double ***_grid;

public:
  CMacIonizeSnapshotDensityFunction(std::string filename, Log *log = nullptr);

  CMacIonizeSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual ~CMacIonizeSnapshotDensityFunction();

  virtual double operator()(CoordinateVector<> position);
};

#endif // CMACIONIZESNAPSHOTDENSITYFUNCTION_HPP
