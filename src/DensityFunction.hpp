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
 * @file DensityFunction.hpp
 *
 * @brief Interface for functors that can be used to fill a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYFUNCTION_HPP
#define DENSITYFUNCTION_HPP

#include "Cell.hpp"
#include "DensityValues.hpp"

/**
 * @brief Interface for functors that can be used to fill a DensityGrid.
 */
class DensityFunction {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~DensityFunction() {}

  /**
   * @brief Perform all computationally expensive initialization that needs to
   * be done before operator() will work.
   *
   * This routine does not need to be implemented by all implementations. It is
   * particularly useful for implementations that read SPH snapshots, since
   * building the Octree used for SPH density sampling can be demanding.
   */
  virtual void initialize() {}

  /**
   * @brief Free up the memory used by the density function. After this,
   * operator() will no longer work.
   */
  virtual void free() {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) = 0;
};

#endif // DENSITYFUNCTION_HPP
