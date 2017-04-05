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
 * @file DensityMask.hpp
 *
 * @brief General interface for masks that can be applied to an existing density
 * grid to alter the structure of the density field.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYMASK_HPP
#define DENSITYMASK_HPP

class DensityGrid;

/**
 * @brief General interface for masks that can be applied to an existing density
 * grid to alter the structure of the density field.
 */
class DensityMask {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~DensityMask() {}

  /**
   * @brief Apply the mask to the given DensityGrid.
   *
   * Note that this operation should not alter the total number of hydrogen
   * atoms in the DensityGrid.
   *
   * @param grid DensityGrid to apply the mask to.
   */
  virtual void apply(DensityGrid &grid) const = 0;
};

#endif // DENSITYMASK_HPP
