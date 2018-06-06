/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file HydroMask.hpp
 *
 * @brief General interface for hydro masks.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROMASK_HPP
#define HYDROMASK_HPP

#include "DensityGrid.hpp"

/**
 * @brief Masked out region where the hydrodynamics is artificially reset to a
 * constant value.
 */
class HydroMask {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~HydroMask() {}

  /**
   * @brief Initialize the mask before the first hydrodynamical time step.
   *
   * @param grid DensityGrid to read from.
   */
  virtual void initialize_mask(DensityGrid &grid) {}

  /**
   * @brief Apply the mask to the given DensityGrid.
   *
   * The primitive and conserved variables of all cells within the mask will be
   * updated, all other cells are left untouched.
   *
   * @param grid DensityGrid to update.
   */
  virtual void apply_mask(DensityGrid &grid) const = 0;
};

#endif // HYDROMASK_HPP
