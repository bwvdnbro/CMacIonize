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
 * @file AMRRefinementScheme.hpp
 *
 * @brief General interface for schemes used to refine an AMRDensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AMRREFINEMENTSCHEME_HPP
#define AMRREFINEMENTSCHEME_HPP

#include "CoordinateVector.hpp"

class DensityValues;

/**
 * @brief General interface for schemes used to refine an AMRDensityGrid.
 *
 * We provide empty implementations for all routines, although in practice
 * every implementation should implement the refine() method, as otherwise the
 * implementation is completely useless.
 */
class AMRRefinementScheme {
public:
  virtual ~AMRRefinementScheme() {}

  /**
   * @brief Decide if the given cell should be refined or not.
   *
   * @param level Current refinement level of the cell.
   * @param midpoint Midpoint of the cell (in m).
   * @param volume Volume of the cell (in m^3).
   * @param cell DensityValues of a cell.
   * @return True if the cell should be refined.
   */
  virtual bool refine(unsigned char level, CoordinateVector<> midpoint,
                      double volume, DensityValues &cell) const {
    return false;
  }
};

#endif // AMRREFINEMENTSCHEME_HPP
