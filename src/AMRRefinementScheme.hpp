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
#include "DensityGrid.hpp"

/**
 * @brief General interface for schemes used to refine an AMRDensityGrid.
 *
 * We provide empty implementations for all routines, although in practice
 * every implementation should implement the refine() method, as otherwise the
 * implementation is completely useless.
 */
class AMRRefinementScheme {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~AMRRefinementScheme() {}

  /**
   * @brief Decide if the given cell should be refine or not.
   *
   * @param level Current refinement level of the cell.
   * @param cell DensityGrid::iterator pointing to a cell.
   * @return True if the cell should be refined.
   */
  virtual bool refine(unsigned char level, DensityGrid::iterator &cell) const {
    return false;
  }

  /**
   * @brief Decide if the given cells should be replaced by a single cell or
   * not.
   *
   * @param level Current refinement level of the cells.
   * @param cells DensityGrid::iterators pointing to the cells.
   * @return True if the cells can be replaced by a single cell on a coarser
   * level.
   */
  virtual bool coarsen(unsigned char level,
                       const DensityGrid::iterator *cells) const {
    return false;
  }
};

#endif // AMRREFINEMENTSCHEME_HPP
