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
 * @file ParallelCartesianDensityGrid.hpp
 *
 * @brief DensityGrid (implementation?) containing a distributed Cartesian grid.
 *
 * This grid is a prototype for DensityGrids that can be distributed across
 * multiple processes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PARALLELCARTESIANDENSITYGRID_HPP
#define PARALLELCARTESIANDENSITYGRID_HPP

#include "Box.hpp"
#include "ParallelCartesianDensitySubGrid.hpp"
#include "Utilities.hpp"

#include <vector>

/**
 * @brief DensityGrid (implementation?) containing a distributed Cartesian grid.
 */
class ParallelCartesianDensityGrid {
private:
  /*! @brief ParallelCartesianDensitySubGrids that make up the actual underlying
   *  grid. */
  std::vector< ParallelCartesianDensitySubGrid * > _subgrids;

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the entire grid.
   * @param numcell Number of cells in each dimension.
   * @param numdomain Total number of subgrids.
   */
  ParallelCartesianDensityGrid(Box box, CoordinateVector< int > numcell,
                               unsigned int numdomain) {
    CoordinateVector< int > blockresolution =
        Utilities::subdivide(numcell, numdomain);
    CoordinateVector<> blockside = box.get_sides();
    blockside[0] /= blockresolution[0];
    blockside[1] /= blockresolution[1];
    blockside[2] /= blockresolution[2];
    (void)blockside;
  }
};

#endif // PARALLELCARTESIANDENSITYGRID_HPP
