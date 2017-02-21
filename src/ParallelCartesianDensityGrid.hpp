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
  std::vector< DensitySubGrid * > _subgrids;

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the entire grid.
   * @param numcell Number of cells in each dimension.
   * @param numdomain Total number of subgrids.
   * @param domain Part of the domain stored on the local process.
   */
  ParallelCartesianDensityGrid(Box box, CoordinateVector< int > numcell,
                               unsigned int numdomain,
                               std::pair< int, int > domain) {
    CoordinateVector< int > blockresolution =
        Utilities::subdivide(numcell, numdomain);
    CoordinateVector<> blocksides = box.get_sides();
    blocksides[0] /= blockresolution[0];
    blocksides[1] /= blockresolution[1];
    blocksides[2] /= blockresolution[2];

    _subgrids.reserve(blockresolution.x() * blockresolution.y() *
                      blockresolution.z());
    for (int ix = 0; ix < blockresolution.x(); ++ix) {
      for (int iy = 0; iy < blockresolution.y(); ++iy) {
        for (int iz = 0; iz < blockresolution.z(); ++iz) {
          int index = ix * blockresolution.y() * blockresolution.z() +
                      iy * blockresolution.z() + iz;
          if (index >= domain.first && index < domain.second) {
            CoordinateVector<> blockanchor = box.get_anchor();
            blockanchor[0] += ix * blocksides[0];
            blockanchor[1] += iy * blocksides[1];
            blockanchor[2] += iz * blocksides[2];
            Box blockbox(blockanchor, blocksides);
            _subgrids.push_back(
                new ParallelCartesianDensitySubGrid(blockbox, blockresolution));
          } else {
            _subgrids.push_back(new GhostDensitySubGrid());
          }
        }
      }
    }

    // set neighbour relations
    for (int ix = 0; ix < blockresolution.x(); ++ix) {
      for (int iy = 0; iy < blockresolution.y(); ++iy) {
        for (int iz = 0; iz < blockresolution.z(); ++iz) {
          int index = ix * blockresolution.y() * blockresolution.z() +
                      iy * blockresolution.z() + iz;
          if (index >= domain.first && index < domain.second) {
            ParallelCartesianDensitySubGrid *block =
                reinterpret_cast< ParallelCartesianDensitySubGrid * >(
                    _subgrids[index]);
            if (ix > 0) {
              int ngbindex =
                  (ix - 1) * blockresolution.y() * blockresolution.z() +
                  iy * blockresolution.z() + iz;
              block->set_neighbour(0, ngbindex);
            }
            if (ix < blockresolution.x() - 1) {
              int ngbindex =
                  (ix + 1) * blockresolution.y() * blockresolution.z() +
                  iy * blockresolution.z() + iz;
              block->set_neighbour(1, ngbindex);
            }

            if (iy > 0) {
              int ngbindex = ix * blockresolution.y() * blockresolution.z() +
                             (iy - 1) * blockresolution.z() + iz;
              block->set_neighbour(2, ngbindex);
            }
            if (iy < blockresolution.y() - 1) {
              int ngbindex = ix * blockresolution.y() * blockresolution.z() +
                             (iy + 1) * blockresolution.z() + iz;
              block->set_neighbour(3, ngbindex);
            }

            if (iz > 0) {
              int ngbindex = ix * blockresolution.y() * blockresolution.z() +
                             iy * blockresolution.z() + iz - 1;
              block->set_neighbour(4, ngbindex);
            }
            if (iz < blockresolution.z() - 1) {
              int ngbindex = ix * blockresolution.y() * blockresolution.z() +
                             iy * blockresolution.z() + iz + 1;
              block->set_neighbour(5, ngbindex);
            }
          }
        }
      }
    }
  }

  /**
   * @brief Destructor.
   *
   * Free memory occupied by the blocks.
   */
  ~ParallelCartesianDensityGrid() {
    for (unsigned int i = 0; i < _subgrids.size(); ++i) {
      delete _subgrids[i];
    }
  }
};

#endif // PARALLELCARTESIANDENSITYGRID_HPP
