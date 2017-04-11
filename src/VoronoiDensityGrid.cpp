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
 * @file VoronoiDensityGrid.cpp
 *
 * @brief VoronoiDensityGrid implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "VoronoiDensityGrid.hpp"

/**
 * @brief Constructor.
 *
 * @param density_function DensityFunction to use to initialize the cell
 * variables.
 * @param box Box containing the entire grid (in m).
 * @param periodic Periodicity flags.
 * @param hydro Flag signaling if hydro is active or not.
 * @param log Log to write logging info to.
 */
VoronoiDensityGrid::VoronoiDensityGrid(DensityFunction &density_function,
                                       Box box,
                                       CoordinateVector< bool > periodic,
                                       bool hydro, Log *log)
    : DensityGrid(density_function, box, periodic, hydro, log),
      _voronoi_grid(box, periodic) {}

/**
 * @brief Get the number of cells in the grid.
 *
 * @return Number of cells in the grid.
 */
unsigned int VoronoiDensityGrid::get_number_of_cells() const { return 0; }

/**
 * @brief Get the index of the cell that contains the given position.
 *
 * @param position Position (in m).
 * @return Index of the cell that contains the position.
 */
unsigned long
VoronoiDensityGrid::get_cell_index(CoordinateVector<> position) const {
  return 0ull;
}

/**
 * @brief Get the midpoint of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Position of the midpoint of that cell (in m).
 */
CoordinateVector<>
VoronoiDensityGrid::get_cell_midpoint(unsigned long index) const {
  return CoordinateVector<>(0.);
}

/**
 * @brief Get the neighbours of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return std::vector containing the neighbours of the cell.
 */
std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                         CoordinateVector<>, double > >
VoronoiDensityGrid::get_neighbours(unsigned long index) {
  std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                           CoordinateVector<>, double > >
      ngbs;
  return ngbs;
}

/**
 * @brief Get the volume of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Volume of that cell (in m^3).
 */
double VoronoiDensityGrid::get_cell_volume(unsigned long index) const {
  return 0.;
}

/**
 * @brief Traverse the given Photon through the grid until the given optical
 * depth is reached (or the Photon leaves the system).
 *
 * @param photon Photon to use.
 * @param optical_depth Target optical depth.
 * @return DensityGrid::iterator to the cell that contains the Photon when it
 * reaches the target optical depth, or VoronoiDensityGrid::end if the Photon
 * leaves the system.
 */
DensityGrid::iterator VoronoiDensityGrid::interact(Photon &photon,
                                                   double optical_depth) {
  return DensityGrid::iterator(0, *this);
}

/**
 * @brief Get an iterator to the first cell in the grid.
 *
 * @return DensityGrid::iterator to the first cell in the grid.
 */
DensityGrid::iterator VoronoiDensityGrid::begin() {
  return DensityGrid::iterator(0, *this);
}

/**
 * @brief Get an iterator to the beyond last cell in the grid.
 *
 * @return DensityGrid::iterator to the beyond last cell in the grid.
 */
DensityGrid::iterator VoronoiDensityGrid::end() {
  return DensityGrid::iterator(0, *this);
}
