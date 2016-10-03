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
 * @file DensityGrid.hpp
 *
 * @brief Density grid: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRID_HPP
#define DENSITYGRID_HPP

#include "Box.hpp"

class DensityFunction;

class DensityGrid {
private:
  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Side lengths of a single cell. */
  CoordinateVector _cellside;

  /*! @brief Maximal cell side among the three dimensions. */
  double _cellside_max;

  /*! @brief Number of cells in 1 dimension. */
  unsigned int _n1D;

  /*! @brief Density grid. */
  double ***_density;

public:
  DensityGrid(Box box, unsigned int n1D, DensityFunction &density_function);

  ~DensityGrid();

  double get_total_mass();

  void get_cell_indices(CoordinateVector position, unsigned int &ix,
                        unsigned int &iy, unsigned int &iz);
  Box get_cell(unsigned int ix, unsigned int iy, unsigned int iz);
  CoordinateVector get_wall_intersection(CoordinateVector &photon_origin,
                                         CoordinateVector &photon_direction,
                                         Box &cell, char &ix, char &iy,
                                         char &iz, double &ds);

  double get_distance(CoordinateVector photon_origin,
                      CoordinateVector photon_direction, double optical_depth);
};

#endif // DENSITYGRID_HPP
