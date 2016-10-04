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

class CrossSections;
class DensityFunction;
class DensityValues;

class DensityGrid {
private:
  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Side lengths of a single cell. */
  CoordinateVector<> _cellside;

  /*! @brief Maximal cell side among the three dimensions. */
  double _cellside_max;

  /*! @brief Number of cells per dimension. Note that by choosing an unsigned
   *  char type, we automatically limit grid sizes to 256^3. */
  CoordinateVector< unsigned char > _ncell;

  /*! @brief Density grid. */
  DensityValues ***_density;

  /*! @brief Cross sections for photoionization. */
  CrossSections &_cross_sections;

public:
  DensityGrid(Box box, CoordinateVector< unsigned char > ncell,
              DensityFunction &density_function, CrossSections &cross_sections);

  ~DensityGrid();

  double get_total_mass();

  CoordinateVector< int > get_cell_indices(CoordinateVector<> position);
  Box get_cell(CoordinateVector< int > index);
  bool is_inside(CoordinateVector< int > index);
  CoordinateVector<> get_wall_intersection(CoordinateVector<> &photon_origin,
                                           CoordinateVector<> &photon_direction,
                                           Box &cell,
                                           CoordinateVector< char > &next_index,
                                           double &ds);

  double get_distance(CoordinateVector<> photon_origin,
                      CoordinateVector<> photon_direction,
                      double optical_depth);
};

#endif // DENSITYGRID_HPP
