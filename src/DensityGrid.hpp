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

class DensityFunction;

class DensityGrid {
private:
  /*! @brief Bottom front left corner of the box containing the grid. */
  double _anchor[3];

  /*! @brief Side lengths of the box containing the grid. */
  double _side[3];

  /*! @brief Number of cells in 1 dimension. */
  unsigned int _n1D;

  /*! @brief Density grid. */
  double ***_density;

public:
  DensityGrid(double anchor_x, double anchor_y, double anchor_z, double side_x,
              double side_y, double side_z, unsigned int n1D,
              DensityFunction &density_function);

  ~DensityGrid();

  double get_total_mass();
};

#endif // DENSITYGRID_HPP
