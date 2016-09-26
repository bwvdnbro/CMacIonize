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
 * @file DensityGrid.cpp
 *
 * @brief Density grid: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DensityGrid.hpp"
#include "DensityFunction.hpp"
using namespace std;

/**
 * @brief Constructor
 *
 * @param anchor_x x coordinate of the anchor of the box containing the grid.
 * @param anchor_y y coordinate of the anchor of the box containing the grid.
 * @param anchor_z z coordinate of the anchor of the box containing the grid.
 * @param side_x Side lenght of the box containing the grid in the x direction.
 * @param side_y Side lenght of the box containing the grid in the y direction.
 * @param side_z Side lenght of the box containing the grid in the z direction.
 * @param n1D Number of cells in 1 dimensions.
 * @param density_function DensityFunction that defines the density field.
 */
DensityGrid::DensityGrid(double anchor_x, double anchor_y, double anchor_z,
                         double side_x, double side_y, double side_z,
                         unsigned int n1D, DensityFunction &density_function) {
  _anchor[0] = anchor_x;
  _anchor[1] = anchor_y;
  _anchor[2] = anchor_z;
  _side[0] = side_x;
  _side[1] = side_y;
  _side[2] = side_z;

  _n1D = n1D;
  _density = new double **[_n1D];
  for (unsigned int i = 0; i < _n1D; ++i) {
    _density[i] = new double *[_n1D];
    for (unsigned int j = 0; j < _n1D; ++j) {
      _density[i][j] = new double[_n1D];
    }
  }

  // fill the density grid
  double cellside_x = _side[0] / _n1D;
  double cellside_y = _side[1] / _n1D;
  double cellside_z = _side[2] / _n1D;
  for (unsigned int i = 0; i < _n1D; ++i) {
    for (unsigned int j = 0; j < _n1D; ++j) {
      for (unsigned int k = 0; k < _n1D; ++k) {
        double x = _anchor[0] + (i + 0.5) * cellside_x;
        double y = _anchor[1] + (j + 0.5) * cellside_y;
        double z = _anchor[2] + (k + 0.5) * cellside_z;
        _density[i][j][k] = density_function(CoordinateVector(x, y, z));
      }
    }
  }
}

/**
 * @brief Destructor
 *
 * Free the memory used by the internal arrays.
 */
DensityGrid::~DensityGrid() {
  for (unsigned int i = 0; i < _n1D; ++i) {
    for (unsigned int j = 0; j < _n1D; ++j) {
      delete[] _density[i][j];
    }
    delete[] _density[i];
  }
  delete[] _density;
}

/**
 * @brief Get the total mass contained in the grid.
 *
 * @return Total mass contained in the grid.
 */
double DensityGrid::get_total_mass() {
  double mtot = 0;
  double cellside_x = _side[0] / _n1D;
  double cellside_y = _side[1] / _n1D;
  double cellside_z = _side[2] / _n1D;
  double cellvolume = cellside_x * cellside_y * cellside_z;

  for (unsigned int i = 0; i < _n1D; ++i) {
    for (unsigned int j = 0; j < _n1D; ++j) {
      for (unsigned int k = 0; k < _n1D; ++k) {
        mtot += _density[i][j][k] * cellvolume;
      }
    }
  }

  return mtot;
}
