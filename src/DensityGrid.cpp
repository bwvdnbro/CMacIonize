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
 * @param box Box containing the grid.
 * @param n1D Number of cells in 1 dimensions.
 * @param density_function DensityFunction that defines the density field.
 */
DensityGrid::DensityGrid(Box box, unsigned int n1D,
                         DensityFunction &density_function)
    : _box(box) {
  _n1D = n1D;
  _density = new double **[_n1D];
  for (unsigned int i = 0; i < _n1D; ++i) {
    _density[i] = new double *[_n1D];
    for (unsigned int j = 0; j < _n1D; ++j) {
      _density[i][j] = new double[_n1D];
    }
  }

  // fill the density grid
  double cellside_x = _box.get_sides().x() / _n1D;
  double cellside_y = _box.get_sides().y() / _n1D;
  double cellside_z = _box.get_sides().z() / _n1D;
  for (unsigned int i = 0; i < _n1D; ++i) {
    for (unsigned int j = 0; j < _n1D; ++j) {
      for (unsigned int k = 0; k < _n1D; ++k) {
        double x = _box.get_anchor().x() + (i + 0.5) * cellside_x;
        double y = _box.get_anchor().y() + (j + 0.5) * cellside_y;
        double z = _box.get_anchor().z() + (k + 0.5) * cellside_z;
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
  double cellside_x = _box.get_sides().x() / _n1D;
  double cellside_y = _box.get_sides().y() / _n1D;
  double cellside_z = _box.get_sides().z() / _n1D;
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

/**
 * @brief Get the distance that needs to be travelled by a photon with the given
 * origin in the given direction to cover the given optical depth.
 *
 * @param photon_origin Origin of the photon.
 * @param photon_direction Direction the photon is moving in.
 * @param optical_depth Optical depth the photon should travel in total.
 * @return Total distance covered by the photon.
 */
double DensityGrid::get_distance(CoordinateVector photon_origin,
                                 CoordinateVector photon_direction,
                                 double optical_depth) {
  double S = 0.;

  // we need these
  double cellside_x = _box.get_sides().x() / _n1D;
  double cellside_y = _box.get_sides().y() / _n1D;
  double cellside_z = _box.get_sides().z() / _n1D;

  double cellside_max = cellside_x;
  if (cellside_y > cellside_max) {
    cellside_max = cellside_y;
  }
  if (cellside_z > cellside_max) {
    cellside_max = cellside_z;
  }

  // while the photon has not exceeded the optical depth and is still in the box
  while (optical_depth > 0.) {

    // find out in which cell the photon is currently hiding
    unsigned int ix = (photon_origin.x() - _box.get_anchor().x()) / cellside_x;
    unsigned int iy = (photon_origin.y() - _box.get_anchor().y()) / cellside_y;
    unsigned int iz = (photon_origin.z() - _box.get_anchor().z()) / cellside_z;

    double cell_xmin = _box.get_anchor().x() + cellside_x * ix;
    double cell_xmax = _box.get_anchor().x() + cellside_x * (ix + 1);
    double cell_ymin = _box.get_anchor().y() + cellside_y * iy;
    double cell_ymax = _box.get_anchor().y() + cellside_y * (iy + 1);
    double cell_zmin = _box.get_anchor().z() + cellside_z * iz;
    double cell_zmax = _box.get_anchor().z() + cellside_z * (iz + 1);

    // find out which cell wall the photon is going to hit next
    CoordinateVector next_x;
    double l;
    if (photon_direction.x() > 0.) {
      // if the photon starts at \vec{o} and travels in the direction \vec{d},
      // the general position of the photon at any later time is given by
      // \vec{o} + l*\vec{d}, with l some positive parameter
      // we know that the photon hits the next x wall when the x-component of
      // this expression equals cell_xmax, so we can solve for l:
      l = (cell_xmax - photon_origin.x()) / photon_direction.x();
    } else if (photon_direction.x() < 0.) {
      l = (cell_xmin - photon_origin.x()) / photon_direction.x();
    } else {
      // we never reach an x wall, since the photon travels parallel with it
      // we just set l to a ridiculous value that will always cause dx
      // to be larger than dy and/or dz
      // we know that at least one of the direction components needs to be
      // larger
      // than 1/\sqrt{3}, since the minimal direction components are found when
      // all three are equal, and we have 3*component_size^2 = 1 (due to the
      // normalization).
      // the largest l values are found for the smallest components, so this
      // gives
      // us a good bound on the denominator in the expression for l
      // the numerator is bound by the cell size: the maximal value is obtained
      // for a cell size cellside_max
      // in other words: if we set l > \sqrt{3}*cellside_max, we are sure dy or
      // dz will always be smaller
      l = 1000. * cellside_max;
    }
    // the y and z coordinates are then trivially found
    next_x = photon_origin + l * photon_direction;
    double dx = (next_x - photon_origin).norm2();

    CoordinateVector next_y;
    if (photon_direction.y() > 0.) {
      l = (cell_ymax - photon_origin.y()) / photon_direction.y();
    } else if (photon_direction.y() < 0.) {
      l = (cell_ymin - photon_origin.y()) / photon_direction.y();
    } else {
      l = 1000. * cellside_max;
    }
    next_y = photon_origin + l * photon_direction;
    double dy = (next_y - photon_origin).norm2();

    CoordinateVector next_z;
    if (photon_direction.z() > 0.) {
      l = (cell_zmax - photon_origin.z()) / photon_direction.z();
    } else if (photon_direction.z() < 0.) {
      l = (cell_zmin - photon_origin.z()) / photon_direction.z();
    } else {
      l = 1000. * cellside_max;
    }
    next_z = photon_origin + l * photon_direction;
    double dz = (next_z - photon_origin).norm2();

    CoordinateVector next_wall;
    double ds;
    // it the equality holds, the photon ends up on a corner and we don't really
    // care which wall is chosen
    if (dx <= dy && dx <= dz) {
      next_wall = next_x;
      ds = dx;
    } else if (dy <= dx && dy <= dz) {
      next_wall = next_y;
      ds = dy;
    } else {
      next_wall = next_z;
      ds = dz;
    }

    // get the optical depth of the path from the current photon location to the
    // cell wall, update S
    S += ds;
    optical_depth -= ds;

    // if the optical depth exceeds or equals the wanted value: exit the loop

    // if the optical depth exceeded the wanted value: find out where in the
    // cell
    // we end up, and correct S
    if (optical_depth < 0.) {
      S += optical_depth;
    }
  }

  return S;
}
