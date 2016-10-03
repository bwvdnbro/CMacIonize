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
  _cellside = CoordinateVector<>(cellside_x, cellside_y, cellside_z);
  for (unsigned int i = 0; i < _n1D; ++i) {
    for (unsigned int j = 0; j < _n1D; ++j) {
      for (unsigned int k = 0; k < _n1D; ++k) {
        double x = _box.get_anchor().x() + (i + 0.5) * _cellside.x();
        double y = _box.get_anchor().y() + (j + 0.5) * _cellside.y();
        double z = _box.get_anchor().z() + (k + 0.5) * _cellside.z();
        _density[i][j][k] = density_function(CoordinateVector<>(x, y, z));
      }
    }
  }

  _cellside_max = _cellside.x();
  if (_cellside.y() > _cellside_max) {
    _cellside_max = _cellside.y();
  }
  if (_cellside.z() > _cellside_max) {
    _cellside_max = _cellside.z();
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
  double cellvolume = _cellside.x() * _cellside.y() * _cellside.z();

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
 * @brief Get the indices of the cell containing the given coordinates.
 *
 * @param position CoordinateVector containing coordinates we want to locate.
 * @param ix Variable to store the x index of the cell in.
 * @param iy Variable to store the y index of the cell in.
 * @param iz Variable to store the z index of the cell in.
 */
void DensityGrid::get_cell_indices(CoordinateVector<> position,
                                   unsigned int &ix, unsigned int &iy,
                                   unsigned int &iz) {
  ix = (position.x() - _box.get_anchor().x()) / _cellside.x();
  iy = (position.y() - _box.get_anchor().y()) / _cellside.y();
  iz = (position.z() - _box.get_anchor().z()) / _cellside.z();
}

/**
 * @brief Get the geometrical box of the cell with the given indices.
 *
 * Since we care more about the bottom front left corner and the upper back
 * right corner, the Box we return does not use a side and anchor, but rather
 * two anchors.
 *
 * @param ix x index of the cell.
 * @param iy y index of the cell.
 * @param iz z index of the cell.
 * @return Box containing the bottom front left corner and the upper back right
 * corner of the cell.
 */
Box DensityGrid::get_cell(unsigned int ix, unsigned int iy, unsigned int iz) {
  double cell_xmin = _box.get_anchor().x() + _cellside.x() * ix;
  double cell_ymin = _box.get_anchor().y() + _cellside.y() * iy;
  double cell_zmin = _box.get_anchor().z() + _cellside.z() * iz;
  return Box(CoordinateVector<>(cell_xmin, cell_ymin, cell_zmin), _cellside);
}

/**
 * @brief Get the intersection point of a photon with one of the walls of a
 * cell.
 *
 * We assume the photon is in the cell and find the closest intersection point
 * with one of the walls of the cell (in the travel direction). We also set
 * appropriate indices to find the neighbouring cell on the other side of the
 * wall.
 *
 * The assumption that the photon is inside the cell is not strict: if it is
 * not, we will still find the intersection point with the closest wall on the
 * photon travel line, but this point could possibly be in a direction opposite
 * to the movement direction of the photon. We exploit this to handle round off:
 * it is possible that a photon ends up very close to an edge or corner of a
 * cell, and hence very close to a wall in the neighbouring cell. Due to round
 * off error, the photon might actually lie on the opposite side of the wall in
 * that neighbouring cell. However, our algorithm will still find the correct
 * wall and will return the correct neighbour on the other side of that wall.
 * We will also find a very small distance covered in the wrong direction in our
 * cell, but this distance is negligible.
 *
 * @param photon_origin Current position of the photon.
 * @param photon_direction Direction the photon is travelling in.
 * @param cell Cell in which the photon currently resides.
 * @param ix x index of the neighbouring cell relative to the current cell.
 * @param iy y index of the neighbouring cell relative to the current cell.
 * @param iz z index of the neighbouring cell relative to the current cell.
 * @param ds Distance covered from the photon position to the intersection
 * point.
 * @return CoordinateVector containing the coordinates of the intersection point
 * of the photon and the closest wall.
 */
CoordinateVector<> DensityGrid::get_wall_intersection(
    CoordinateVector<> &photon_origin, CoordinateVector<> &photon_direction,
    Box &cell, char &ix, char &iy, char &iz, double &ds) {
  CoordinateVector<> cell_bottom_anchor = cell.get_anchor();
  CoordinateVector<> cell_top_anchor = cell.get_top_anchor();

  // find out which cell wall the photon is going to hit next
  CoordinateVector<> next_x;
  double l;
  if (photon_direction.x() > 0.) {
    // if the photon starts at \vec{o} and travels in the direction \vec{d},
    // the general position of the photon at any later time is given by
    // \vec{o} + l*\vec{d}, with l some positive parameter
    // we know that the photon hits the next x wall when the x-component of
    // this expression equals cell_xmax, so we can solve for l:
    l = (cell_top_anchor.x() - photon_origin.x()) / photon_direction.x();
    ix = 1;
  } else if (photon_direction.x() < 0.) {
    l = (cell_bottom_anchor.x() - photon_origin.x()) / photon_direction.x();
    ix = -1;
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
    l = 1000. * _cellside_max;
    ix = 0;
  }
  // the y and z coordinates are then trivially found
  next_x = photon_origin + l * photon_direction;
  double dx = (next_x - photon_origin).norm2();

  CoordinateVector<> next_y;
  if (photon_direction.y() > 0.) {
    l = (cell_top_anchor.y() - photon_origin.y()) / photon_direction.y();
    iy = 1;
  } else if (photon_direction.y() < 0.) {
    l = (cell_bottom_anchor.y() - photon_origin.y()) / photon_direction.y();
    iy = -1;
  } else {
    l = 1000. * _cellside_max;
    iy = 0;
  }
  next_y = photon_origin + l * photon_direction;
  double dy = (next_y - photon_origin).norm2();

  CoordinateVector<> next_z;
  if (photon_direction.z() > 0.) {
    l = (cell_top_anchor.z() - photon_origin.z()) / photon_direction.z();
    iz = 1;
  } else if (photon_direction.z() < 0.) {
    l = (cell_bottom_anchor.z() - photon_origin.z()) / photon_direction.z();
    iz = -1;
  } else {
    l = 1000. * _cellside_max;
    iz = 0;
  }
  next_z = photon_origin + l * photon_direction;
  double dz = (next_z - photon_origin).norm2();

  CoordinateVector<> next_wall;
  if (dx < dy && dx < dz) {
    next_wall = next_x;
    ds = dx;
    iy = 0;
    iz = 0;
  } else if (dy < dx && dy < dz) {
    next_wall = next_y;
    ds = dy;
    ix = 0;
    iz = 0;
  } else if (dz < dx && dz < dy) {
    next_wall = next_z;
    ds = dz;
    ix = 0;
    iy = 0;
  } else {
    // special cases: at least two of the smallest values are equal
    if (dx == dy && dx < dz) {
      // it does not matter which values we pick, they will be the same
      next_wall = next_x;
      ds = dx;
      iz = 0;
    } else if (dx == dz && dx < dy) {
      next_wall = next_x;
      ds = dx;
      iy = 0;
    } else if (dy == dz && dy < dx) {
      next_wall = next_y;
      ds = dy;
      ix = 0;
    } else {
      // all values are equal, we sit on a corner of the box
      next_wall = next_x;
      ds = dx;
    }
  }

  // ds contains the squared norm, take the square root
  ds = sqrt(ds);

  return next_wall;
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
double DensityGrid::get_distance(CoordinateVector<> photon_origin,
                                 CoordinateVector<> photon_direction,
                                 double optical_depth) {
  double S = 0.;

  // find out in which cell the photon is currently hiding
  unsigned int ix, iy, iz;
  get_cell_indices(photon_origin, ix, iy, iz);

  // while the photon has not exceeded the optical depth and is still in the box
  while (optical_depth > 0.) {
    Box cell = get_cell(ix, iy, iz);

    double ds;
    char next_ix, next_iy, next_iz;
    CoordinateVector<> next_wall = get_wall_intersection(
        photon_origin, photon_direction, cell, next_ix, next_iy, next_iz, ds);

    // get the optical depth of the path from the current photon location to the
    // cell wall, update S
    S += ds;
    optical_depth -= ds;

    photon_origin = next_wall;
    ix += next_ix;
    iy += next_iy;
    iz += next_iz;

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
