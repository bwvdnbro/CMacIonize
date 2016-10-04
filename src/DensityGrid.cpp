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
#include "CrossSections.hpp"
#include "DensityFunction.hpp"
#include "DensityValues.hpp"
#include "Photon.hpp"
using namespace std;

/**
 * @brief Constructor
 *
 * @param box Box containing the grid.
 * @param n1D Number of cells in 1 dimensions.
 * @param density_function DensityFunction that defines the density field.
 * @param cross_sections Photoionization cross sections.
 */
DensityGrid::DensityGrid(Box box, CoordinateVector< unsigned char > ncell,
                         DensityFunction &density_function,
                         CrossSections &cross_sections)
    : _box(box), _ncell(ncell), _cross_sections(cross_sections) {
  _density = new DensityValues **[_ncell.x()];
  for (unsigned int i = 0; i < _ncell.x(); ++i) {
    _density[i] = new DensityValues *[_ncell.y()];
    for (unsigned int j = 0; j < _ncell.y(); ++j) {
      _density[i][j] = new DensityValues[_ncell.z()];
    }
  }

  // fill the density grid
  double cellside_x = _box.get_sides().x() / _ncell.x();
  double cellside_y = _box.get_sides().y() / _ncell.y();
  double cellside_z = _box.get_sides().z() / _ncell.z();
  _cellside = CoordinateVector<>(cellside_x, cellside_y, cellside_z);
  for (unsigned int i = 0; i < _ncell.x(); ++i) {
    for (unsigned int j = 0; j < _ncell.y(); ++j) {
      for (unsigned int k = 0; k < _ncell.z(); ++k) {
        double x = _box.get_anchor().x() + (i + 0.5) * _cellside.x();
        double y = _box.get_anchor().y() + (j + 0.5) * _cellside.y();
        double z = _box.get_anchor().z() + (k + 0.5) * _cellside.z();
        _density[i][j][k].set_total_density(
            density_function(CoordinateVector<>(x, y, z)));
        // initialize the neutral fractions to very low values
        _density[i][j][k].set_neutral_fraction_H(1.e-6);
        _density[i][j][k].set_neutral_fraction_He(1.e-6);
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
  for (unsigned int i = 0; i < _ncell.x(); ++i) {
    for (unsigned int j = 0; j < _ncell.y(); ++j) {
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

  for (unsigned int i = 0; i < _ncell.x(); ++i) {
    for (unsigned int j = 0; j < _ncell.y(); ++j) {
      for (unsigned int k = 0; k < _ncell.z(); ++k) {
        mtot += _density[i][j][k].get_total_density() * cellvolume;
      }
    }
  }

  return mtot;
}

/**
 * @brief Get the indices of the cell containing the given coordinates.
 *
 * @param position CoordinateVector containing coordinates we want to locate.
 * @return CoordinateVector<unsigned int> containing the three indices of the
 * cell.
 */
CoordinateVector< int >
DensityGrid::get_cell_indices(CoordinateVector<> position) {
  int ix = (position.x() - _box.get_anchor().x()) / _cellside.x();
  int iy = (position.y() - _box.get_anchor().y()) / _cellside.y();
  int iz = (position.z() - _box.get_anchor().z()) / _cellside.z();
  return CoordinateVector< int >(ix, iy, iz);
}

/**
 * @brief Get the geometrical box of the cell with the given indices.
 *
 * Since we care more about the bottom front left corner and the upper back
 * right corner, the Box we return does not use a side and anchor, but rather
 * two anchors.
 *
 * @param index Indices of the cell.
 * @return Box containing the bottom front left corner and the upper back right
 * corner of the cell.
 */
Box DensityGrid::get_cell(CoordinateVector< int > index) {
  double cell_xmin = _box.get_anchor().x() + _cellside.x() * index.x();
  double cell_ymin = _box.get_anchor().y() + _cellside.y() * index.y();
  double cell_zmin = _box.get_anchor().z() + _cellside.z() * index.z();
  return Box(CoordinateVector<>(cell_xmin, cell_ymin, cell_zmin), _cellside);
}

/**
 * @brief Check whether the given index points to a valid cell.
 *
 * @param index Indices of the cell.
 * @return True if the indices are valid, false otherwise.
 */
bool DensityGrid::is_inside(CoordinateVector< int > index) {
  bool inside = (index.x() >= 0 && index.x() < _ncell.x());
  inside &= (index.y() >= 0 && index.y() < _ncell.y());
  inside &= (index.z() >= 0 && index.z() < _ncell.z());
  return inside;
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
 * @param next_index Index of the neighbouring cell, relative with respect to
 * the current cell.
 * @param ds Distance covered from the photon position to the intersection
 * point.
 * @return CoordinateVector containing the coordinates of the intersection point
 * of the photon and the closest wall.
 */
CoordinateVector<> DensityGrid::get_wall_intersection(
    CoordinateVector<> &photon_origin, CoordinateVector<> &photon_direction,
    Box &cell, CoordinateVector< char > &next_index, double &ds) {
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
    next_index[0] = 1;
  } else if (photon_direction.x() < 0.) {
    l = (cell_bottom_anchor.x() - photon_origin.x()) / photon_direction.x();
    next_index[0] = -1;
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
    next_index[0] = 0;
  }
  // the y and z coordinates are then trivially found
  next_x = photon_origin + l * photon_direction;
  double dx = (next_x - photon_origin).norm2();

  CoordinateVector<> next_y;
  if (photon_direction.y() > 0.) {
    l = (cell_top_anchor.y() - photon_origin.y()) / photon_direction.y();
    next_index[1] = 1;
  } else if (photon_direction.y() < 0.) {
    l = (cell_bottom_anchor.y() - photon_origin.y()) / photon_direction.y();
    next_index[1] = -1;
  } else {
    l = 1000. * _cellside_max;
    next_index[1] = 0;
  }
  next_y = photon_origin + l * photon_direction;
  double dy = (next_y - photon_origin).norm2();

  CoordinateVector<> next_z;
  if (photon_direction.z() > 0.) {
    l = (cell_top_anchor.z() - photon_origin.z()) / photon_direction.z();
    next_index[2] = 1;
  } else if (photon_direction.z() < 0.) {
    l = (cell_bottom_anchor.z() - photon_origin.z()) / photon_direction.z();
    next_index[2] = -1;
  } else {
    l = 1000. * _cellside_max;
    next_index[2] = 0;
  }
  next_z = photon_origin + l * photon_direction;
  double dz = (next_z - photon_origin).norm2();

  CoordinateVector<> next_wall;
  if (dx < dy && dx < dz) {
    next_wall = next_x;
    ds = dx;
    next_index[1] = 0;
    next_index[2] = 0;
  } else if (dy < dx && dy < dz) {
    next_wall = next_y;
    ds = dy;
    next_index[0] = 0;
    next_index[2] = 0;
  } else if (dz < dx && dz < dy) {
    next_wall = next_z;
    ds = dz;
    next_index[0] = 0;
    next_index[1] = 0;
  } else {
    // special cases: at least two of the smallest values are equal
    if (dx == dy && dx < dz) {
      // it does not matter which values we pick, they will be the same
      next_wall = next_x;
      ds = dx;
      next_index[2] = 0;
    } else if (dx == dz && dx < dy) {
      next_wall = next_x;
      ds = dx;
      next_index[1] = 0;
    } else if (dy == dz && dy < dx) {
      next_wall = next_y;
      ds = dy;
      next_index[0] = 0;
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
 * @brief Let the given Photon travel through the density grid until the given
 * optical depth is reached.
 *
 * @param photon Photon.
 * @param optical_depth Optical depth the photon should travel in total.
 * @return True if the Photon is still in the box after the optical depth has
 * been reached, false otherwise.
 */
bool DensityGrid::interact(Photon &photon, double optical_depth) {
  double S = 0.;

  CoordinateVector<> photon_origin = photon.get_position();
  CoordinateVector<> photon_direction = photon.get_direction();

  // find out in which cell the photon is currently hiding
  CoordinateVector< int > index = get_cell_indices(photon_origin);

  double xsecH =
      _cross_sections.get_cross_section(ELEMENT_H, photon.get_energy());
  double xsecHe =
      _cross_sections.get_cross_section(ELEMENT_He, photon.get_energy());
  // Helium mass fraction. Should be a parameter.
  double AHe = 0.1;

  // while the photon has not exceeded the optical depth and is still in the box
  while (is_inside(index) && optical_depth > 0.) {
    Box cell = get_cell(index);

    double ds;
    CoordinateVector< char > next_index;
    CoordinateVector<> next_wall = get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);

    // get the optical depth of the path from the current photon location to the
    // cell wall, update S
    DensityValues &density = _density[index.x()][index.y()][index.z()];
    double tau = ds * density.get_total_density() *
                 (xsecH * density.get_neutral_fraction_H() +
                  AHe * xsecHe * density.get_neutral_fraction_He());
    optical_depth -= tau;

    index += next_index;

    // if the optical depth exceeds or equals the wanted value: exit the loop

    // if the optical depth exceeded the wanted value: find out where in the
    // cell
    // we end up, and correct S
    if (optical_depth < 0.) {
      double Scorr = ds * optical_depth / tau;
      // order is important here!
      photon_origin += (next_wall - photon_origin) * (ds + Scorr) / ds;
      ds += Scorr;
    } else {
      photon_origin = next_wall;
    }

    // ds is now the actual distance travelled in the cell
    // update contributions to mean intensity integrals
    density.increase_mean_intensity_H(ds * xsecH);
    density.increase_mean_intensity_He(ds * xsecHe);

    S += ds;
  }

  photon.set_position(photon_origin);

  return is_inside(index);
}
