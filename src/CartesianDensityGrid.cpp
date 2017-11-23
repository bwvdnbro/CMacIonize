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
 * @file CartesianDensityGrid.cpp
 *
 * @brief Cartesian density grid: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CartesianDensityGrid.hpp"
#include "DensityFunction.hpp"
#include "DensityValues.hpp"
#include "IonizationStateCalculator.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Photon.hpp"
#include "RecombinationRates.hpp"
#include "SimulationBox.hpp"
#include "Timer.hpp"
#include <cfloat>
#include <sstream>

/**
 * @brief Constructor
 *
 * @param simulation_box Simulation box (in m).
 * @param ncell Number of cells for each dimension.
 * @param periodic Periodicity flags.
 * @param hydro Hydro flag.
 * @param log Log to write log messages to.
 */
CartesianDensityGrid::CartesianDensityGrid(
    const Box<> &simulation_box, CoordinateVector< int_fast32_t > ncell,
    CoordinateVector< bool > periodic, bool hydro, Log *log)
    : DensityGrid(simulation_box, periodic, hydro, log), _box(simulation_box),
      _periodicity_flags(periodic), _ncell(ncell), _log(log) {

  if (_log) {
    _log->write_status(
        "Creating CartesianDensityGrid of ", _ncell.x(), " x ", _ncell.y(),
        " x ", _ncell.z(), " inside a box with anchor [", _box.get_anchor().x(),
        " m, ", _box.get_anchor().y(), " m, ", _box.get_anchor().z(),
        " m] and sides [", _box.get_sides().x(), " m, ", _box.get_sides().y(),
        " m, ", _box.get_sides().z(), " m]...");
    if (_periodicity_flags.x()) {
      _log->write_status("x boundary is periodic.");
    } else {
      _log->write_status("x boundary is not periodic.");
    }
    if (_periodicity_flags.y()) {
      _log->write_status("y boundary is periodic.");
    } else {
      _log->write_status("y boundary is not periodic.");
    }
    if (_periodicity_flags.z()) {
      _log->write_status("z boundary is periodic.");
    } else {
      _log->write_status("z boundary is not periodic.");
    }
  }

  const cellsize_t totnumcell = _ncell.x() * _ncell.y() * _ncell.z();
  allocate_memory(totnumcell);

  double cellside_x = _box.get_sides().x() / _ncell.x();
  double cellside_y = _box.get_sides().y() / _ncell.y();
  double cellside_z = _box.get_sides().z() / _ncell.z();
  _cellside = CoordinateVector<>(cellside_x, cellside_y, cellside_z);
  _inverse_cellside[0] = 1. / _cellside.x();
  _inverse_cellside[1] = 1. / _cellside.y();
  _inverse_cellside[2] = 1. / _cellside.z();

  if (_log) {
    _log->write_info("Cell size is ", _cellside.x(), " m x ", _cellside.y(),
                     " m x ", _cellside.z(), " m.");
    _log->write_status("Done creating grid.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - number of cells: Number of cells in the grid (default: [64, 64, 64])
 *
 * @param simulation_box SimulationBox.
 * @param parameters ParameterFile to read.
 * @param hydro Is hydrodynamics enabled?
 * @param log Log to write log messages to.
 */
CartesianDensityGrid::CartesianDensityGrid(const SimulationBox &simulation_box,
                                           ParameterFile &parameters,
                                           bool hydro, Log *log)
    : CartesianDensityGrid(
          simulation_box.get_box(),
          parameters.get_value< CoordinateVector< int_fast32_t > >(
              "DensityGrid:number of cells",
              CoordinateVector< int_fast32_t >(64)),
          simulation_box.get_periodicity(), hydro, log) {}

/**
 * @brief Initialize the cells in the grid.
 *
 * @param block Block that should be initialized by this MPI process.
 * @param density_function DensityFunction to use.
 */
void CartesianDensityGrid::initialize(
    std::pair< cellsize_t, cellsize_t > &block,
    DensityFunction &density_function) {
  DensityGrid::initialize(block, density_function);
  DensityGrid::set_densities(block, density_function);
}

/**
 * @brief Get the total number of cells in this grid.
 *
 * @return Total number of cells.
 */
cellsize_t CartesianDensityGrid::get_number_of_cells() const {
  return _ncell.x() * _ncell.y() * _ncell.z();
}

/**
 * @brief Get the indices of the cell containing the given coordinates.
 *
 * @param position CoordinateVector containing coordinates we want to locate.
 * @return CoordinateVector<unsigned int> containing the three indices of the
 * cell.
 */
CoordinateVector< int_fast32_t >
CartesianDensityGrid::get_cell_indices(CoordinateVector<> position) const {
  int_fast32_t ix =
      (position.x() - _box.get_anchor().x()) * _inverse_cellside.x();
  int_fast32_t iy =
      (position.y() - _box.get_anchor().y()) * _inverse_cellside.y();
  int_fast32_t iz =
      (position.z() - _box.get_anchor().z()) * _inverse_cellside.z();
  return CoordinateVector< int_fast32_t >(ix, iy, iz);
}

/**
 * @brief Get the geometrical box of the cell with the given indices.
 *
 * @param index Indices of the cell.
 * @return Box containing the bottom front left corner and the upper back right
 * corner of the cell (in m).
 */
Box<> CartesianDensityGrid::get_cell(
    CoordinateVector< int_fast32_t > index) const {
  double cell_xmin = _box.get_anchor().x() + _cellside.x() * index.x();
  double cell_ymin = _box.get_anchor().y() + _cellside.y() * index.y();
  double cell_zmin = _box.get_anchor().z() + _cellside.z() * index.z();
  return Box<>(CoordinateVector<>(cell_xmin, cell_ymin, cell_zmin), _cellside);
}

/**
 * @brief Check whether the given index points to a valid cell.
 *
 * This method also applies periodic boundary conditions (if applicable).
 *
 * @param index Indices of the cell.
 * @param position Current position of the photon.
 * @return True if the indices are valid, false otherwise.
 */
bool CartesianDensityGrid::is_inside(CoordinateVector< int_fast32_t > &index,
                                     CoordinateVector<> &position) const {
  bool inside = true;
  if (!_periodicity_flags.x()) {
    inside &= (index.x() >= 0 && index.x() < _ncell.x());
  } else {
    if (index.x() < 0) {
      index[0] = _ncell.x() - 1;
      position[0] += _box.get_sides().x();
    }
    if (index.x() >= _ncell.x()) {
      index[0] = 0;
      position[0] -= _box.get_sides().x();
    }
  }
  if (!_periodicity_flags.y()) {
    inside &= (index.y() >= 0 && index.y() < _ncell.y());
  } else {
    if (index.y() < 0) {
      index[1] = _ncell.y() - 1;
      position[1] += _box.get_sides().y();
    }
    if (index.y() >= _ncell.y()) {
      index[1] = 0;
      position[1] -= _box.get_sides().y();
    }
  }
  if (!_periodicity_flags.z()) {
    inside &= (index.z() >= 0 && index.z() < _ncell.z());
  } else {
    if (index.z() < 0) {
      index[2] = _ncell.z() - 1;
      position[2] += _box.get_sides().z();
    }
    if (index.z() >= _ncell.z()) {
      index[2] = 0;
      position[2] -= _box.get_sides().z();
    }
  }
  return inside;
}

/**
 * @brief Check whether the given index points to a valid cell.
 *
 * This version ignores the boundary flags and always assumes open boundaries.
 *
 * @param index Indices of the cell.
 * @param position Current position of the photon.
 * @return True if the indices are valid, false otherwise.
 */
bool CartesianDensityGrid::is_inside_non_periodic(
    CoordinateVector< int_fast32_t > &index,
    CoordinateVector<> &position) const {
  bool inside = true;
  inside &= (index.x() >= 0 && index.x() < _ncell.x());
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
 * @param photon_origin Current position of the photon (in m).
 * @param photon_direction Direction the photon is travelling in.
 * @param inverse_photon_direction Inverse of the photon travel direction.
 * @param cell Cell in which the photon currently resides.
 * @param next_index Index of the neighbouring cell, relative with respect to
 * the current cell.
 * @param ds Distance covered from the photon position to the intersection
 * point (in m).
 * @return CoordinateVector containing the coordinates of the intersection point
 * of the photon and the closest wall (in m).
 */
CoordinateVector<> CartesianDensityGrid::get_wall_intersection(
    const CoordinateVector<> &photon_origin,
    const CoordinateVector<> &photon_direction,
    const CoordinateVector<> &inverse_photon_direction, const Box<> &cell,
    CoordinateVector< int_fast8_t > &next_index, double &ds) {

  const CoordinateVector<> cell_bottom_anchor = cell.get_anchor();
  const CoordinateVector<> cell_top_anchor = cell.get_top_anchor();

  const double dx = (photon_direction.x() > 0.)
                        ? (cell_top_anchor.x() - photon_origin.x()) *
                              inverse_photon_direction.x()
                        : ((photon_direction.x() < 0.)
                               ? (cell_bottom_anchor.x() - photon_origin.x()) *
                                     inverse_photon_direction.x()
                               : DBL_MAX);
  const double dy = (photon_direction.y() > 0.)
                        ? (cell_top_anchor.y() - photon_origin.y()) *
                              inverse_photon_direction.y()
                        : ((photon_direction.y() < 0.)
                               ? (cell_bottom_anchor.y() - photon_origin.y()) *
                                     inverse_photon_direction.y()
                               : DBL_MAX);
  const double dz = (photon_direction.z() > 0.)
                        ? (cell_top_anchor.z() - photon_origin.z()) *
                              inverse_photon_direction.z()
                        : ((photon_direction.z() < 0.)
                               ? (cell_bottom_anchor.z() - photon_origin.z()) *
                                     inverse_photon_direction.z()
                               : DBL_MAX);

  ds = std::min(dx, std::min(dy, dz));

  next_index[0] = (dx == ds) ? ((photon_direction.x() > 0.) ? 1 : -1) : 0;
  next_index[1] = (dy == ds) ? ((photon_direction.y() > 0.) ? 1 : -1) : 0;
  next_index[2] = (dz == ds) ? ((photon_direction.z() > 0.) ? 1 : -1) : 0;

  return photon_origin + ds * photon_direction;
}

/**
 * @brief Get the total optical depth traversed by the given Photon until it
 * reaches the boundaries of the simulation box.
 *
 * @param photon Photon.
 * @return Total optical depth along the photon's path before it reaches the
 * boundaries of the simulation box.
 */
double CartesianDensityGrid::integrate_optical_depth(const Photon &photon) {

  double optical_depth = 0.;

  CoordinateVector<> photon_origin = photon.get_position();
  const CoordinateVector<> photon_direction = photon.get_direction();
  const CoordinateVector<> inverse_photon_direction =
      photon.get_inverse_direction();

  // find out in which cell the photon is currently hiding
  CoordinateVector< int_fast32_t > index = get_cell_indices(photon_origin);

  // while the photon is still in the box
  while (is_inside(index, photon_origin)) {
    Box<> cell = get_cell(index);

    double ds;
    CoordinateVector< int_fast8_t > next_index;
    CoordinateVector<> next_wall =
        get_wall_intersection(photon_origin, photon_direction,
                              inverse_photon_direction, cell, next_index, ds);

    // get the optical depth of the path from the current photon location to the
    // cell wall, update S
    DensityGrid::iterator it(get_long_index(index), *this);

    // Helium abundance. Should be a parameter.
    optical_depth +=
        get_optical_depth(ds, it.get_ionization_variables(), photon);

    photon_origin = next_wall;
    index += next_index;
  }

  return optical_depth;
}

/**
 * @brief Let the given Photon travel through the density grid until the given
 * optical depth is reached.
 *
 * @param photon Photon.
 * @param optical_depth Optical depth the photon should travel in total
 * (dimensionless).
 * @return DensityGrid::iterator pointing to the cell the photon was last in,
 * or DensityGrid::end() if the photon left the box.
 */
DensityGrid::iterator CartesianDensityGrid::interact(Photon &photon,
                                                     double optical_depth) {

  double S = 0.;

  CoordinateVector<> photon_origin = photon.get_position();
  const CoordinateVector<> photon_direction = photon.get_direction();
  const CoordinateVector<> inverse_photon_direction =
      photon.get_inverse_direction();

  // find out in which cell the photon is currently hiding
  CoordinateVector< int_fast32_t > index = get_cell_indices(photon_origin);

  uint_fast32_t ncell = 0;
  DensityGrid::iterator last_cell = end();
  // while the photon has not exceeded the optical depth and is still in the box
  while (is_inside(index, photon_origin) && optical_depth > 0.) {
    ++ncell;
    Box<> cell = get_cell(index);

    double ds;
    CoordinateVector< int_fast8_t > next_index;
    CoordinateVector<> next_wall =
        get_wall_intersection(photon_origin, photon_direction,
                              inverse_photon_direction, cell, next_index, ds);

    // get the optical depth of the path from the current photon location to the
    // cell wall, update S
    DensityGrid::iterator it(get_long_index(index), *this);
    last_cell = it;

    // Helium abundance. Should be a parameter.
    double tau = get_optical_depth(ds, it.get_ionization_variables(), photon);
    optical_depth -= tau;

    // if the optical depth exceeds or equals the wanted value: exit the loop

    // if the optical depth exceeded the wanted value: find out where in the
    // cell we end up, and correct S
    if (optical_depth < 0.) {
      double Scorr = ds * optical_depth / tau;
      // order is important here!
      photon_origin += (next_wall - photon_origin) * (ds + Scorr) / ds;
      ds += Scorr;
    } else {
      photon_origin = next_wall;
      // we don't want to do this if the photon does not actually leave the
      // cell, since this might trigger a periodic boundary position change:
      // the position of the photon is adapted for a traversal of the periodic
      // boundaries, while it does not traverse them
      // if there are no periodic boundaries, it does not matter where we
      // update the index, since the index will not be used any more after
      // leaving this loop
      index += next_index;
    }

    // ds is now the actual distance travelled in the cell
    // update contributions to mean intensity integrals
    update_integrals(ds, it, photon);

    S += ds;
  }

  if (ncell == 0 && optical_depth > 0.) {
    cmac_error("Photon leaves the system immediately (position: %g %g %g, "
               "direction: %g %g %g)!",
               photon_origin.x(), photon_origin.y(), photon_origin.z(),
               photon_direction.x(), photon_direction.y(),
               photon_direction.z());
  }

  photon.set_position(photon_origin);

  if (!is_inside(index, photon_origin)) {
    last_cell = end();
  }

  return last_cell;
}

/**
 * @brief Get the total line emission along a ray with the given origin and
 * direction.
 *
 * @param origin Origin of the ray (in m).
 * @param direction Direction of the ray.
 * @param line EmissionLine name of the line to trace.
 * @return Accumulated emission along the ray (in J m^-2 s^-1).
 */
double CartesianDensityGrid::get_total_emission(CoordinateVector<> origin,
                                                CoordinateVector<> direction,
                                                EmissionLine line) {

  const CoordinateVector<> inverse_direction(
      1. / direction.x(), 1. / direction.y(), 1. / direction.z());

  double S = 0.;

  // find out in which cell the origin lies
  CoordinateVector< int_fast32_t > index = get_cell_indices(origin);

  uint_fast32_t ncell = 0;
  // while the photon has not exceeded the optical depth and is still in the box
  while (is_inside_non_periodic(index, origin)) {
    ++ncell;
    Box<> cell = get_cell(index);

    double ds;
    CoordinateVector< int_fast8_t > next_index;
    CoordinateVector<> next_wall = get_wall_intersection(
        origin, direction, inverse_direction, cell, next_index, ds);

    // get the optical depth of the path from the current photon location to the
    // cell wall, update S
    DensityGrid::iterator it(get_long_index(index), *this);

    S += ds * it.get_emissivities()->get_emissivity(line);

    origin = next_wall;
    index += next_index;
  }

  return S;
}

/**
 * @brief Get the neighbours of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return std::vector containing the neighbours of the cell.
 */
std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                         CoordinateVector<>, double, CoordinateVector<> > >
CartesianDensityGrid::get_neighbours(cellsize_t index) {

  std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                           CoordinateVector<>, double, CoordinateVector<> > >
      ngbs;

  double sidelength[3];
  sidelength[0] = _box.get_sides().x() / _ncell.x();
  sidelength[1] = _box.get_sides().y() / _ncell.y();
  sidelength[2] = _box.get_sides().z() / _ncell.z();
  double surface_area[3];
  surface_area[0] = sidelength[1] * sidelength[2];
  surface_area[1] = sidelength[0] * sidelength[2];
  surface_area[2] = sidelength[0] * sidelength[1];
  const CoordinateVector< int_fast32_t > cellindices = get_indices(index);
  const CoordinateVector<> cell_midpoint = get_cell_midpoint(cellindices);
  for (uint_fast32_t i = 0; i < 3; ++i) {
    if (cellindices[i] > 0) {
      CoordinateVector< int_fast32_t > ngb_low(cellindices);
      ngb_low[i] -= 1;
      CoordinateVector<> correction;
      correction[i] -= 0.5 * sidelength[i];
      CoordinateVector<> midpoint = cell_midpoint + correction;
      CoordinateVector<> normal;
      normal[i] = -1.;
      DensityGrid::iterator ngb_it =
          DensityGrid::iterator(get_long_index(ngb_low), *this);
      ngbs.push_back(
          std::make_tuple(ngb_it, midpoint, normal, surface_area[i],
                          ngb_it.get_cell_midpoint() - cell_midpoint));
    } else {
      if (_periodicity_flags[i]) {
        CoordinateVector< int_fast32_t > ngb_low(cellindices);
        ngb_low[i] = _ncell[i] - 1;
        CoordinateVector<> correction;
        correction[i] -= 0.5 * sidelength[i];
        CoordinateVector<> midpoint = cell_midpoint + correction;
        CoordinateVector<> normal;
        normal[i] = -1.;
        DensityGrid::iterator ngb_it(get_long_index(ngb_low), *this);
        CoordinateVector<> rel_pos = ngb_it.get_cell_midpoint() - cell_midpoint;
        rel_pos[i] -= _box.get_sides()[i];
        ngbs.push_back(std::make_tuple(ngb_it, midpoint, normal,
                                       surface_area[i], rel_pos));
      } else {
        // mirror the cell: reflective boundaries
        CoordinateVector<> correction;
        correction[i] -= 0.5 * sidelength[i];
        CoordinateVector<> midpoint = cell_midpoint + correction;
        CoordinateVector<> normal;
        normal[i] = -1.;
        ngbs.push_back(std::make_tuple(end(), midpoint, normal, surface_area[i],
                                       2. * correction));
      }
    }

    if (cellindices[i] < _ncell[i] - 1) {
      CoordinateVector< int_fast32_t > ngb_high(cellindices);
      ngb_high[i] += 1;
      CoordinateVector<> correction;
      correction[i] += 0.5 * sidelength[i];
      CoordinateVector<> midpoint = cell_midpoint + correction;
      CoordinateVector<> normal;
      normal[i] = 1.;
      DensityGrid::iterator ngb_it(get_long_index(ngb_high), *this);
      ngbs.push_back(
          std::make_tuple(ngb_it, midpoint, normal, surface_area[i],
                          ngb_it.get_cell_midpoint() - cell_midpoint));
    } else {
      if (_periodicity_flags[i]) {
        CoordinateVector< int_fast32_t > ngb_high(cellindices);
        ngb_high[i] = 0;
        CoordinateVector<> correction;
        correction[i] += 0.5 * sidelength[i];
        CoordinateVector<> midpoint = cell_midpoint + correction;
        CoordinateVector<> normal;
        normal[i] = 1.;
        DensityGrid::iterator ngb_it(get_long_index(ngb_high), *this);
        CoordinateVector<> rel_pos = ngb_it.get_cell_midpoint() - cell_midpoint;
        rel_pos[i] += _box.get_sides()[i];
        ngbs.push_back(std::make_tuple(ngb_it, midpoint, normal,
                                       surface_area[i], rel_pos));
      } else {
        // mirror the cell: reflective boundaries
        CoordinateVector<> correction;
        correction[i] += 0.5 * sidelength[i];
        CoordinateVector<> midpoint = cell_midpoint + correction;
        CoordinateVector<> normal;
        normal[i] = 1.;
        ngbs.push_back(std::make_tuple(end(), midpoint, normal, surface_area[i],
                                       2. * correction));
      }
    }
  }

  return ngbs;
}

/**
 * @brief Get the faces of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Faces of the cell.
 */
std::vector< Face > CartesianDensityGrid::get_faces(cellsize_t index) const {

  const double sidelength[3] = {_box.get_sides().x() / _ncell.x(),
                                _box.get_sides().y() / _ncell.y(),
                                _box.get_sides().z() / _ncell.z()};
  const CoordinateVector< int_fast32_t > cellindices = get_indices(index);
  const CoordinateVector<> cell_midpoint = get_cell_midpoint(cellindices);

  std::vector< Face > faces;

  std::vector< CoordinateVector<> > vertices(4);

  // negative x face
  vertices[0] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[1] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[2] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  vertices[3] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  faces.push_back(
      Face(CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                              cell_midpoint.y(), cell_midpoint.z()),
           vertices));
  // positive x face
  vertices[0] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[1] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[2] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  vertices[3] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  faces.push_back(
      Face(CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                              cell_midpoint.y(), cell_midpoint.z()),
           vertices));

  // negative y face
  vertices[0] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[1] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[2] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  vertices[3] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  faces.push_back(
      Face(CoordinateVector<>(cell_midpoint.x(),
                              cell_midpoint.y() - 0.5 * sidelength[1],
                              cell_midpoint.z()),
           vertices));
  // positive y face
  vertices[0] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[1] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[2] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  vertices[3] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  faces.push_back(
      Face(CoordinateVector<>(cell_midpoint.x(),
                              cell_midpoint.y() + 0.5 * sidelength[1],
                              cell_midpoint.z()),
           vertices));

  // negative z face
  vertices[0] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[1] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[2] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  vertices[3] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() - 0.5 * sidelength[2]);
  faces.push_back(
      Face(CoordinateVector<>(cell_midpoint.x(), cell_midpoint.y(),
                              cell_midpoint.z() - 0.5 * sidelength[2]),
           vertices));
  // positive z face
  vertices[0] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  vertices[1] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() - 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  vertices[2] = CoordinateVector<>(cell_midpoint.x() + 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  vertices[3] = CoordinateVector<>(cell_midpoint.x() - 0.5 * sidelength[0],
                                   cell_midpoint.y() + 0.5 * sidelength[1],
                                   cell_midpoint.z() + 0.5 * sidelength[2]);
  faces.push_back(
      Face(CoordinateVector<>(cell_midpoint.x(), cell_midpoint.y(),
                              cell_midpoint.z() + 0.5 * sidelength[2]),
           vertices));

  return faces;
}
