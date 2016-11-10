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
 * @file AMRDensityGrid.hpp
 *
 * @brief AMR density grid: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AMRDENSITYGRID_HPP
#define AMRDENSITYGRID_HPP

#include "AMRGrid.hpp"
#include "AMRRefinementSchemeFactory.hpp"
#include "DensityGrid.hpp"
#include "ParameterFile.hpp"

#include <cfloat>
#include <ostream>

/**
 * @brief AMR density grid.
 */
class AMRDensityGrid : public DensityGrid {
private:
  /*! @brief AMRGrid used as grid. */
  AMRGrid< DensityValues > _grid;

  /*! @brief DensityFunction used to initialize densities. */
  DensityFunction &_density_function;

  /*! @brief AMRRefinementScheme used to refine cells. */
  AMRRefinementScheme *_refinement_scheme;

  /**
   * @brief Get the largest odd factor of the given number.
   *
   * This is the number you get by iteratively dividing the number by two, until
   * the result is no longer even.
   *
   * @param number Number to decompose.
   * @return Largest odd factor of the number.
   */
  inline int get_largest_odd_factor(int number) {
    while ((number % 2) == 0) {
      number >>= 1;
    }
    return number;
  }

  /**
   * @brief Get the largest power of two factor of the given number.
   *
   * This is the factor you get by multiplying all factors of two contained
   * in the number.
   *
   * @param number Number to decompose.
   * @return Largest power of two factor.
   */
  inline int get_power_of_two(int number) {
    return number / get_largest_odd_factor(number);
  }

  /**
   * @brief Refine the given cell using the given refinement scheme.
   *
   * @param refinement_scheme AMRRefinementScheme to apply.
   * @param index Index of the cell.
   * @param level Refinement level of the cell.
   * @param midpoint CoordinateVector<> specifying the midpoint of the cell.
   * @param values DensityValues of the cell.
   * @param density_function DensityFunction used to initialize newly created
   * cells.
   * @param next_index Index of the next cell at the current level.
   */
  inline void refine_cell(AMRRefinementScheme *refinement_scheme,
                          unsigned long index, unsigned char level,
                          CoordinateVector<> midpoint, DensityValues &values,
                          DensityFunction &density_function,
                          unsigned long next_index) {
    if (refinement_scheme->refine(level, midpoint, values)) {
      // initialize and check refinement criterion for children
      unsigned long child_index = _grid.refine_cell(index);
      while (child_index != next_index) {
        DensityValues &cell = _grid[child_index];
        CoordinateVector<> child_midpoint = _grid.get_midpoint(child_index);
        cell.set_total_density(density_function(child_midpoint));
        // copy all other values from parent
        cell.set_neutral_fraction_H(values.get_neutral_fraction_H());
        cell.set_neutral_fraction_He(values.get_neutral_fraction_He());
        cell.set_temperature(values.get_temperature());
        cell.set_helium_abundance(values.get_helium_abundance());
        cell.set_old_neutral_fraction_H(values.get_old_neutral_fraction_H());
        // set reemission probabilities
        set_reemission_probabilities(cell.get_temperature(), cell);
        unsigned long next_child_index = _grid.get_next_key(child_index);
        refine_cell(refinement_scheme, child_index, level + 1, child_midpoint,
                    cell, density_function, next_child_index);
        child_index = next_child_index;
      }
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the grid.
   * @param ncell Number of cells in the low resolution grid.
   * @param helium_abundance Helium abundance (relative w.r.t. hydrogen).
   * @param initial_temperature Initial temperature of the gas (in K).
   * @param density_function DensityFunction that defines the density field.
   * @param refinement_scheme Refinement scheme used to refine cells. Memory
   * management for this pointer is taken over by this class.
   * @param periodic Periodicity flags.
   * @param log Log to write logging info to.
   */
  AMRDensityGrid(
      Box box, CoordinateVector< int > ncell, double helium_abundance,
      double initial_temperature, DensityFunction &density_function,
      AMRRefinementScheme *refinement_scheme = nullptr,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      Log *log = nullptr)
      : DensityGrid(box, periodic, log), _density_function(density_function),
        _refinement_scheme(refinement_scheme) {
    // find the smallest number of blocks that fits the requested top level grid
    // for one dimension, this is the largest odd factor in that dimension
    // for all three dimensions, this is the factor you get when you divide the
    // requested number of cells in that dimensions by the smallest common power
    // of two of all three dimensions.
    int power_of_2_x = get_power_of_two(ncell.x());
    int power_of_2_y = get_power_of_two(ncell.y());
    int power_of_2_z = get_power_of_two(ncell.z());
    int power_of_2 = std::min(power_of_2_x, power_of_2_y);
    power_of_2 = std::min(power_of_2, power_of_2_z);
    CoordinateVector< int > nblock = ncell / power_of_2;
    _grid = AMRGrid< DensityValues >(box, nblock);

    // find out how many cells each block should have at the lowest level
    // this is just the power in power_of_2
    unsigned char level = 0;
    while (power_of_2 > 1) {
      power_of_2 >>= 1;
      ++level;
    }
    _grid.create_all_cells(level);

    if (_log) {
      int levelint = level;
      _log->write_status("Created AMRGrid with ", nblock.x(), "x", nblock.y(),
                         "x", nblock.z(), " top level blocks, going ", levelint,
                         " levels deep.");
    }

    initialize(initial_temperature, helium_abundance, density_function);

    // apply mesh refinement
    if (_refinement_scheme) {
      if (_log) {
        _log->write_status("Applying refinement.");
      }

      // refining a cell could affect the iterator
      // we therefore take care to update the iterator before refining the cell
      // this way, we already have the next cell in the iteration before we
      // start playing around with the current one.
      auto it = begin();
      while (it != end()) {
        unsigned long index = it.get_index();
        unsigned char level = _grid.get_level(index);
        DensityValues &values = it.get_values();
        CoordinateVector<> midpoint = it.get_cell_midpoint();
        ++it;
        refine_cell(_refinement_scheme, index, level, midpoint, values,
                    _density_function, it.get_index());
      }

      if (_log) {
        _log->write_status("Number of cells after refinement: ",
                           _grid.get_number_of_cells());
      }
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read.
   * @param density_function DensityFunction used to set the densities in each
   * cell.
   * @param log Log to write log messages to.
   */
  AMRDensityGrid(ParameterFile &params, DensityFunction &density_function,
                 Log *log)
      : AMRDensityGrid(
            Box(params.get_physical_vector< QUANTITY_LENGTH >(
                    "densityfunction.box_anchor", "[0. m, 0. m, 0. m]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densityfunction.box_sides", "[1. m, 1. m, 1. m]")),
            params.get_value< CoordinateVector< int > >(
                "densitygrid.ncell", CoordinateVector< int >(64)),
            params.get_value< double >("helium_abundance", 0.1),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "densitygrid.initial_temperature", "8000. K"),
            density_function, AMRRefinementSchemeFactory::generate(params, log),
            params.get_value< CoordinateVector< bool > >(
                "densitygrid.periodicity", CoordinateVector< bool >(false)),
            log) {}

  virtual ~AMRDensityGrid() {
    if (_refinement_scheme != nullptr) {
      delete _refinement_scheme;
    }
  }

  /**
   * @brief Reset the mean intensity counters, update the reemission
   * probabilities, and reapply the refinement scheme to all cells.
   */
  virtual void reset_grid() {
    auto it = begin();
    while (it != end()) {
      unsigned long index = it.get_index();
      unsigned char level = _grid.get_level(index);
      DensityValues &values = it.get_values();
      values.reset_mean_intensities();
      CoordinateVector<> midpoint = it.get_cell_midpoint();
      ++it;
      refine_cell(_refinement_scheme, index, level, midpoint, values,
                  _density_function, it.get_index());
    }
  }

  /**
   * @brief Get the number of (lowest level) cells in the grid.
   *
   * The lowest level cells are all the cells that have no children and that,
   * together, cover the whole box exactly once.
   *
   * @return Number of lowest level AMR cells.
   */
  virtual unsigned int get_number_of_cells() {
    return _grid.get_number_of_cells();
  }

  /**
   * @brief Get the index of the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return Index of the cell containing that position.
   */
  virtual unsigned long get_cell_index(CoordinateVector<> position) {
    return _grid.get_key(position);
  }

  /**
   * @brief Get the midpoint of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Midpoint of that cell (in m).
   */
  virtual CoordinateVector<> get_cell_midpoint(unsigned long index) {
    return _grid.get_midpoint(index);
  }

  /**
   * @brief Get the values stored in the cell with the given index.
   *
   * @param index Index of a cell.
   * @return DensityValues stored in that cell.
   */
  virtual DensityValues &get_cell_values(unsigned long index) {
    return _grid[index];
  }

  /**
   * @brief Get the volume of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Volume of that cell (in m^3).
   */
  virtual double get_cell_volume(unsigned long index) {
    return _grid.get_volume(index);
  }

  /**
   * @brief Get the intersection point of a photon with one of the walls of a
   * cell.
   *
   * We assume the photon is in the cell and find the closest intersection point
   * with one of the walls of the cell (in the travel direction). We also set
   * the index of the cell at the other side of the wall.
   *
   * @param photon_origin Current position of the photon (in m).
   * @param photon_direction Direction the photon is travelling in.
   * @param cell Geometry of the cell in which the photon currently resides.
   * @param index Index of the cell that contains the photon. Is updated to the
   * new neighbouring cell index.
   * @param ds Distance covered from the photon position to the intersection
   * point (in m).
   * @param periodic_correction CoordinateVector used to store periodic
   * correction terms that will be applied to the photon position if we jump to
   * the neighbouring cell.
   * @return CoordinateVector containing the coordinates of the intersection
   * point of the photon and the closest wall (in m).
   */
  inline CoordinateVector<>
  get_wall_intersection(CoordinateVector<> &photon_origin,
                        CoordinateVector<> &photon_direction, Box &cell,
                        unsigned long &index, double &ds,
                        CoordinateVector<> &periodic_correction) {
    CoordinateVector<> cell_bottom_anchor = cell.get_anchor();
    CoordinateVector<> cell_top_anchor = cell.get_top_anchor();

    CoordinateVector< char > next_direction;

    // find out which cell wall the photon is going to hit next
    CoordinateVector<> next_x;
    double l;
    if (photon_direction.x() > 0.) {
      l = (cell_top_anchor.x() - photon_origin.x()) / photon_direction.x();
      next_direction[0] = 1;
    } else if (photon_direction.x() < 0.) {
      l = (cell_bottom_anchor.x() - photon_origin.x()) / photon_direction.x();
      next_direction[0] = -1;
    } else {
      l = DBL_MAX;
      next_direction[0] = 0;
    }
    // the y and z coordinates are then trivially found
    next_x = photon_origin + l * photon_direction;
    double dx = (next_x - photon_origin).norm2();

    CoordinateVector<> next_y;
    if (photon_direction.y() > 0.) {
      l = (cell_top_anchor.y() - photon_origin.y()) / photon_direction.y();
      next_direction[1] = 1;
    } else if (photon_direction.y() < 0.) {
      l = (cell_bottom_anchor.y() - photon_origin.y()) / photon_direction.y();
      next_direction[1] = -1;
    } else {
      l = DBL_MAX;
      next_direction[1] = 0;
    }
    next_y = photon_origin + l * photon_direction;
    double dy = (next_y - photon_origin).norm2();

    CoordinateVector<> next_z;
    if (photon_direction.z() > 0.) {
      l = (cell_top_anchor.z() - photon_origin.z()) / photon_direction.z();
      next_direction[2] = 1;
    } else if (photon_direction.z() < 0.) {
      l = (cell_bottom_anchor.z() - photon_origin.z()) / photon_direction.z();
      next_direction[2] = -1;
    } else {
      l = DBL_MAX;
      next_direction[2] = 0;
    }
    next_z = photon_origin + l * photon_direction;
    double dz = (next_z - photon_origin).norm2();

    CoordinateVector<> next_wall;
    if (dx < dy && dx < dz) {
      next_wall = next_x;
      ds = dx;
      next_direction[1] = 0;
      next_direction[2] = 0;
    } else if (dy < dx && dy < dz) {
      next_wall = next_y;
      ds = dy;
      next_direction[0] = 0;
      next_direction[2] = 0;
    } else if (dz < dx && dz < dy) {
      next_wall = next_z;
      ds = dz;
      next_direction[0] = 0;
      next_direction[1] = 0;
    } else {
      // special cases: at least two of the smallest values are equal
      if (dx == dy && dx < dz) {
        // it does not matter which values we pick, they will be the same
        next_wall = next_x;
        ds = dx;
        next_direction[2] = 0;
      } else if (dx == dz && dx < dy) {
        next_wall = next_x;
        ds = dx;
        next_direction[1] = 0;
      } else if (dy == dz && dy < dx) {
        next_wall = next_y;
        ds = dy;
        next_direction[0] = 0;
      } else {
        // all values are equal, we sit on a corner of the box
        next_wall = next_x;
        ds = dx;
      }
    }

    // ds contains the squared norm, take the square root
    ds = sqrt(ds);

    // find the next index
    // next_direction stores the index of the cell AT THE SAME LEVEL, relative
    // w.r.t. the current cell
    // to find the index of the next cell, we need to
    // - find that cell or the lowest lying parent cell
    // - if that cell has children: find the child that contains the
    //   intersection point
    // If the next cell is outside the box, we set the next index to
    // AMRGRID_MAXINDEX. However, if the box is periodic in that dimension, we
    // need to figure out what the index of the next cell at the other side of
    // the boundary is.
    unsigned long new_index =
        _grid.get_neighbour(index, next_direction, next_wall);
    if (index == AMRGRID_MAXKEY) {
      // we are outside the box. Check if periodic boundaries should be applied
      // for now, we assume we move in one direction only, because moving in
      // more directions at the same time requires extra checks...
      unsigned int num_0 = (next_direction.x() != 0);
      num_0 += (next_direction.y() != 0);
      num_0 += (next_direction.z() != 0);
      if (num_0 != 2) {
        error("Not supported yet!");
      }
      CoordinateVector<> new_position = next_wall;
      if (next_direction.x() != 0 && _periodic.x()) {
        new_position[0] -= next_direction.x() * _box.get_sides().x();
        periodic_correction[0] = -next_direction.x() * _box.get_sides().x();
        new_index = _grid.get_first_key(next_direction, new_position);
      }
      if (next_direction.y() != 0 && _periodic.y()) {
        new_position[1] -= next_direction.y() * _box.get_sides().y();
        periodic_correction[1] = -next_direction.y() * _box.get_sides().y();
        new_index = _grid.get_first_key(next_direction, new_position);
      }
      if (next_direction.z() != 0 && _periodic.z()) {
        new_position[2] -= next_direction.z() * _box.get_sides().z();
        periodic_correction[2] = -next_direction.z() * _box.get_sides().z();
        new_index = _grid.get_first_key(next_direction, new_position);
      }
    }
    index = new_index;

    return next_wall;
  }

  /**
   * @brief Let the given Photon travel through the density grid until the given
   * optical depth is reached.
   *
   * @param photon Photon.
   * @param optical_depth Optical depth the photon should travel in total
   * (dimensionless).
   * @return True if the Photon is still in the box after the optical depth has
   * been reached, false otherwise.
   */
  virtual bool interact(Photon &photon, double optical_depth) {
    CoordinateVector<> photon_origin = photon.get_position();
    CoordinateVector<> photon_direction = photon.get_direction();

    // find out in which cell the photon is currently hiding
    unsigned long index = get_cell_index(photon_origin);

    // while the photon has not exceeded the optical depth and is still in the
    // box
    while (index != AMRGRID_MAXKEY && optical_depth > 0.) {
      Box cell = _grid.get_geometry(index);

      double ds = 0.;
      unsigned long old_index = index;
      CoordinateVector<> periodic_correction;
      CoordinateVector<> next_wall =
          get_wall_intersection(photon_origin, photon_direction, cell, index,
                                ds, periodic_correction);

      // get the optical depth of the path from the current photon location to
      // the
      // cell wall, update S
      DensityValues &density = get_cell_values(old_index);

      // Helium abundance. Should be a parameter.
      double tau = get_optical_depth(ds, density, photon);
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
        // apply periodic boundaries if necessary
        photon_origin += periodic_correction;
      }

      // ds is now the actual distance travelled in the cell
      // update contributions to mean intensity integrals
      update_integrals(ds, density, photon);
    }

    photon.set_position(photon_origin);

    return index != AMRGRID_MAXKEY;
  }

  /**
   * @brief Increment the iterator index.
   *
   * In this case, the index encodes a lot of extra information and we cannot
   * simply increment it by 1.
   *
   * @param index Index to increment.
   */
  virtual void increase_index(unsigned long &index) {
    index = _grid.get_next_key(index);
  }

  /**
   * @brief Get an iterator to the first cell in the grid.
   *
   * @return Iterator to the first cell in the grid.
   */
  virtual DensityGrid::iterator begin() {
    return iterator(_grid.get_first_key(), *this);
  }

  /**
   * @brief Get an iterator to the last cell in the grid.
   *
   * @return Iterator to the last cell in the grid.
   */
  virtual DensityGrid::iterator end() {
    return iterator(_grid.get_max_key(), *this);
  }

  /**
   * @brief Print the grid to the given stream for visual inspection.
   *
   * @param stream std::ostream to write to.
   */
  inline void print(std::ostream &stream) { _grid.print(stream); }
};

#endif // AMRDENSITYGRID_HPP
