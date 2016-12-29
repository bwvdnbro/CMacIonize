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
#include "Abundances.hpp"
#include "DensityGrid.hpp"
#include "ParameterFile.hpp"

#include <cfloat>
#include <ostream>
#include <vector>

/**
 * @brief AMR density grid.
 */
class AMRDensityGrid : public DensityGrid {
private:
  /*! @brief AMRGrid used as grid. */
  AMRGrid< DensityValues > _grid;

  /*! @brief Convenient cell list used for faster cell indexing. */
  std::vector< AMRGridCell< DensityValues > * > _cells;

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
  inline static int get_largest_odd_factor(int number) {
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
  inline static int get_power_of_two(int number) {
    return number / get_largest_odd_factor(number);
  }

  /**
   * @brief Refine the given cell using the given refinement scheme to decide
   * how deep the refinement should be.
   *
   * @param refinement_scheme AMRRefinementScheme.
   * @param cell AMRGridCell to refine.
   * @param density_function DensityFunction that sets the value of the density
   * in the cell.
   */
  inline void refine_cell(AMRRefinementScheme &refinement_scheme,
                          AMRGridCell< DensityValues > &cell,
                          DensityFunction &density_function) {
    unsigned char level = cell.get_level();
    CoordinateVector<> midpoint = cell.get_midpoint();
    double volume = cell.get_volume();
    DensityValues values = cell.value();
    if (refinement_scheme.refine(level, midpoint, volume, values)) {
      cell.create_all_cells(level, level + 1);
      for (unsigned int ic = 0; ic < 8; ++ic) {
        AMRChildPosition child = static_cast< AMRChildPosition >(ic);
        AMRGridCell< DensityValues > *childcell = cell.get_child(child);
        DensityValues &childvalues = childcell->value();
        // we only set the density based on the density function, as all other
        // variables are only initial conditions
        childvalues.set_total_density(
            density_function(childcell->get_midpoint()).get_total_density());
        for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
          IonName ion = static_cast< IonName >(i);
          childvalues.set_ionic_fraction(ion, values.get_ionic_fraction(ion));
        }
        childvalues.set_temperature(values.get_temperature());
        childvalues.set_old_neutral_fraction_H(
            values.get_old_neutral_fraction_H());
        // set reemission probabilities
        set_reemission_probabilities(childvalues.get_temperature(),
                                     childvalues);
        // recursively refine further
        refine_cell(refinement_scheme, *childcell, density_function);
      }
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the grid.
   * @param ncell Number of cells in the low resolution grid.
   * @param density_function DensityFunction that defines the density field.
   * @param refinement_scheme Refinement scheme used to refine cells. Memory
   * management for this pointer is taken over by this class.
   * @param periodic Periodicity flags.
   * @param log Log to write logging info to.
   */
  inline AMRDensityGrid(
      Box box, CoordinateVector< int > ncell, DensityFunction &density_function,
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

    // construct the cell list
    _cells.resize(_grid.get_number_of_cells());
    unsigned int index = 0;
    unsigned long key = _grid.get_first_key();
    while (key != _grid.get_max_key()) {
      _cells[index] = &_grid[key];
      ++index;
      key = _grid.get_next_key(key);
    }

    if (_log) {
      int levelint = level;
      _log->write_status("Created AMRGrid with ", nblock.x(), "x", nblock.y(),
                         "x", nblock.z(), " top level blocks, going ", levelint,
                         " levels deep, in a box with origin [",
                         _box.get_anchor().x(), " m, ", _box.get_anchor().y(),
                         " m, ", _box.get_anchor().z(), " m], and sides [",
                         _box.get_sides().x(), " m, ", _box.get_sides().y(),
                         " m, ", _box.get_sides().z(), " m].");
    }

    initialize(density_function);

    // apply mesh refinement
    if (_refinement_scheme) {
      if (_log) {
        _log->write_status("Applying refinement.");
      }

      for (unsigned int i = 0; i < _cells.size(); ++i) {
        refine_cell(*_refinement_scheme, *_cells[i], _density_function);
      }

      for (unsigned int i = 0; i < _cells.size(); ++i) {
        AMRGridCell< DensityValues > *parent = _cells[i]->get_parent();
        if (parent && parent->get_child(AMRCHILDPOSITION_LFB) == _cells[i]) {
          // we only check the first child, since all other children have the
          // same parent
          CoordinateVector<> midpoints[8];
          double volumes[8];
          DensityValues cells[8];
          bool lowest_level = true;
          for (int ichild = 0; ichild < 8; ++ichild) {
            AMRGridCell< DensityValues > *child =
                parent->get_child(static_cast< AMRChildPosition >(ichild));
            if (child->is_single_cell()) {
              midpoints[ichild] = child->get_midpoint();
              volumes[ichild] = child->get_volume();
              cells[ichild] = child->value();
            } else {
              lowest_level = false;
            }
          }
          if (lowest_level &&
              _refinement_scheme->coarsen(parent->get_level() + 1, midpoints,
                                          volumes, cells)) {
            cmac_status("Coarsening needed");
          }
        }
      }

      if (_log) {
        _log->write_status(
            "Number of cells after refinement: ", _grid.get_number_of_cells());
      }
    }

    // finalize grid: set neighbour relations
    _grid.set_ngbs(_periodic);
    // reconstruct the cell list, it might have changed due to refinement
    _cells.resize(_grid.get_number_of_cells());
    index = 0;
    key = _grid.get_first_key();
    while (key != _grid.get_max_key()) {
      _cells[index] = &_grid[key];
      ++index;
      key = _grid.get_next_key(key);
    }

    // make sure all values are correctly initialized (also in refined cells)
    initialize(density_function);
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read.
   * @param density_function DensityFunction used to set the densities in each
   * cell.
   * @param log Log to write log messages to.
   */
  inline AMRDensityGrid(ParameterFile &params,
                        DensityFunction &density_function, Log *log)
      : AMRDensityGrid(
            Box(params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_anchor", "[0. m, 0. m, 0. m]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_sides", "[1. m, 1. m, 1. m]")),
            params.get_value< CoordinateVector< int > >(
                "densitygrid:ncell", CoordinateVector< int >(64)),
            density_function, AMRRefinementSchemeFactory::generate(params, log),
            params.get_value< CoordinateVector< bool > >(
                "densitygrid:periodicity", CoordinateVector< bool >(false)),
            log) {}

  /**
   * @brief Destructor.
   *
   * Deletes the AMRRefinementScheme (if present).
   */
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
    if (_log) {
      _log->write_status("Resetting grid...");
    }
    for (unsigned int i = 0; i < _cells.size(); ++i) {
      _cells[i]->value().reset_mean_intensities();
      if (_refinement_scheme) {
        refine_cell(*_refinement_scheme, *_cells[i], _density_function);
      }
    }
    if (_log) {
      _log->write_status(
          "Number of cells after reset: ", _grid.get_number_of_cells(), ".");
    }
    // reset the ngbs
    _grid.set_ngbs(_periodic);
    // reset the cell list
    _cells.resize(_grid.get_number_of_cells());
    unsigned int index = 0;
    unsigned long key = _grid.get_first_key();
    while (key != _grid.get_max_key()) {
      _cells[index] = &_grid[key];
      ++index;
      key = _grid.get_next_key(key);
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
  virtual unsigned int get_number_of_cells() const { return _cells.size(); }

  /**
   * @brief Get the index of the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return Index of the cell containing that position.
   */
  virtual unsigned long get_cell_index(CoordinateVector<> position) const {
    return _grid.get_key(position);
  }

  /**
   * @brief Get the midpoint of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Midpoint of that cell (in m).
   */
  virtual CoordinateVector<> get_cell_midpoint(unsigned long index) const {
    return _cells[index]->get_midpoint();
  }

  /**
   * @brief Get the values stored in the cell with the given index.
   *
   * @param index Index of a cell.
   * @return DensityValues stored in that cell.
   */
  virtual DensityValues &get_cell_values(unsigned long index) const {
    return _cells[index]->value();
  }

  /**
   * @brief Get the values stored in the cell which contains the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return DensityValues of the cell containing that position (in SI units).
   */
  virtual DensityValues &get_cell_values(CoordinateVector<> position) const {
    return _grid.get_cell(position);
  }

  /**
   * @brief Get the volume of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Volume of that cell (in m^3).
   */
  virtual double get_cell_volume(unsigned long index) const {
    return _cells[index]->get_volume();
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
   * @param box Geometry of the cell in which the photon currently resides.
   * @param cell Pointer to the current cell, is updated to the
   * neighbouring cell that is next in the algorithm.
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
                        CoordinateVector<> &photon_direction, Box &box,
                        AMRGridCell< DensityValues > *&cell, double &ds,
                        CoordinateVector<> &periodic_correction) {
    CoordinateVector<> cell_bottom_anchor = box.get_anchor();
    CoordinateVector<> cell_top_anchor = box.get_top_anchor();

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
    AMRNgbPosition ngbposition;
    if (dx < dy && dx < dz) {
      next_wall = next_x;
      ds = dx;
      next_direction[1] = 0;
      next_direction[2] = 0;
      if (next_direction[0] < 0.) {
        ngbposition = AMRNGBPOSITION_LEFT;
      } else {
        ngbposition = AMRNGBPOSITION_RIGHT;
      }
    } else if (dy < dx && dy < dz) {
      next_wall = next_y;
      ds = dy;
      next_direction[0] = 0;
      next_direction[2] = 0;
      if (next_direction[1] < 0.) {
        ngbposition = AMRNGBPOSITION_FRONT;
      } else {
        ngbposition = AMRNGBPOSITION_BACK;
      }
    } else if (dz < dx && dz < dy) {
      next_wall = next_z;
      ds = dz;
      next_direction[0] = 0;
      next_direction[1] = 0;
      if (next_direction[2] < 0.) {
        ngbposition = AMRNGBPOSITION_BOTTOM;
      } else {
        ngbposition = AMRNGBPOSITION_TOP;
      }
    } else {
      // special cases: at least two of the smallest values are equal
      // find out which values are equal, and pick one of the two cells as next
      // cell
      if (dx == dy || dx == dz) {
        next_wall = next_x;
        ds = dx;
        next_direction[1] = 0;
        next_direction[2] = 0;
        if (next_direction[0] < 0.) {
          ngbposition = AMRNGBPOSITION_LEFT;
        } else {
          ngbposition = AMRNGBPOSITION_RIGHT;
        }
      } else {
        next_wall = next_y;
        ds = dy;
        next_direction[0] = 0;
        next_direction[2] = 0;
        if (next_direction[1] < 0.) {
          ngbposition = AMRNGBPOSITION_FRONT;
        } else {
          ngbposition = AMRNGBPOSITION_BACK;
        }
      }
    }

    // ds contains the squared norm, take the square root
    ds = sqrt(ds);

    AMRGridCell< DensityValues > *next_cell = cell->get_ngb(ngbposition);
    if (next_cell != nullptr) {
      // calculate periodic boundary corrections (if any)
      if (_periodic.x()) {
        CoordinateVector<> nm = next_cell->get_geometry().get_anchor();
        if (next_direction[0] > 0. && nm.x() < cell_bottom_anchor.x()) {
          periodic_correction[0] = -_box.get_sides().x();
        } else if (next_direction[0] < 0. && nm.x() > cell_bottom_anchor.x()) {
          periodic_correction[0] = _box.get_sides().x();
        }
      }
      if (_periodic.y()) {
        CoordinateVector<> nm = next_cell->get_geometry().get_anchor();
        if (next_direction[1] > 0. && nm.y() < cell_bottom_anchor.y()) {
          periodic_correction[1] = -_box.get_sides().y();
        } else if (next_direction[1] < 0. && nm.y() > cell_bottom_anchor.y()) {
          periodic_correction[1] = _box.get_sides().y();
        }
      }
      if (_periodic.z()) {
        CoordinateVector<> nm = next_cell->get_geometry().get_anchor();
        if (next_direction[2] > 0. && nm.z() < cell_bottom_anchor.z()) {
          periodic_correction[2] = -_box.get_sides().z();
        } else if (next_direction[2] < 0. && nm.z() > cell_bottom_anchor.z()) {
          periodic_correction[2] = _box.get_sides().z();
        }
      }
      // find the child cell containing the new position
      while (!next_cell->is_single_cell()) {
        next_cell = next_cell->get_child(next_wall);
      }
    }
    cell = next_cell;

    return next_wall;
  }

  /**
   * @brief Let the given Photon travel through the density grid until the given
   * optical depth is reached.
   *
   * @param photon Photon.
   * @param optical_depth Optical depth the photon should travel in total
   * (dimensionless).
   * @return A pointer to the values of the last cell the photon was in, nullptr
   * if the photon left the box.
   */
  virtual DensityValues *interact(Photon &photon, double optical_depth) {
    CoordinateVector<> photon_origin = photon.get_position();
    CoordinateVector<> photon_direction = photon.get_direction();

    unsigned long index = get_cell_index(photon_origin);
    AMRGridCell< DensityValues > *current_cell = &_grid[index];

    // while the photon has not exceeded the optical depth and is still in the
    // box
    DensityValues *last_cell = nullptr;
    while (current_cell != nullptr && optical_depth > 0.) {
      Box cell = current_cell->get_geometry();

      double ds = 0.;
      AMRGridCell< DensityValues > *old_cell = current_cell;
      CoordinateVector<> periodic_correction;
      CoordinateVector<> next_wall =
          get_wall_intersection(photon_origin, photon_direction, cell,
                                current_cell, ds, periodic_correction);

      // get the optical depth of the path from the current photon location to
      // the
      // cell wall, update S
      DensityValues &density = old_cell->value();
      last_cell = &density;

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

    if (current_cell == nullptr) {
      last_cell = nullptr;
    }

    return last_cell;
  }

  /**
   * @brief Increment the iterator index.
   *
   * In this case, the index encodes a lot of extra information and we cannot
   * simply increment it by 1.
   *
   * @param index Index to increment.
   */
  virtual void increase_index(unsigned long &index) { ++index; }

  /**
   * @brief Get an iterator to the first cell in the grid.
   *
   * @return Iterator to the first cell in the grid.
   */
  virtual DensityGrid::iterator begin() { return iterator(0, *this); }

  /**
   * @brief Get an iterator to the last cell in the grid.
   *
   * @return Iterator to the last cell in the grid.
   */
  virtual DensityGrid::iterator end() { return iterator(_cells.size(), *this); }

  /**
   * @brief Get begin and end iterators to a chunk of the grid with given begin
   * and end fractions.
   *
   * @param begin Fraction of the total grid where we want the chunk to begin.
   * @param end Fraction of the total grid where we want the chunk to end.
   * @return std::pair of iterators pointing to the begin and end of the chunk.
   */
  virtual std::pair< DensityGrid::iterator, DensityGrid::iterator >
  get_chunk(double begin, double end) {
    unsigned int npart = _cells.size();
    unsigned int ibegin = begin * npart;
    unsigned int iend = end * npart;
    return std::make_pair(iterator(ibegin, *this), iterator(iend, *this));
  }

  /**
   * @brief Print the grid to the given stream for visual inspection.
   *
   * @param stream std::ostream to write to.
   */
  inline void print(std::ostream &stream) { _grid.print(stream); }
};

#endif // AMRDENSITYGRID_HPP
