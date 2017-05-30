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
  AMRGrid< unsigned long > _grid;

  /*! @brief Convenient cell list used for faster cell indexing. */
  std::vector< AMRGridCell< unsigned long > * > _cells;

  /*! @brief AMRRefinementScheme used to refine cells. */
  AMRRefinementScheme *_refinement_scheme;

  /*! @brief Refinement interval: number of grid resets before the grid is
      refined for the first time. */
  unsigned char _refinement_interval;

  /*! @brief Number of times the reset() method was called. */
  unsigned char _reset_count;

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
   * @brief Check if the cell with the given index should be refined, using the
   * given AMRRefinementScheme. Apply the refinement if necessary, using the
   * given DensityFunction to recalculate the densities in the refined cells.
   *
   * @param refinement_scheme AMRRefinementScheme.
   * @param index Index of the cell that should be checked.
   * @param density_function DensityFunction used to recalculate densities for
   * refined cells.
   */
  inline void refine_cell(AMRRefinementScheme &refinement_scheme,
                          unsigned long index,
                          DensityFunction &density_function) {
    AMRGridCell< unsigned long > &cell = *_cells[index];
    DensityGrid::iterator it(index, *this);
    unsigned char level = cell.get_level();

    if (refinement_scheme.refine(level, it)) {

      // make a local copy of the contents of the cell
      // the contents of the new cells is initialized to these values
      double old_ionic_fractions[NUMBER_OF_IONNAMES];
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        old_ionic_fractions[i] = _ionic_fraction[i][index];
      }
      double old_temperature = _temperature[index];
      double old_hydrogen_reemission_probability =
          _hydrogen_reemission_probability[index];
      double old_helium_reemission_probability[4];
      for (int i = 0; i < 4; ++i) {
        old_helium_reemission_probability[i] =
            _helium_reemission_probability[i][index];
      }
      // we do not copy the mean intensity integrals from the old cell, as these
      // will be reset before the refined cells are used
      double old_neutral_fraction_H_old = _neutral_fraction_H_old[index];
      // we will not copy the heating terms for the same reasons
      // nor the EmissivityValues
      // nor the Lock, since that has to be unique

      cell.create_all_cells(level, level + 1);
      for (unsigned int ic = 0; ic < 8; ++ic) {
        AMRChildPosition child = static_cast< AMRChildPosition >(ic);
        AMRGridCell< unsigned long > *childcell = cell.get_child(child);
        DensityValues funcvalue = density_function(childcell->get_midpoint());
        // the first child replaces the old cell
        // the other children are added to the end of the internal lists
        if (ic == 0) {
          _number_density[index] = funcvalue.get_number_density();
          for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
            _ionic_fraction[i][index] = old_ionic_fractions[i];
          }
          _temperature[index] = old_temperature;
          _hydrogen_reemission_probability[index] =
              old_hydrogen_reemission_probability;
          for (int i = 0; i < 4; ++i) {
            _helium_reemission_probability[i][index] =
                old_helium_reemission_probability[i];
          }
          for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
            _mean_intensity[i][index] = 0.;
          }
          _mean_intensity_H_old[index] = 0.;
          _neutral_fraction_H_old[index] = old_neutral_fraction_H_old;
          _heating_H[index] = 0.;
          _heating_He[index] = 0.;
          _emissivities[index] = nullptr;
          _cells[index] = childcell;
          childcell->value() = index;
        } else {
          _number_density.push_back(funcvalue.get_number_density());
          for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
            _ionic_fraction[i].push_back(old_ionic_fractions[i]);
          }
          _temperature.push_back(old_temperature);
          _hydrogen_reemission_probability.push_back(
              old_hydrogen_reemission_probability);
          for (int i = 0; i < 4; ++i) {
            _helium_reemission_probability[i].push_back(
                old_helium_reemission_probability[i]);
          }
          for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
            _mean_intensity[i].push_back(0.);
          }
          _mean_intensity_H_old.push_back(0.);
          _neutral_fraction_H_old.push_back(old_neutral_fraction_H_old);
          _heating_H.push_back(0.);
          _heating_He.push_back(0.);
          _emissivities.push_back(nullptr);
#ifndef USE_LOCKFREE
          _lock.push_back(Lock());
#endif
          _cells.push_back(childcell);
          childcell->value() = _cells.size() - 1;
        }
        // recursively refine further
        refine_cell(refinement_scheme, childcell->value(), density_function);
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
   * @param refinement_interval Number of grid resets before the grid is
   * refined for the first time.
   * @param periodic Periodicity flags.
   * @param hydro Hydro flag.
   * @param log Log to write logging info to.
   */
  inline AMRDensityGrid(
      Box box, CoordinateVector< int > ncell, DensityFunction &density_function,
      AMRRefinementScheme *refinement_scheme = nullptr,
      unsigned char refinement_interval = 5,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      bool hydro = false, Log *log = nullptr)
      : DensityGrid(density_function, box, periodic, hydro, log),
        _refinement_scheme(refinement_scheme),
        _refinement_interval(refinement_interval), _reset_count(0) {

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
    _grid = AMRGrid< unsigned long >(box, nblock);

    // find out how many cells each block should have at the lowest level
    // this is just the power in power_of_2
    unsigned char level = 0;
    while (power_of_2 > 1) {
      power_of_2 >>= 1;
      ++level;
    }
    _grid.create_all_cells(level);

    // construct the cell list and set the contents of the cells to the correct
    // index values
    _cells.resize(_grid.get_number_of_cells());
    unsigned int index = 0;
    unsigned long key = _grid.get_first_key();
    while (key != _grid.get_max_key()) {
      _cells[index] = &_grid[key];
      _cells[index]->value() = index;
      ++index;
      key = _grid.get_next_key(key);
    }

    allocate_memory(_grid.get_number_of_cells());

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
            params.get_value< unsigned char >("densitygrid:refinement_interval",
                                              5),
            params.get_value< CoordinateVector< bool > >(
                "densitygrid:periodicity", CoordinateVector< bool >(false)),
            params.get_value< bool >("hydro:active", false), log) {}

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
   * @brief Initialize the cells of the grid.
   *
   * @param block Block that should be initialized by this MPI process.
   */
  virtual void initialize(std::pair< unsigned long, unsigned long > &block) {
    DensityGrid::initialize(block);
    DensityGrid::initialize(block, _density_function);

    // apply mesh refinement
    if (_refinement_scheme) {
      if (_log) {
        _log->write_status("Applying refinement.");
      }

      /// WE HAVE TO DO SPECIAL THINGS FOR MPI HERE!
      // we only refine the cells that were already in the grid
      // the new cells that are added during refinement are recursively refined
      // within the refinement routine
      const unsigned int cell2size = _cells.size();
      for (unsigned int i = 0; i < cell2size; ++i) {
        refine_cell(*_refinement_scheme, i, _density_function);
      }

      if (_log) {
        _log->write_status("Number of cells after refinement: ",
                           _grid.get_number_of_cells());
      }
    }

    // finalize grid: set neighbour relations
    _grid.set_ngbs(_periodic);

    // make sure all values are correctly initialized (also in refined cells)
    // the refinement procedure itself only reads the density from the density
    // function, as it is also used in reset_grid()
    // at this point, we want to read all values from the density function,
    // since it also contains the initial temperature etc.
    DensityGrid::initialize(block, _density_function);
  }

  /**
   * @brief Reset the mean intensity counters, update the reemission
   * probabilities, and reapply the refinement scheme to all cells.
   */
  virtual void reset_grid() {
    if (_log) {
      _log->write_status("Resetting grid...");
    }

    ++_reset_count;
    // we only refine the cells that are already present
    // newly added refined cells are recursively refined within the refinement
    // routine
    if (_refinement_scheme != nullptr && _reset_count >= _refinement_interval) {
      const unsigned int cells2size = _cells.size();
      for (unsigned int i = 0; i < cells2size; ++i) {
        refine_cell(*_refinement_scheme, i, _density_function);
      }

      if (_log) {
        _log->write_status("Number of cells after reset: ",
                           _grid.get_number_of_cells(), ".");
      }

      // reset the ngbs
      _grid.set_ngbs(_periodic);
    }

    // make sure all cells are correctly reset (also the new ones, if any)
    DensityGrid::reset_grid();
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
    return _grid.get_cell(position);
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
                        AMRGridCell< unsigned long > *&cell, double &ds,
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

    AMRGridCell< unsigned long > *next_cell = cell->get_ngb(ngbposition);
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
   * @return DensityGrid::iterator pointing to the cell the photon was last in,
   * or DensityGrid::end() if the photon left the box.
   */
  virtual DensityGrid::iterator interact(Photon &photon, double optical_depth) {
    CoordinateVector<> photon_origin = photon.get_position();
    CoordinateVector<> photon_direction = photon.get_direction();

    unsigned long index = get_cell_index(photon_origin);
    AMRGridCell< unsigned long > *current_cell = _cells[index];

    // while the photon has not exceeded the optical depth and is still in the
    // box
    DensityGrid::iterator last_cell = end();
    while (current_cell != nullptr && optical_depth > 0.) {
      Box cell = current_cell->get_geometry();

      double ds = 0.;
      AMRGridCell< unsigned long > *old_cell = current_cell;
      CoordinateVector<> periodic_correction;
      CoordinateVector<> next_wall =
          get_wall_intersection(photon_origin, photon_direction, cell,
                                current_cell, ds, periodic_correction);

      // get the optical depth of the path from the current photon location to
      // the cell wall, update S
      DensityGrid::iterator it(old_cell->value(), *this);
      last_cell = it;

      // Helium abundance. Should be a parameter.
      double tau = get_optical_depth(ds, it, photon);
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
      update_integrals(ds, it, photon);
    }

    photon.set_position(photon_origin);

    if (current_cell == nullptr) {
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
  virtual double get_total_emission(CoordinateVector<> origin,
                                    CoordinateVector<> direction,
                                    EmissionLine line) {
    double S = 0.;

    unsigned long index = get_cell_index(origin);
    AMRGridCell< unsigned long > *current_cell = _cells[index];

    while (current_cell != nullptr) {
      Box cell = current_cell->get_geometry();

      double ds = 0.;
      AMRGridCell< unsigned long > *old_cell = current_cell;
      CoordinateVector<> periodic_correction;
      CoordinateVector<> next_wall = get_wall_intersection(
          origin, direction, cell, current_cell, ds, periodic_correction);

      DensityGrid::iterator it(old_cell->value(), *this);

      origin = next_wall;

      S += it.get_emissivities()->get_emissivity(line);

      if (periodic_correction.norm2() > 0.) {
        break;
      }
    }

    return S;
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
   * @brief Get the neighbours of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return std::vector containing the neighbours of the cell.
   */
  virtual std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                                   CoordinateVector<>, double > >
  get_neighbours(unsigned long index) {
    std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                             CoordinateVector<>, double > >
        ngbs;

    return ngbs;
  }

  /**
   * @brief Print the grid to the given stream for visual inspection.
   *
   * @param stream std::ostream to write to.
   */
  inline void print(std::ostream &stream) { _grid.print(stream); }
};

#endif // AMRDENSITYGRID_HPP
