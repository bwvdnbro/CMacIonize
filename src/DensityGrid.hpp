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
 * @brief General interface for density grids.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRID_HPP
#define DENSITYGRID_HPP

#include "Abundances.hpp"
#include "Box.hpp"
#include "Cell.hpp"
#include "Configuration.hpp"
#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "DensityValues.hpp"
#include "EmissivityValues.hpp"
#include "HydroVariables.hpp"
#include "IonizationVariables.hpp"
#include "Lock.hpp"
#include "Log.hpp"
#include "Photon.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "Timer.hpp"
#include "UnitConverter.hpp"
#include "WorkDistributor.hpp"

#ifdef USE_LOCKFREE
#include "Atomic.hpp"
#endif

#include <cmath>
#include <tuple>

/*! @brief Size of the variables storing cell indices; this should be big enough
 *  to store at least the number of cells. */
typedef size_t cellsize_t;

/**
 * @brief General interface for density grids.
 */
class DensityGrid {
public:
  class iterator;

protected:
  /*! @brief Box containing the grid. */
  const Box<> _box;

  /*! @brief Periodicity flags. */
  const CoordinateVector< bool > _periodicity_flags;

  /*! @brief Ionization energy of hydrogen (in Hz). */
  const double _ionization_energy_H;

  /*! @brief Ionization energy of helium (in Hz). */
  const double _ionization_energy_He;

  /*! @brief Flag indicating whether hydro is active or not. */
  const bool _has_hydro;

  /*! @brief Ionization calculation variables. */
  std::vector< IonizationVariables > _ionization_variables;

  /// hydro

  /*! @brief Hydrodynamic variables. */
  std::vector< HydroVariables > _hydro_variables;

  /// end hydro

  /*! @brief EmissivityValues for the cells. */
  std::vector< EmissivityValues * > _emissivities;

#ifndef USE_LOCKFREE
  /*! @brief Locks to ensure safe write access to the cell data. */
  std::vector< Lock > _lock;
#endif

  /*! @brief Log to write log messages to. */
  Log *_log;

  /**
   * @brief Get the optical depth for a photon travelling the given path in the
   * given cell.
   *
   * @param ds Path length the photon traverses (in m).
   * @param ionization_variables IonizationVariables of the cell.
   * @param photon Photon.
   * @return Optical depth.
   */
  inline static double
  get_optical_depth(double ds, const IonizationVariables &ionization_variables,
                    const Photon &photon) {
    return ds * ionization_variables.get_number_density() *
           (photon.get_cross_section(ION_H_n) *
                ionization_variables.get_ionic_fraction(ION_H_n) +
            photon.get_cross_section_He_corr() *
                ionization_variables.get_ionic_fraction(ION_He_n));
  }

  /**
   * @brief Update the contributions to the mean intensity integrals due to the
   * given photon travelling the given path length in the given cell.
   *
   * @param ds Path length the photon traverses (in m).
   * @param cell DensityValues of the cell the photon travels in.
   * @param photon Photon.
   */
  inline void update_integrals(double ds, DensityGrid::iterator &cell,
                               const Photon &photon) const {

    IonizationVariables &ionization_variables = cell.get_ionization_variables();
    if (ionization_variables.get_number_density() > 0.) {
      // we tried speeding things up by using lock-free addition, but it turns
      // out that the overhead caused by doing this is larger than the overhead
      // of using a single lock
      // this is mainly because we have to do a large number of additions
      // to minimize collisions (two threads trying to access the same cell at
      // the same time), we first calculate all terms that need to be added, and
      // then lock the cell and do all additions as fast as possible
      double dmean_intensity[NUMBER_OF_IONNAMES];
      const double dsw = ds * photon.get_weight();
      for (int ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        dmean_intensity[ion] = dsw * photon.get_cross_section(ion);
      }
      const double dheating_H = dmean_intensity[ION_H_n] *
                                (photon.get_energy() - _ionization_energy_H);
      const double dheating_He = dmean_intensity[ION_He_n] *
                                 (photon.get_energy() - _ionization_energy_He);
#ifndef USE_LOCKFREE
      cell.lock();
#endif
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        ionization_variables.increase_mean_intensity(ion, dmean_intensity[ion]);
      }
      ionization_variables.increase_heating(HEATINGTERM_H, dheating_H);
      ionization_variables.increase_heating(HEATINGTERM_He, dheating_He);
      Tracker *tracker = ionization_variables.get_tracker();
      if (tracker != nullptr) {
        tracker->count_photon(photon);
      }
#ifndef USE_LOCKFREE
      cell.unlock();
#endif
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * The constructor of a DensityGrid implementation should not initialize the
   * entire grid, but should only initialize the internal variables that are
   * necessary to perform the initialization at a later stage. We postpone the
   * computationally expensive part of the initialization to the initialize()
   * routine.
   *
   * @param box Box containing the grid.
   * @param periodic Periodicity flags.
   * @param hydro Hydro flag.
   * @param log Log to write log messages to.
   */
  DensityGrid(
      Box<> box,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      bool hydro = false, Log *log = nullptr)
      : _box(box), _periodicity_flags(periodic),
        _ionization_energy_H(
            UnitConverter::to_SI< QUANTITY_FREQUENCY >(13.6, "eV")),
        _ionization_energy_He(
            UnitConverter::to_SI< QUANTITY_FREQUENCY >(24.6, "eV")),
        _has_hydro(hydro), _log(log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~DensityGrid() {}

  /**
   * @brief Allocate memory for the given number of cells.
   *
   * @param numcell Number of cells that will be stored in the grid.
   */
  inline void allocate_memory(cellsize_t numcell) {
    if (_log) {
      _log->write_status(
          "Allocating memory for ", numcell, " cells (",
          Utilities::human_readable_bytes(numcell * sizeof(DensityValues)),
          ")...");
    }
    // we allocate memory for the cells, so that --dry-run can already check the
    // available memory
    _ionization_variables.resize(numcell);
    _emissivities.resize(numcell, nullptr);
#ifndef USE_LOCKFREE
    _lock.resize(numcell);
#endif

    if (_log) {
      _log->write_status("Done allocating memory.");
    }
  }

  /**
   * @brief Routine that does the actual initialization of the grid.
   *
   * This routine should do all the computationally intensive work that needs to
   * be done to initialize the grid. This work should not be done in the
   * constructor.
   *
   * @param block Block that should be initialized by this MPI process.
   * @param density_function DensityFunction to use.
   */
  virtual void initialize(std::pair< cellsize_t, cellsize_t > &block,
                          DensityFunction &density_function) {
    if (_has_hydro) {
      if (_log) {
        _log->write_status("Initializing hydro arrays...");
      }
      const cellsize_t numcell = get_number_of_cells();
      _hydro_variables.resize(numcell);
      if (_log) {
        _log->write_status("Done.");
      }
    }
  }

  /**
   * @brief Get the total number of cells in the grid.
   *
   * @return Number of cells in the grid.
   */
  virtual cellsize_t get_number_of_cells() const = 0;

  /**
   * @brief Get the Box containing the grid.
   *
   * @return Box containing the grid (in m).
   */
  inline const Box<> get_box() const { return _box; }

  /**
   * @brief Get the number of periodic boundaries of this grid.
   *
   * @return Number of periodic boundaries of this grid (between 0 and 3).
   */
  inline uint_fast8_t get_number_of_periodic_boundaries() const {
    uint_fast8_t numperiodic = 0;
    for (uint_fast8_t i = 0; i < 3; ++i) {
      numperiodic += _periodicity_flags[i];
    }
    return numperiodic;
  }

  /**
   * @brief Get the index of the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return Index of the cell containing that position.
   */
  virtual cellsize_t get_cell_index(CoordinateVector<> position) const = 0;

  /**
   * @brief Get the midpoint of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Midpoint of that cell (in m).
   */
  virtual CoordinateVector<> get_cell_midpoint(cellsize_t index) const = 0;

  /**
   * @brief Check if hydro is active.
   *
   * @return True if hydro is active.
   */
  inline bool has_hydro() const { return _has_hydro; }

  /**
   * @brief Get the neighbours of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return std::vector containing iterators to the neighbours, together with
   * the midpoint, surface normal and surface area of the boundary face between
   * the cell and this neighbour, and the relative position of the neighbour
   * w.r.t. the cell.
   */
  virtual std::vector<
      std::tuple< iterator, CoordinateVector<>, CoordinateVector<>, double,
                  CoordinateVector<> > >
  get_neighbours(cellsize_t index) = 0;

  /**
   * @brief Get the faces of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Faces of the cell.
   */
  virtual std::vector< Face > get_faces(cellsize_t index) const = 0;

  /**
   * @brief Get an iterator to the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return DensityGrid::iterator pointing to the cell containing that
   * position.
   */
  inline iterator get_cell(CoordinateVector<> position) {
    return iterator(get_cell_index(position), *this);
  }

  /**
   * @brief Get the volume of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Volume of that cell (in m^3).
   */
  virtual double get_cell_volume(cellsize_t index) const = 0;

  /**
   * @brief Get the total optical depth traversed by the given Photon until it
   * reaches the boundaries of the simulation box.
   *
   * @param photon Photon.
   * @return Total optical depth along the photon's path before it reaches the
   * boundaries of the simulation box.
   */
  virtual double integrate_optical_depth(const Photon &photon) = 0;

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
  virtual DensityGrid::iterator interact(Photon &photon,
                                         double optical_depth) = 0;

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
                                    EmissionLine line) = 0;

  /**
   * @brief Index increment used in the iterator.
   *
   * More sofisticated grids (like the AMR grid) might implement their own
   * version.
   *
   * @param index Index to increase.
   * @param increment Increment (default = 1).
   */
  virtual void increase_index(cellsize_t &index, cellsize_t increment = 1) {
    index += increment;
  }

  /**
   * @brief Iterator to loop over the cells in the grid.
   */
  class iterator : public Cell {
  private:
    /*! @brief Index of the cell the iterator is currently pointing to. */
    cellsize_t _index;

    /*! @brief Pointer to the DensityGrid over which we iterate (we cannot use a
     *  reference, since then things like it = it would not work). */
    DensityGrid *_grid;

  public:
    /**
     * @brief Constructor.
     *
     * @param index Index of the cell the iterator is currently pointing to.
     * @param grid DensityGrid over which we iterate.
     */
    inline iterator(cellsize_t index, DensityGrid &grid)
        : _index(index), _grid(&grid) {}

    /**
     * @brief Get the midpoint of the cell the iterator is pointing to.
     *
     * @return Cell midpoint (in m).
     */
    virtual CoordinateVector<> get_cell_midpoint() const {
      return _grid->get_cell_midpoint(_index);
    }

    /**
     * @brief Get the faces of the cell.
     *
     * @return Faces of the cell.
     */
    virtual std::vector< Face > get_faces() const {
      return _grid->get_faces(_index);
    }

    /**
     * @brief Reset the mean intensity counters for the cell the iterator is
     * currently pointing to.
     */
    inline void reset_mean_intensities() {
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        _grid->_ionization_variables[_index].set_mean_intensity(ion, 0.);
      }
      for (int_fast32_t name = 0; name < NUMBER_OF_HEATINGTERMS; ++name) {
        _grid->_ionization_variables[_index].set_heating(name, 0.);
      }
    }

    /**
     * @brief Dereference operator.
     *
     * Needed to work with the fancy MPI communication functions.
     *
     * @return Simply returns a reference to this iterator.
     */
    inline iterator &operator*() { return *this; }

    /**
     * @brief Get the mean intensity for the given ion for the cell the iterator
     * is currently pointing to.
     *
     * @param ion IonName.
     * @return Mean intensity (in m^3).
     */
    inline double get_mean_intensity(int_fast32_t ion) const {
      return _grid->_ionization_variables[_index].get_mean_intensity(ion);
    }

    /**
     * @brief Set the mean intensity for the given ion for the cell the iterator
     * is currently pointing to.
     *
     * @param mean_intensity New mean intensity value (in m^3).
     * @param ion IonName.
     */
    inline void set_mean_intensity(double mean_intensity, IonName ion) {
      _grid->_ionization_variables[_index].set_mean_intensity(ion,
                                                              mean_intensity);
    }

    /**
     * @brief Get read only access to the ionization variables stored in this
     * cell.
     *
     * @return Read only access to the ionization variables.
     */
    inline const IonizationVariables &get_ionization_variables() const {
      return _grid->_ionization_variables[_index];
    }

    /**
     * @brief Get read/write access to the ionization variables stored in this
     * cell.
     *
     * @return Read/write access to the ionization variables.
     */
    inline IonizationVariables &get_ionization_variables() {
      return _grid->_ionization_variables[_index];
    }

    /**
     * @brief Get read only access to the hydrodynamical variables stored in
     * this cell.
     *
     * @return Read only access to the hydrodynamical variables.
     */
    inline const HydroVariables &get_hydro_variables() const {
      return _grid->_hydro_variables[_index];
    }

    /**
     * @brief Get read/write access to the hydrodynamical variables stored in
     * this cell.
     *
     * @return Read/write access to the hydrodynamical variables.
     */
    inline HydroVariables &get_hydro_variables() {
      return _grid->_hydro_variables[_index];
    }

    /**
     * @brief Get the neighbours of the cell the iterator is currently pointing
     * to.
     *
     * @return std::vector containing iterators to the neighbours of the cell.
     */
    inline std::vector<
        std::tuple< iterator, CoordinateVector<>, CoordinateVector<>, double,
                    CoordinateVector<> > >
    get_neighbours() const {
      return _grid->get_neighbours(_index);
    }

    /**
     * @brief Get the EmissivityValues for the cell the iterator is currently
     * pointing to.
     *
     * @return EmissivityValues.
     */
    inline EmissivityValues *get_emissivities() const {
      return _grid->_emissivities[_index];
    }

    /**
     * @brief Set the EmissivityValues for the cell the iterator is currently
     * pointing to.
     *
     * @param emissivities EmissivityValues.
     */
    inline void set_emissivities(EmissivityValues *emissivities) {
      _grid->_emissivities[_index] = emissivities;
    }

#ifndef USE_LOCKFREE
    /**
     * @brief Lock the cell the iterator is pointing to.
     */
    inline void lock() { _grid->_lock[_index].lock(); }

    /**
     * @brief Unlock the cell the iterator is pointing to.
     */
    inline void unlock() { _grid->_lock[_index].unlock(); }
#endif

    /**
     * @brief Get the volume of the cell the iterator is pointing to.
     *
     * @return Volume of the cell (in m^3).
     */
    virtual double get_volume() const { return _grid->get_cell_volume(_index); }

    /**
     * @brief Increment operator.
     *
     * We only implemented the pre-increment version, since the post-increment
     * version creates a new object and is computationally more expensive.
     *
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator++() {
      _grid->increase_index(_index);
      return *this;
    }

    /**
     * @brief Increment operator.
     *
     * @param increment Increment to add.
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator+=(cellsize_t increment) {
      _grid->increase_index(_index, increment);
      return *this;
    }

    /**
     * @brief Free addition operator.
     *
     * @param increment Increment to add to the iterator.
     * @return Incremented iterator.
     */
    inline iterator operator+(cellsize_t increment) {
      iterator it(*this);
      it += increment;
      return it;
    }

    /**
     * @brief Get the index of the cell the iterator is currently pointing to.
     *
     * @return Index of the current cell.
     */
    inline cellsize_t get_index() const { return _index; }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators point to the same cell of the same grid.
     */
    inline bool operator==(iterator it) const {
      return (_grid == it._grid && _index == it._index);
    }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators do not point to the same cell of the same
     * grid.
     */
    inline bool operator!=(iterator it) const { return !(*this == it); }
  };

  /**
   * @brief Get an iterator to the first cell in the grid.
   *
   * @return Iterator to the first cell in the grid.
   */
  virtual iterator begin() = 0;

  /**
   * @brief Get an iterator to the last cell in the grid.
   *
   * @return Iterator to the last cell in the grid.
   */
  virtual iterator end() = 0;

  /**
   * @brief Get the velocity of the interface between the two given cells.
   *
   * For static grids, there is no need to implement this function, as the
   * interface velocity will always be zero in this case.
   *
   * @param left iterator to the left cell.
   * @param right iterator to the right cell.
   * @param interface_midpoint Coordinates of the midpoint of the interface
   * (in m).
   * @return Velocity of the interface (in m s^-1).
   */
  virtual CoordinateVector<>
  get_interface_velocity(const iterator left, const iterator right,
                         const CoordinateVector<> interface_midpoint) const {
    return CoordinateVector<>(0.);
  }

  /**
   * @brief Get begin and end iterators to a chunk of the grid with given begin
   * and end index.
   *
   * @param begin Start index.
   * @param end End index.
   * @return std::pair of iterators pointing to the begin and end of the chunk.
   */
  inline std::pair< iterator, iterator > get_chunk(cellsize_t begin,
                                                   cellsize_t end) {
    return std::make_pair(iterator(begin, *this), iterator(end, *this));
  }

  /**
   * @brief Functor class used to initialize the DensityGrid.
   */
  class DensityGridInitializationFunction {
  private:
    /*! @brief DensityFunction that sets the density for each cell in the grid.
     */
    DensityFunction &_function;

    /*! @brief Do we need to initialize hydro variables? */
    bool _hydro;

  public:
    /**
     * @brief Constructor.
     *
     * @param function DensityFunction that set the density for each cell in the
     * grid.
     * @param hydro Do we need to initialize hydro variables?
     */
    DensityGridInitializationFunction(DensityFunction &function, bool hydro)
        : _function(function), _hydro(hydro) {}

    /**
     * @brief Routine that sets the density for a single cell in the grid.
     *
     * @param it DensityGrid::iterator pointing to a single cell in the grid.
     */
    inline void operator()(iterator it) {

      DensityValues vals = _function(it);
      IonizationVariables &ionization_variables = it.get_ionization_variables();
      ionization_variables.set_number_density(vals.get_number_density());
      ionization_variables.set_temperature(vals.get_temperature());
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        ionization_variables.set_ionic_fraction(ion,
                                                vals.get_ionic_fraction(ion));
      }
      ionization_variables.set_cosmic_ray_factor(vals.get_cosmic_ray_factor());
      if (_hydro) {
        const CoordinateVector<> v = vals.get_velocity();
        it.get_hydro_variables().set_primitives_velocity(v);
      }
    }
  };

  void set_densities(std::pair< cellsize_t, cellsize_t > &block,
                     DensityFunction &function, int_fast32_t worksize = -1);

  /**
   * @brief Reset the mean intensity counters and update the reemission
   * probabilities for all cells.
   *
   * @param density_function DensityFunction to use to set the density in newly
   * created cells.
   */
  virtual void reset_grid(DensityFunction &density_function) {
    for (auto it = begin(); it != end(); ++it) {
      it.reset_mean_intensities();
    }
  }

  /**
   * @brief Evolve the grid in time with the given timestep.
   *
   * This method should only be implemented for moving grids.
   *
   * @param timestep Timestep with which to move the grid (in s).
   */
  virtual void evolve(double timestep) {}

  /**
   * @brief Set the velocity for the grid movement.
   *
   * This method should only be implemented for moving grids.
   *
   * @param gamma Polytropic index of the gas.
   * @param velocity_unit_in_SI Conversion factor from internal velocity unit to
   * SI units (in m s^-1).
   */
  virtual void set_grid_velocity(const double gamma,
                                 const double velocity_unit_in_SI) {}

  /**
   * @brief Get the total number of hydrogen atoms contained in the grid.
   *
   * This method is used in the unit tests to check whether the grid contains
   * the correct density field.
   *
   * @return Total number of hydrogen atoms contained in the grid.
   */
  inline double get_total_hydrogen_number() {
    double ntot = 0.;
    for (auto it = begin(); it != end(); ++it) {
      ntot +=
          it.get_ionization_variables().get_number_density() * it.get_volume();
    }
    return ntot;
  }

  /**
   * @brief Get the average temperature throughout the grid.
   *
   * This method is used in the unit tests to check whether the grid contains
   * the correct temperature field.
   *
   * @return Average temperature in the grid (in K).
   */
  inline double get_average_temperature() {
    double temperature = 0.;
    double ntot = 0.;
    for (auto it = begin(); it != end(); ++it) {
      double n =
          it.get_ionization_variables().get_number_density() * it.get_volume();
      temperature += n * it.get_ionization_variables().get_temperature();
      ntot += n;
    }
    return temperature / ntot;
  }

  /**
   * @brief Write the general grid variables to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    _box.write_restart_file(restart_writer);
    _periodicity_flags.write_restart_file(restart_writer);
    restart_writer.write(_ionization_energy_H);
    restart_writer.write(_ionization_energy_He);
    restart_writer.write(_has_hydro);
    {
      const auto size = _ionization_variables.size();
      restart_writer.write(size);
      for (std::vector< IonizationVariables >::size_type i = 0; i < size; ++i) {
        _ionization_variables[i].write_restart_file(restart_writer);
      }
    }
    {
      const auto size = _hydro_variables.size();
      restart_writer.write(size);
      for (std::vector< HydroVariables >::size_type i = 0; i < size; ++i) {
        _hydro_variables[i].write_restart_file(restart_writer);
      }
    }
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   * @param log Log to write logging info to.
   */
  inline DensityGrid(RestartReader &restart_reader, Log *log = nullptr)
      : _box(restart_reader), _periodicity_flags(restart_reader),
        _ionization_energy_H(restart_reader.read< double >()),
        _ionization_energy_He(restart_reader.read< double >()),
        _has_hydro(restart_reader.read< bool >()), _log(log) {

    {
      const std::vector< IonizationVariables >::size_type size =
          restart_reader
              .read< std::vector< IonizationVariables >::size_type >();
      _ionization_variables.resize(size);
      for (std::vector< IonizationVariables >::size_type i = 0; i < size; ++i) {
        _ionization_variables[i] = IonizationVariables(restart_reader);
      }
    }
    {
      const std::vector< HydroVariables >::size_type size =
          restart_reader.read< std::vector< HydroVariables >::size_type >();
      _hydro_variables.resize(size);
      for (std::vector< HydroVariables >::size_type i = 0; i < size; ++i) {
        _hydro_variables[i] = HydroVariables(restart_reader);
      }
    }

#ifndef USE_LOCKFREE
    _lock.resize(_ionization_variables.size());
#endif
  }
};

#endif // DENSITYGRID_HPP
