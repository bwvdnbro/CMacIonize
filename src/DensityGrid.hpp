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
#include "Timer.hpp"
#include "UnitConverter.hpp"
#include "WorkDistributor.hpp"

#ifdef USE_LOCKFREE
#include "Atomic.hpp"
#endif

#include <cmath>
#include <tuple>

/**
 * @brief General interface for density grids.
 */
class DensityGrid {
public:
  class iterator;

protected:
  /*! @brief DensityFunction defining the density field. */
  DensityFunction &_density_function;

  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Periodicity flags. */
  CoordinateVector< bool > _periodic;

  /*! @brief Ionization energy of hydrogen (in Hz). */
  double _ionization_energy_H;

  /*! @brief Ionization energy of helium (in Hz). */
  double _ionization_energy_He;

  /*! @brief Ionization calculation variables. */
  std::vector< IonizationVariables > _ionization_variables;

  /*! @brief Number densities stored in the grid (in m^-3). */
  std::vector< double > _number_density;

  /*! @brief Ionic fractions. For hydrogen and helium, these are the neutral
   *  fractions. For other elements, they are the fraction of the end product
   *  of ionization (e.g. _ionic_fraction[ION_C_p1] is the fraction of C that
   *  is in the form of C++). */
  std::vector< double > _ionic_fraction[NUMBER_OF_IONNAMES];

  /*! @brief Temperatures stored in the grid (in K). */
  std::vector< double > _temperature;

  /*! @brief Probabilities for re-emitting an ionizing photon after absorption
   *  by hydrogen. */
  std::vector< double > _hydrogen_reemission_probability;

  /*! @brief Probabilities for re-emitting an ionizing photon after absorption
   *  by helium. */
  std::vector< double > _helium_reemission_probability[4];

  /*! @brief Mean intensity integrals of ionizing radiation (without
   *  normalization factor, in m^3). */
  std::vector< double > _mean_intensity[NUMBER_OF_IONNAMES];

  /*! @brief Mean intensity of hydrogen ionizing radiation during the previous
   *  sub-step (in m^3 s^-1). */
  std::vector< double > _mean_intensity_H_old;

  /*! @brief Hydrogen neutral fraction during the previous iteration. */
  std::vector< double > _neutral_fraction_H_old;

  /*! @brief Heating due to the ionization of hydrogen (without normalization
   *  factor, in m^3 s^-1). */
  std::vector< double > _heating_H;

  /*! @brief Heating due to the ionization of helium (without normalization
   *  factor, in m^3 s^-1). */
  std::vector< double > _heating_He;

  /// hydro

  /*! @brief Flag indicating whether hydro is active or not. */
  bool _hydro;

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
   * @param cell DensityGrid::iterator pointing to the cell the photon travels
   * in.
   * @param photon Photon.
   * @return Optical depth.
   */
  inline static double get_optical_depth(double ds,
                                         const DensityGrid::iterator &cell,
                                         const Photon &photon) {
    return ds * cell.get_number_density() *
           (photon.get_cross_section(ION_H_n) *
                cell.get_ionic_fraction(ION_H_n) +
            photon.get_cross_section_He_corr() *
                cell.get_ionic_fraction(ION_He_n));
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
    if (cell.get_number_density() > 0.) {
      // we tried speeding things up by using lock-free addition, but it turns
      // out that the overhead caused by doing this is larger than the overhead
      // of using a single lock
      // this is mainly because we have to do a large number of additions
      // to minimize collisions (two threads trying to access the same cell at
      // the same time), we first calculate all terms that need to be added, and
      // then lock the cell and do all additions as fast as possible
      double dmean_intensity[NUMBER_OF_IONNAMES];
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        IonName ion = static_cast< IonName >(i);
        dmean_intensity[i] =
            ds * photon.get_weight() * photon.get_cross_section(ion);
      }
      double dheating_H = ds * photon.get_weight() *
                          photon.get_cross_section(ION_H_n) *
                          (photon.get_energy() - _ionization_energy_H);
      double dheating_He = ds * photon.get_weight() *
                           photon.get_cross_section(ION_He_n) *
                           (photon.get_energy() - _ionization_energy_He);
#ifndef USE_LOCKFREE
      cell.lock();
#endif
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        IonName ion = static_cast< IonName >(i);
        cell.increase_mean_intensity(ion, dmean_intensity[i]);
      }
      cell.increase_heating_H(dheating_H);
      cell.increase_heating_He(dheating_He);
#ifndef USE_LOCKFREE
      cell.unlock();
#endif
    }
  }

protected:
  /**
   * @brief Set the re-emission probabilities for the given cell for the given
   * temperature.
   *
   * These quantities are all dimensionless.
   *
   * @param T Temperature (in K).
   * @param it DensityGrid::iterator pointing to a cell.
   */
  inline static void set_reemission_probabilities(double T,
                                                  DensityGrid::iterator &it) {
    // reemission probabilities
    double alpha_1_H = 1.58e-13 * std::pow(T * 1.e-4, -0.53);
    double alpha_A_agn = 4.18e-13 * std::pow(T * 1.e-4, -0.7);
    it.set_hydrogen_reemission_probability(alpha_1_H / alpha_A_agn);
    double alpha_1_He = 1.54e-13 * std::pow(T * 1.e-4, -0.486);
    double alpha_e_2tS = 2.1e-13 * std::pow(T * 1.e-4, -0.381);
    double alpha_e_2sS = 2.06e-14 * std::pow(T * 1.e-4, -0.451);
    double alpha_e_2sP = 4.17e-14 * std::pow(T * 1.e-4, -0.695);
    // We make sure the sum of all probabilities is 1...
    double alphaHe = alpha_1_He + alpha_e_2tS + alpha_e_2sS + alpha_e_2sP;

    it.set_helium_reemission_probability(0, alpha_1_He / alphaHe);
    it.set_helium_reemission_probability(
        1, it.get_helium_reemission_probability(0) + alpha_e_2tS / alphaHe);
    it.set_helium_reemission_probability(
        2, it.get_helium_reemission_probability(1) + alpha_e_2sS / alphaHe);
    it.set_helium_reemission_probability(
        3, it.get_helium_reemission_probability(2) + alpha_e_2sP / alphaHe);
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
   * @param density_function DensityFunction that defines the density field.
   * @param box Box containing the grid.
   * @param periodic Periodicity flags.
   * @param hydro Hydro flag.
   * @param log Log to write log messages to.
   */
  DensityGrid(
      DensityFunction &density_function, Box box,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      bool hydro = false, Log *log = nullptr)
      : _density_function(density_function), _box(box), _periodic(periodic),
        _hydro(hydro), _log(log) {

    _ionization_energy_H =
        UnitConverter::to_SI< QUANTITY_FREQUENCY >(13.6, "eV");
    _ionization_energy_He =
        UnitConverter::to_SI< QUANTITY_FREQUENCY >(24.6, "eV");
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~DensityGrid() {}

  /**
   * @brief Allocate memory for the given number of cells.
   *
   * @param numcell Number of cells that will be stored in the grid.
   */
  inline void allocate_memory(unsigned long numcell) {
    if (_log) {
      _log->write_status(
          "Allocating memory for ", numcell, " cells (",
          Utilities::human_readable_bytes(numcell * sizeof(DensityValues)),
          ")...");
    }
    // we allocate memory for the cells, so that --dry-run can already check the
    // available memory
    _ionization_variables.resize(numcell);
    _number_density.resize(numcell);
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _ionic_fraction[i].resize(numcell);
    }
    _temperature.resize(numcell);
    _hydrogen_reemission_probability.resize(numcell);
    for (int i = 0; i < 4; ++i) {
      _helium_reemission_probability[i].resize(numcell);
    }
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _mean_intensity[i].resize(numcell);
    }
    _mean_intensity_H_old.resize(numcell);
    _neutral_fraction_H_old.resize(numcell);
    _heating_H.resize(numcell);
    _heating_He.resize(numcell);
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
   */
  virtual void initialize(std::pair< unsigned long, unsigned long > &block) {
    if (_log) {
      _log->write_status("Initializing DensityFunction...");
    }
    _density_function.initialize();
    if (_log) {
      _log->write_status("Done.");
    }
    if (_hydro) {
      if (_log) {
        _log->write_status("Initializing hydro arrays...");
      }
      const unsigned int numcell = get_number_of_cells();
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
  virtual unsigned int get_number_of_cells() const = 0;

  /**
   * @brief Get the Box containing the grid.
   *
   * @return Box containing the grid (in m).
   */
  inline Box get_box() const { return _box; }

  /**
   * @brief Get the number of periodic boundaries of this grid.
   *
   * @return Number of periodic boundaries of this grid (between 0 and 3).
   */
  inline unsigned int get_number_of_periodic_boundaries() const {
    unsigned int numperiodic = 0;
    for (unsigned int i = 0; i < 3; ++i) {
      numperiodic += _periodic[i];
    }
    return numperiodic;
  }

  /**
   * @brief Get the index of the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return Index of the cell containing that position.
   */
  virtual unsigned long get_cell_index(CoordinateVector<> position) const = 0;

  /**
   * @brief Get the midpoint of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Midpoint of that cell (in m).
   */
  virtual CoordinateVector<> get_cell_midpoint(unsigned long index) const = 0;

  /**
   * @brief Get a handle to the number density vector that can be used in MPI
   * communications.
   *
   * @return Reference to the internal std::vector.
   */
  inline std::vector< double > &get_number_density_handle() {
    return _number_density;
  }

  /**
   * @brief Get a handle to the temperature vector that can be used in MPI
   * communications.
   *
   * @return Reference to the internal std::vector.
   */
  inline std::vector< double > &get_temperature_handle() {
    return _temperature;
  }

  /**
   * @brief Get a handle to the ionic fraction vector for the given ion that can
   * be used in MPI communications.
   *
   * @param ion IonName.
   * @return Reference to the internal std::vector.
   */
  inline std::vector< double > &get_ionic_fraction_handle(IonName ion) {
    return _ionic_fraction[ion];
  }

  /**
   * @brief Get a handle to the mean intensity vector for the given ion that can
   * be used in MPI communications.
   *
   * @param ion IonName.
   * @return Reference to the internal std::vector.
   */
  inline std::vector< double > &get_mean_intensity_handle(IonName ion) {
    return _mean_intensity[ion];
  }

  /**
   * @brief Get a handle to the hydrogen ionization heating vector that can be
   * used in MPI communications.
   *
   * @return Reference to the internal std::vector.
   */
  inline std::vector< double > &get_heating_H_handle() { return _heating_H; }

  /**
   * @brief Get a handle to the helium ionization heating vector that can be
   * used in MPI communications.
   *
   * @return Reference to the internal std::vector.
   */
  inline std::vector< double > &get_heating_He_handle() { return _heating_He; }

  /**
   * @brief Check if hydro is active.
   *
   * @return True if hydro is active.
   */
  inline bool has_hydro() const { return _hydro; }

  /**
   * @brief Get the neighbours of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return std::vector containing iterators to the neighbours, together with
   * the midpoint, surface normal and surface area of the boundary face between
   * the cell and this neighbour.
   */
  virtual std::vector<
      std::tuple< iterator, CoordinateVector<>, CoordinateVector<>, double > >
  get_neighbours(unsigned long index) = 0;

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
  virtual double get_cell_volume(unsigned long index) const = 0;

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
  virtual void increase_index(unsigned long &index,
                              unsigned long increment = 1) {
    index += increment;
  }

  /**
   * @brief Iterator to loop over the cells in the grid.
   */
  class iterator {
  private:
    /*! @brief Index of the cell the iterator is currently pointing to. */
    unsigned long _index;

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
    inline iterator(unsigned long index, DensityGrid &grid)
        : _index(index), _grid(&grid) {}

    /**
     * @brief Get the midpoint of the cell the iterator is pointing to.
     *
     * @return Cell midpoint (in m).
     */
    inline CoordinateVector<> get_cell_midpoint() const {
      return _grid->get_cell_midpoint(_index);
    }

    /**
     * @brief Get the number density of hydrogen for the cell the iterator is
     * currently pointing to.
     *
     * @return Number density of hydrogen (in m^-3).
     */
    inline double get_number_density() const {
      return _grid->_number_density[_index];
    }

    /**
     * @brief set the number density of hydrogen for the cell the iterator is
     * currently pointing to.
     *
     * @param number_density Number density of hydrogen (in m^-3).
     */
    inline void set_number_density(double number_density) {
      _grid->_number_density[_index] = number_density;
    }

    /**
     * @brief Get the temperature of the cell the iterator is currently pointing
     * to.
     *
     * @return Temperature (in K).
     */
    inline double get_temperature() const {
      return _grid->_temperature[_index];
    }

    /**
     * @brief Set the temperature of the cell the iterator is currently pointing
     * to.
     *
     * @param temperature Temperature (in K).
     */
    inline void set_temperature(double temperature) {
      _grid->_temperature[_index] = temperature;
    }

    /**
     * @brief Get the ionic fraction for the given IonName for the cell the
     * iterator is currently pointing to.
     *
     * @param ion IonName.
     * @return Ionic fraction for that ion.
     */
    inline double get_ionic_fraction(IonName ion) const {
      return _grid->_ionic_fraction[ion][_index];
    }

    /**
     * @brief Set the ionic fraction for the given IonName for the cell the
     * iterator is currently pointing to.
     *
     * @param ion IonName.
     * @param ionic_fraction Ionic fraction for that ion.
     */
    inline void set_ionic_fraction(IonName ion, double ionic_fraction) {
      _grid->_ionic_fraction[ion][_index] = ionic_fraction;
    }

    /**
     * @brief Get the neutral fraction of hydrogen during the previous iteration
     * for the cell the iterator is currently pointing to.
     *
     * @return Neutral fraction of hydrogen during the previous iteration.
     */
    inline double get_neutral_fraction_H_old() const {
      return _grid->_neutral_fraction_H_old[_index];
    }

    /**
     * @brief Set the neutral fraction of hydrogen during the previous iteration
     * for the cell the iterator is currently pointing to.
     *
     * @param neutral_fraction_H_old Neutral fraction of hydrogen during the
     * previous iteration.
     */
    inline void set_neutral_fraction_H_old(double neutral_fraction_H_old) {
      _grid->_neutral_fraction_H_old[_index] = neutral_fraction_H_old;
    }

    /**
     * @brief Get the mean intensity of ionizing radiation for the given IonName
     * for the cell the iterator is currently pointing to.
     *
     * @param ion IonName.
     * @return Mean intensity of ionizing radiation for that ion (without
     * normalization factor, in m^3).
     */
    inline double get_mean_intensity(IonName ion) const {
      return _grid->_mean_intensity[ion][_index];
    }

    /**
     * @brief Increase the mean ionizing intensity for the given ion for the
     * cell the iterator is currently pointing to.
     *
     * @param ion IonName.
     * @param mean_intensity_increment Mean ionizing intensity increment
     * (without normalization factor, in m^3).
     */
    inline void increase_mean_intensity(IonName ion,
                                        double mean_intensity_increment) {
#ifdef USE_LOCKFREE
      Atomic::add(_grid->_mean_intensity[ion][_index],
                  mean_intensity_increment);
#else
      _grid->_mean_intensity[ion][_index] += mean_intensity_increment;
#endif
    }

    /**
     * @brief Get the heating by ionization of hydrogen for the cell the
     * iterator is currently pointing to.
     *
     * @return Heating by ionization of hydrogen (in m^3 s^-1).
     */
    inline double get_heating_H() const { return _grid->_heating_H[_index]; }

    /**
     * @brief Increase the ionizing heating for hydrogen for the cell the
     * iterator is currently pointing to.
     *
     * @param heating_H_increment Heating by ionization of hydrogen (in m^3
     * s^-1).
     */
    inline void increase_heating_H(double heating_H_increment) {
#ifdef USE_LOCKFREE
      Atomic::add(_grid->_heating_H[_index], heating_H_increment);
#else
      _grid->_heating_H[_index] += heating_H_increment;
#endif
    }

    /**
     * @brief Get the heating by ionization of helium for the cell the iterator
     * is currently pointing to.
     *
     * @return Heating by ionization of helium (in m^3 s^-1).
     */
    inline double get_heating_He() const { return _grid->_heating_He[_index]; }

    /**
     * @brief Increase the ionizing heating for helium for the cell the iterator
     * is currently pointing to.
     *
     * @param heating_He_increment Heating by ionization of helium (in m^3
     * s^-1).
     */
    inline void increase_heating_He(double heating_He_increment) {
#ifdef USE_LOCKFREE
      Atomic::add(_grid->_heating_He[_index], heating_He_increment);
#else
      _grid->_heating_He[_index] += heating_He_increment;
#endif
    }

    /**
     * @brief Get the mean intensity of hydrogen ionizing radiation during the
     * previous iteration for the cell the iterator is currently pointing to.
     *
     * @return Mean intensity of hydrogen ionizing radiation during the previous
     * iteration (without normalization factor, in m^3).
     */
    inline double get_mean_intensity_H_old() const {
      return _grid->_mean_intensity_H_old[_index];
    }

    /**
     * @brief Set the mean intensity of hydrogen ionizing radiation during the
     * previous iteration for the cell the iterator is currently pointing to.
     *
     * @param mean_intensity_H_old Mean intensity of hydrogen ionizing radiation
     * during the previous iteration (without normalization factor, in m^3).
     */
    inline void set_mean_intensity_H_old(double mean_intensity_H_old) {
      _grid->_mean_intensity_H_old[_index] = mean_intensity_H_old;
    }

    /**
     * @brief Get the probability for hydrogen of reemitting an ionizing photon
     * for the cell the iterator is currently pointing to.
     *
     * @return Hydrogen ionizing reemission probability.
     */
    inline double get_hydrogen_reemission_probability() const {
      return _grid->_hydrogen_reemission_probability[_index];
    }

    /**
     * @brief Set the probability for hydrogen of reemitting an ionizing photon
     * for the cell the iterator is currently pointing to.
     *
     * @param hydrogen_reemission_probability Hydrogen ionizing reemission
     * probability.
     */
    inline void set_hydrogen_reemission_probability(
        double hydrogen_reemission_probability) {
      _grid->_hydrogen_reemission_probability[_index] =
          hydrogen_reemission_probability;
    }

    /**
     * @brief Get the helium reemission probability in the given channel for the
     * cell the iterator is currently pointing to.
     *
     * @param channel Channel in which the radiation is emitted.
     * @return Helium reemission probability in that channel.
     */
    inline double
    get_helium_reemission_probability(unsigned char channel) const {
      return _grid->_helium_reemission_probability[channel][_index];
    }

    /**
     * @brief Set the helium reemission probability in the given channel for the
     * cell the iterator is currently pointing to.
     *
     * @param channel Channel in which the radiation is emitted.
     * @param helium_reemission_probability Helium reemission probability in
     * that channel.
     */
    inline void
    set_helium_reemission_probability(unsigned char channel,
                                      double helium_reemission_probability) {
      _grid->_helium_reemission_probability[channel][_index] =
          helium_reemission_probability;
    }

    /**
     * @brief Reset the mean intensity counters for the cell the iterator is
     * currently pointing to.
     */
    inline void reset_mean_intensities() {
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        _grid->_mean_intensity[i][_index] = 0.;
      }
      _grid->_mean_intensity_H_old[_index] = 0.;
      _grid->_heating_H[_index] = 0.;
      _grid->_heating_He[_index] = 0.;
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
        std::tuple< iterator, CoordinateVector<>, CoordinateVector<>, double > >
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
    inline double get_volume() const { return _grid->get_cell_volume(_index); }

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
    inline iterator &operator+=(unsigned long increment) {
      _grid->increase_index(_index, increment);
      return *this;
    }

    /**
     * @brief Free addition operator.
     *
     * @param increment Increment to add to the iterator.
     * @return Incremented iterator.
     */
    inline iterator operator+(unsigned long increment) {
      iterator it(*this);
      it += increment;
      return it;
    }

    /**
     * @brief Get the index of the cell the iterator is currently pointing to.
     *
     * @return Index of the current cell.
     */
    inline unsigned long get_index() const { return _index; }

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
  inline std::pair< iterator, iterator > get_chunk(unsigned long begin,
                                                   unsigned long end) {
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
      DensityValues vals = _function(it.get_cell_midpoint());
      it.set_number_density(vals.get_number_density());
      it.set_temperature(vals.get_temperature());
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        IonName ion = static_cast< IonName >(i);
        it.set_ionic_fraction(ion, vals.get_ionic_fraction(ion));
      }
      if (_hydro) {
        const CoordinateVector<> v = vals.get_velocity();
        it.get_hydro_variables().set_primitives_velocity(v);
      }
      set_reemission_probabilities(it.get_temperature(), it);
    }
  };

  void initialize(std::pair< unsigned long, unsigned long > &block,
                  DensityFunction &function, int worksize = -1);

  /**
   * @brief Reset the mean intensity counters and update the reemission
   * probabilities for all cells.
   */
  virtual void reset_grid() {
    for (auto it = begin(); it != end(); ++it) {
      set_reemission_probabilities(it.get_temperature(), it);
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
   */
  virtual void set_grid_velocity() {}

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
      ntot += it.get_number_density() * it.get_volume();
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
      double n = it.get_number_density() * it.get_volume();
      temperature += n * it.get_temperature();
      ntot += n;
    }
    return temperature / ntot;
  }
};

#endif // DENSITYGRID_HPP
