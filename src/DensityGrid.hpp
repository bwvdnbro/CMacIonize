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
#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "DensityValues.hpp"
#include "Log.hpp"
#include "Photon.hpp"
#include "Timer.hpp"
#include "UnitConverter.hpp"

#include <cmath>

/**
 * @brief General interface for density grids.
 */
class DensityGrid {
protected:
  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Periodicity flags. */
  CoordinateVector< bool > _periodic;

  /*! @brief Ionization energy of hydrogen (in Hz). */
  double _ionization_energy_H;

  /*! @brief Ionization energy of helium (in Hz). */
  double _ionization_energy_He;

  /*! @brief Log to write log messages to. */
  Log *_log;

  /**
   * @brief Get the optical depth for a photon travelling the given path in the
   * given cell.
   *
   * @param ds Path length the photon traverses (in m).
   * @param cell DensityValues of the cell the photon travels in.
   * @param photon Photon.
   * @return Optical depth.
   */
  inline static double get_optical_depth(const double ds, DensityValues &cell,
                                         Photon &photon) {
    return ds * cell.get_total_density() *
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
  inline void update_integrals(const double ds, DensityValues &cell,
                               Photon &photon) const {
    if (cell.get_total_density() > 0.) {
      cell.lock();
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        IonName ion = static_cast< IonName >(i);
        cell.increase_mean_intensity(
            ion, ds * photon.get_weight() * photon.get_cross_section(ion));
      }
      cell.increase_heating_H(ds * photon.get_weight() *
                              photon.get_cross_section(ION_H_n) *
                              (photon.get_energy() - _ionization_energy_H));
      cell.increase_heating_He(ds * photon.get_weight() *
                               photon.get_cross_section(ION_He_n) *
                               (photon.get_energy() - _ionization_energy_He));
      cell.unlock();
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
   * @param cell DensityValues of the cell.
   */
  inline static void set_reemission_probabilities(double T,
                                                  DensityValues &cell) {
    double alpha_1_H = 1.58e-13 * std::pow(T * 1.e-4, -0.53);
    double alpha_A_agn = 4.18e-13 * std::pow(T * 1.e-4, -0.7);
    cell.set_pHion(alpha_1_H / alpha_A_agn);

    double alpha_1_He = 1.54e-13 * std::pow(T * 1.e-4, -0.486);
    double alpha_e_2tS = 2.1e-13 * std::pow(T * 1.e-4, -0.381);
    double alpha_e_2sS = 2.06e-14 * std::pow(T * 1.e-4, -0.451);
    double alpha_e_2sP = 4.17e-14 * std::pow(T * 1.e-4, -0.695);
    // We make sure the sum of all probabilities is 1...
    double alphaHe = alpha_1_He + alpha_e_2tS + alpha_e_2sS + alpha_e_2sP;

    cell.set_pHe_em(0, alpha_1_He / alphaHe);
    cell.set_pHe_em(1, cell.get_pHe_em(0) + alpha_e_2tS / alphaHe);
    cell.set_pHe_em(2, cell.get_pHe_em(1) + alpha_e_2sS / alphaHe);
    cell.set_pHe_em(3, cell.get_pHe_em(2) + alpha_e_2sP / alphaHe);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the grid.
   * @param periodic Periodicity flags.
   * @param log Log to write log messages to.
   */
  DensityGrid(
      Box box,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      Log *log = nullptr)
      : _box(box), _periodic(periodic), _log(log) {

    _ionization_energy_H =
        UnitConverter::to_SI< QUANTITY_FREQUENCY >(13.6, "eV");
    _ionization_energy_He =
        UnitConverter::to_SI< QUANTITY_FREQUENCY >(24.6, "eV");
  }

  virtual ~DensityGrid() {}

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
   * @brief Get the values stored in the cell with the given index.
   *
   * @param index Index of a cell.
   * @return DensityValues stored in that cell.
   */
  virtual DensityValues &get_cell_values(unsigned long index) const = 0;

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
   * @return Pointer to the values of the cell where the photon currently
   * resides, a nullptr if the photon left the box.
   */
  virtual DensityValues *interact(Photon &photon, double optical_depth) = 0;

  /**
   * @brief Index increment used in the iterator.
   *
   * More sofisticated grids (like the AMR grid) might implement their own
   * version.
   *
   * @param index Index to increase.
   */
  virtual void increase_index(unsigned long &index) { ++index; }

  /**
   * @brief Iterator to loop over the cells in the grid.
   */
  class iterator {
  private:
    /*! @brief Index of the cell the iterator is currently pointing to. */
    unsigned long _index;

    /*! @brief Reference to the DensityGrid over which we iterate. */
    DensityGrid &_grid;

  public:
    /**
     * @brief Constructor.
     *
     * @param index Index of the cell the iterator is currently pointing to.
     * @param grid DensityGrid over which we iterate.
     */
    inline iterator(unsigned long index, DensityGrid &grid)
        : _index(index), _grid(grid) {}

    /**
     * @brief Get the midpoint of the cell the iterator is pointing to.
     *
     * @return Cell midpoint (in m).
     */
    inline CoordinateVector<> get_cell_midpoint() const {
      return _grid.get_cell_midpoint(_index);
    }

    /**
     * @brief Get the DensityValues of the cell the iterator is pointing to.
     *
     * @return DensityValue the iterator is pointing to.
     */
    inline DensityValues &get_values() const {
      return _grid.get_cell_values(_index);
    }

    /**
     * @brief Get the volume of the cell the iterator is pointing to.
     *
     * @return Volume of the cell (in m^3).
     */
    inline double get_volume() const { return _grid.get_cell_volume(_index); }

    /**
     * @brief Increment operator.
     *
     * We only implemented the pre-increment version, since the post-increment
     * version creates a new object and is computationally more expensive.
     *
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator++() {
      _grid.increase_index(_index);
      return *this;
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
      return (&_grid == &it._grid && _index == it._index);
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
   * @brief Initialize the cells in the grid.
   *
   * All implementations should call this method in their constructor, after the
   * grid itself has been set up.
   *
   * @param function DensityFunction that sets the density.
   * @param ntot Total number of cells to initialize.
   */
  void initialize(DensityFunction &function, unsigned int ntot) {
    unsigned int nguess = 0.01 * ntot;
    unsigned int ninfo = 0.1 * ntot;
    unsigned int ndone = 0;
    Timer guesstimer;
    for (auto it = begin(); it != end(); ++it) {
      DensityValues &cell = it.get_values();
      DensityValues vals = function(it.get_cell_midpoint());
      cell.set_total_density(vals.get_total_density());
      cell.set_temperature(vals.get_temperature());
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        IonName ion = static_cast< IonName >(i);
        cell.set_ionic_fraction(ion, vals.get_ionic_fraction(ion));
      }
      set_reemission_probabilities(vals.get_temperature(), cell);
      ++ndone;
      if (_log) {
        if (ndone == nguess) {
          unsigned int tguess = round(99. * guesstimer.stop());
          _log->write_status("Filling grid will take approximately ", tguess,
                             " seconds.");
        }
        if (ndone % ninfo == 0) {
          unsigned int pdone = round(100. * ndone / ntot);
          _log->write_info("Did ", pdone, " percent.");
        }
      }
    }
  }

  /**
   * @brief Reset the mean intensity counters and update the reemission
   * probabilities for all cells.
   */
  virtual void reset_grid() {
    for (auto it = begin(); it != end(); ++it) {
      DensityValues &cell = it.get_values();
      set_reemission_probabilities(cell.get_temperature(), cell);
      cell.reset_mean_intensities();
    }
  }

  /**
   * @brief Get the total number of hydrogen atoms contained in the grid.
   *
   * This method is used in the unit tests to check whether the grid contains
   * the correct density field.
   *
   * @return Total number of hydrogen atoms contained in the grid.
   */
  inline double get_total_hydrogen_number() {
    double mtot = 0.;
    for (auto it = begin(); it != end(); ++it) {
      mtot += it.get_values().get_total_density() * it.get_volume();
    }
    return mtot;
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
    double mtot = 0.;
    for (auto it = begin(); it != end(); ++it) {
      double m = it.get_values().get_total_density() * it.get_volume();
      temperature += m * it.get_values().get_temperature();
      mtot += m;
    }
    return temperature / mtot;
  }
};

#endif // DENSITYGRID_HPP
