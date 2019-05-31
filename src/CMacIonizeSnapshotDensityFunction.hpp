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
 * @file CMacIonizeSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density grid from a CMacIonize snapshot.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CMACIONIZESNAPSHOTDENSITYFUNCTION_HPP
#define CMACIONIZESNAPSHOTDENSITYFUNCTION_HPP

#include "AMRGrid.hpp"
#include "Box.hpp"
#include "DensityFunction.hpp"
#include "PointLocations.hpp"

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads a density grid from a CMacIonize snapshot.
 */
class CMacIonizeSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Name of the snapshot file. */
  const std::string _filename;

  /*! @brief Use the mass density instead of the number density? */
  const bool _use_density;

  /*! @brief Use the pressure instead of the temperature? */
  const bool _use_pressure;

  /*! @brief Initial value for the neutral fractions (if not present). */
  const double _initial_neutral_fraction;

  /*! @brief Log to write logging info to. */
  Log *_log;

  /*! @brief Box containing the grid. */
  Box<> _box;

  /*! @brief Number of cells in each dimension. */
  CoordinateVector< uint_fast32_t > _ncell;

  /*! @brief Cartesian density grid (if applicable). */
  DensityValues ***_cartesian_grid;

  /*! @brief AMR density grid (if applicable). */
  AMRGrid< DensityValues > *_amr_grid;

  /*! @brief Locations of the Voronoi grid generators (if applicable). */
  std::vector< CoordinateVector<> > _voronoi_generators;

  /*! @brief PointLocations object used to query the Voronoi grid generator
   *  locations (if applicable). */
  PointLocations *_voronoi_pointlocations;

  /*! @brief DensityValues stored in the Voronoi grid (if applicable). */
  std::vector< DensityValues > _voronoi_densityvalues;

  /**
   * @brief Get the largest odd factor of the given number.
   *
   * This is the number you get by iteratively dividing the number by two, until
   * the result is no longer even.
   *
   * @param number Number to decompose.
   * @return Largest odd factor of the number.
   */
  inline static uint_fast32_t get_largest_odd_factor(uint_fast32_t number) {
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
  inline static uint_fast32_t get_power_of_two(uint_fast32_t number) {
    return number / get_largest_odd_factor(number);
  }

public:
  CMacIonizeSnapshotDensityFunction(std::string filename,
                                    const bool use_density,
                                    const bool use_pressure,
                                    const double initial_neutral_fraction,
                                    Log *log = nullptr);

  CMacIonizeSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual ~CMacIonizeSnapshotDensityFunction();

  virtual void initialize();

  virtual void free();

  virtual DensityValues operator()(const Cell &cell) const;
};

#endif // CMACIONIZESNAPSHOTDENSITYFUNCTION_HPP
