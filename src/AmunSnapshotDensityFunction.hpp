/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AmunSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density field from an Amun HDF5 snapshot
 * file: header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AMUNSNAPSHOTDENSITYFUNCTION_HPP
#define AMUNSNAPSHOTDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"

#include <string>

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads a density field from a Gadget snapshot.
 */
class AmunSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Box dimensions (in m). */
  const Box<> _box;

  /*! @brief Position shift (in fractions of the box size). */
  const CoordinateVector<> _shift;

  /*! @brief Number of cells in each dimension. */
  CoordinateVector< uint_fast32_t > _number_of_cells;

  /*! @brief Initial neutral fraction. */
  const double _initial_neutral_fraction;

  /*! @brief Number densities (in m^-3). */
  std::vector< double > _number_densities;

  /*! @brief Velocities (in km s^-1). */
  std::vector< CoordinateVector<> > _velocities;

  /*! @brief Pressures (in kg m^-1 s^-2). */
  std::vector< double > _temperatures;

public:
  AmunSnapshotDensityFunction(
      const std::string folder, const std::string prefix,
      const uint_fast32_t padding, const uint_fast32_t number_of_files,
      const Box<> box, const double number_density, const double sound_speed,
      const double temperature, const double initial_neutral_fraction,
      const CoordinateVector<> shift);

  AmunSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual ~AmunSnapshotDensityFunction();

  virtual DensityValues operator()(const Cell &cell) const;
};

#endif // AMUNSNAPSHOTDENSITYFUNCTION_HPP
