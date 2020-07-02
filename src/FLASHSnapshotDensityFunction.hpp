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
 * @file FLASHSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads densities from a FLASH snapshot.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FLASHSNAPSHOTDENSITYFUNCTION_HPP
#define FLASHSNAPSHOTDENSITYFUNCTION_HPP

#include "AMRGrid.hpp"
#include "DensityFunction.hpp"

#include <string>

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads densities from a FLASH snapshot.
 */
class FLASHSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief AMRGrid containing the snapshot file contents. */
  AMRGrid< DensityValues > *_grid;

  /*! @brief Flag indicating if cosmic ray heating variables should be read or
   *  not. */
  const bool _read_cosmic_ray_heating;

  /*! @brief Log to write logging info to. */
  Log *_log;

public:
  FLASHSnapshotDensityFunction(std::string filename, double temperature = -1.,
                               bool read_cosmic_ray_heating = false,
                               Log *log = nullptr);
  FLASHSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual void free();

  virtual DensityValues operator()(const Cell &cell);
};

#endif // FLASHSNAPSHOTDENSITYFUNCTION_HPP
