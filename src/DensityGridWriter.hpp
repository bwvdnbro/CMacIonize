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
 * @file DensityGridWriter.hpp
 *
 * @brief Snapshot file writer for the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDWRITER_HPP
#define DENSITYGRIDWRITER_HPP

#include <cstdlib>

class DensityGrid;
class Log;
class ParameterFile;

/**
 * @brief Snapshot file writer for the DensityGrid.
 */
class DensityGridWriter {
protected:
  /*! @brief DensityGrid containing the data to write. */
  DensityGrid &_grid;

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  /**
   * @brief Constructor.
   *
   * @param grid Grid to write out.
   * @param log Log to write logging information to.
   */
  DensityGridWriter(DensityGrid &grid, Log *log = nullptr)
      : _grid(grid), _log(log) {}

  virtual ~DensityGridWriter() {}

  /**
   * @brief Write a snapshot.
   *
   * @param iteration Iteration number to use in the snapshot file name(s).
   * @param params ParameterFile containing the run parameters that should be
   * written to the file.
   */
  virtual void write(unsigned int iteration, ParameterFile &params) = 0;
};

#endif // DENSITYGRIDWRITER_HPP
