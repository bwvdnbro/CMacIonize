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
 * @brief HDF5-file writer for the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDWRITER_HPP
#define DENSITYGRIDWRITER_HPP

#include <string>

class DensityGrid;
class Log;
class ParameterFile;

/**
 * @brief HDF5-file writer for the DensityGrid.
 */
class DensityGridWriter {
private:
  /*! @brief Prefix of the name for the file to write. */
  std::string _prefix;

  /*! @brief DensityGrid containing the data to write. */
  DensityGrid &_grid;

  /*! @brief Number of digits used for the counter in the filenames. */
  unsigned char _padding;

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  DensityGridWriter(std::string prefix, DensityGrid &grid, Log *log = NULL,
                    unsigned char padding = 3);
  DensityGridWriter(ParameterFile &params, DensityGrid &grid, Log *log = NULL);

  void write(unsigned int iteration);
};

#endif // DENSITYGRIDWRITER_HPP
