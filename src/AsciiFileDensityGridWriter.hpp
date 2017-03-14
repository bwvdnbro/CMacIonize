/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AsciiFileDensityGridWriter.hpp
 *
 * @brief DensityGridWriter instance that writes an ASCII file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ASCIIFILEDENSITYGRIDWRITER_HPP
#define ASCIIFILEDENSITYGRIDWRITER_HPP

#include "DensityGridWriter.hpp"

/**
 * @brief DensityGridWriter instance that writes an ASCII file.
 */
class AsciiFileDensityGridWriter : public DensityGridWriter {
private:
  /*! @brief Prefix of snapshot file names. */
  std::string _prefix;

public:
  AsciiFileDensityGridWriter(std::string prefix, DensityGrid &grid,
                             std::string output_folder, Log *log = nullptr);

  AsciiFileDensityGridWriter(ParameterFile &params, DensityGrid &grid,
                             Log *log = nullptr);

  virtual void write(unsigned int iteration, ParameterFile &params);
};

#endif // ASCIIFILEDENSITYGRIDWRITER_HPP
