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
  AsciiFileDensityGridWriter(std::string prefix, std::string output_folder,
                             Log *log = nullptr);

  AsciiFileDensityGridWriter(std::string output_folder, ParameterFile &params,
                             Log *log = nullptr);

  virtual void write(DensityGrid &grid, uint_fast32_t iteration,
                     ParameterFile &params, double time = 0.,
                     const InternalHydroUnits *hydro_units = nullptr);

  virtual void write(DensitySubGridCreator< DensitySubGrid > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.);

  /**
   * @brief Write a snapshot for a split grid with hydro.
   *
   * @param grid_creator Grid.
   * @param counter Counter value to add to the snapshot file name.
   * @param params ParameterFile containing the run parameters that should be
   * written to the file.
   * @param time Simulation time (in s).
   */
  virtual void write(DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.) {}
};

#endif // ASCIIFILEDENSITYGRIDWRITER_HPP
