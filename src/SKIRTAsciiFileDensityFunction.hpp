/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SKIRTAsciiFileDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density grid from an ASCII text file that
 * is also compatible with SKIRT.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 */
#ifndef SKIRTASCIIFILEDENSITYFUNCTION_HPP
#define SKIRTASCIIFILEDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"

#include <string>

class Log;
class ParameterFile;
class Octree;

/**
 * @brief DensityFunction that reads a density grid from an ASCII text file that
 * is also compatible with SKIRT.
 *
 * The ASCII file should contain positions and densities. The DensityFunction
 * will return the density for the closest position.
 */
class SKIRTAsciiFileDensityFunction : public DensityFunction {
private:
  /*! @brief Positions in the file (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Number densities in the file (in m^-3). */
  std::vector< double > _number_densities;

  /*! @brief Octree used to speed up neighbour searching. */
  Octree *_octree;

public:
  SKIRTAsciiFileDensityFunction(const std::string filename,
                                const std::string xname,
                                const std::string yname,
                                const std::string zname,
                                const std::string nHname, Log *log = nullptr);
  SKIRTAsciiFileDensityFunction(ParameterFile &params, Log *log = nullptr);
  ~SKIRTAsciiFileDensityFunction();

  virtual DensityValues operator()(const Cell &cell);
};

#endif // SKIRTASCIIFILEDENSITYFUNCTION_HPP
