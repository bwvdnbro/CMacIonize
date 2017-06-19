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
 * @file AsciiFileDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density grid from an ASCII text file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ASCIIFILEDENSITYFUNCTION_HPP
#define ASCIIFILEDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"

#include <string>

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads a density grid from an ASCII text file.
 */
class AsciiFileDensityFunction : public DensityFunction {
private:
  /*! @brief Density grid (in m^-3). */
  double ***_grid;

  /*! @brief Dimensions of the grid. */
  CoordinateVector< int > _ncell;

  /*! @brief Box containing the grid. */
  Box<> _box;

  /*! @brief Initial temperature of the ISM (in K). */
  double _temperature;

  /*! @brief Log to write logging info to. */
  Log *_log;

public:
  AsciiFileDensityFunction(std::string filename, CoordinateVector< int > ncell,
                           Box<> box, double temperature,
                           double length_unit_in_SI = 1.,
                           double density_unit_in_SI = 1., Log *log = nullptr);
  AsciiFileDensityFunction(ParameterFile &params, Log *log = nullptr);
  ~AsciiFileDensityFunction();

  virtual DensityValues operator()(CoordinateVector<> position) const;

  double get_total_hydrogen_number() const;
};

#endif // ASCIIFILEDENSITYFUNCTION_HPP
