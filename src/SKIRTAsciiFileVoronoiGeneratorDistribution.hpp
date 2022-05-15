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
 * @file SKIRTAsciiFileVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution that reads generator positions from an
 * ASCII text file that is also compatible with SKIRT.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 */
#ifndef SKIRTASCIIFILEVORONOIGENERATORDISTRIBUTION_HPP
#define SKIRTASCIIFILEVORONOIGENERATORDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "SKIRTAsciiFile.hpp"
#include "VoronoiGeneratorDistribution.hpp"

#include <string>

/**
 * @brief VoronoiGeneratorDistribution that reads generator positions from an
 * ASCII text file that is also compatible with SKIRT.
 *
 * The ASCII file should at least contain positions.
 */
class SKIRTAsciiFileVoronoiGeneratorDistribution
    : public VoronoiGeneratorDistribution {
private:
  /*! @brief Positions in the file (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Next generator index. */
  size_t _next_index;

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the ASCII text file to read.
   * @param xname Name of the x positions column.
   * @param yname Name of the y positions column.
   * @param zname Name of the z positions column.
   * @param log Log to write logging info to.
   */
  SKIRTAsciiFileVoronoiGeneratorDistribution(const std::string filename,
                                             const std::string xname,
                                             const std::string yname,
                                             const std::string zname, Log *log)
      : _next_index(0) {

    if (log) {
      log->write_info("Initialising SKIRTAsciiFileVoronoiGeneratorDistribution "
                      "from file \"",
                      filename, "\", reading fields \"", xname, "\", \"", yname,
                      "\" and \"", zname, "\".");
    }

    SKIRTAsciiFile file(filename);

    if (!file.has_column(xname) || !file.is_quantity(xname, QUANTITY_LENGTH)) {
      cmac_error("Could not initiate x positions from column \"%s\"!",
                 xname.c_str());
    }
    if (!file.has_column(yname) || !file.is_quantity(yname, QUANTITY_LENGTH)) {
      cmac_error("Could not initiate y positions from column \"%s\"!",
                 yname.c_str());
    }
    if (!file.has_column(zname) || !file.is_quantity(zname, QUANTITY_LENGTH)) {
      cmac_error("Could not initiate z positions from column \"%s\"!",
                 zname.c_str());
    }

    const std::vector< double > &x = file.get_column(xname);
    const std::vector< double > &y = file.get_column(yname);
    const std::vector< double > &z = file.get_column(zname);

    _positions.resize(file.number_of_rows());
    for (size_t i = 0; i < file.number_of_rows(); ++i) {
      _positions[i][0] = x[i];
      _positions[i][1] = y[i];
      _positions[i][2] = z[i];
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - filename: Name of the ASCII file (required)
   *  - xname: Name of the x positions column (default: x-coordinate)
   *  - yname: Name of the y positions column (default: y-coordinate)
   *  - zname: Name of the z positions column (default: z-coordinate)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SKIRTAsciiFileVoronoiGeneratorDistribution(ParameterFile &params,
                                             Log *log = nullptr)
      : SKIRTAsciiFileVoronoiGeneratorDistribution(
            params.get_filename(
                "DensityGrid:VoronoiGeneratorDistribution:filename"),
            params.get_value< std::string >(
                "DensityGrid:VoronoiGeneratorDistribution:xname",
                "x-coordinate"),
            params.get_value< std::string >(
                "DensityGrid:VoronoiGeneratorDistribution:yname",
                "y-coordinate"),
            params.get_value< std::string >(
                "DensityGrid:VoronoiGeneratorDistribution:zname",
                "z-coordinate"),
            log) {}

  /**
   * @brief Get the number of positions that this distribution generates.
   *
   * @return Number of positions in the file.
   */
  virtual generatornumber_t get_number_of_positions() const {
    return _positions.size();
  }

  /**
   * @brief Get the next generator position.
   *
   * @return Next generator position (in m).
   */
  virtual CoordinateVector<> get_position() {

    const CoordinateVector<> pos = _positions[_next_index];

    ++_next_index;
    cmac_assert(_next_index <= get_number_of_positions());

    return pos;
  }
};

#endif // SKIRTASCIIFILEVORONOIGENERATORDISTRIBUTION_HPP
