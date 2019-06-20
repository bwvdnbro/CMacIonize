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
 * @file AsciiFilePhotonSourceDistribution.cpp
 *
 * @brief AsciiFilePhotonSourceDistribution implementation.
 *
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#include "AsciiFilePhotonSourceDistribution.hpp"
#include "CoordinateVector.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
//#include "UnitConverter.hpp"

#include <cinttypes>
#include <fstream>

/**
 * @brief Constructor.
 *
 * @param filename Name of the ASCII text file to read.
 * @param log Log to write logging info to.
 */
AsciiFilePhotonSourceDistribution::AsciiFilePhotonSourceDistribution(
    std::string filename, Log *log)
      : _log(log) {

    std::ifstream file(filename);
    if (!file.is_open()) {
      cmac_error("Could not open file \"%s\"!", filename.c_str());
    }

    uint_fast32_t index_i = 0;
    std::string line;
    photonsourcenumber_t number_of_positions;

    while (getline(file, line)) {
      if (line[0] != '#') {
        std::stringstream linestream(line);
        linestream >> number_of_positions;
        _positions.resize(number_of_positions);
        _weights.resize(number_of_positions);
        break;
      }
    }

    while (getline(file, line)) {
      if (line[0] != '#') {
        double total_luminosity;
        std::stringstream linestream1(line);
        linestream1 >> total_luminosity;
        _total_luminosity = total_luminosity;
        break;
      }
    }

    while (getline(file, line)) {
      if (index_i == number_of_positions) {
        cmac_warning(
            "The file %s exceeds the number of number of sources provided.\n",
            filename.c_str());
        break;
      }
      if (line[0] != '#') {
        double sx, sy, sz, w;
        std::stringstream linestream2(line);
        // read source coordinates and luminosity weight
        linestream2 >> sx >> sy >> sz >> w;
        _positions[index_i][0] = sx;
        _positions[index_i][1] = sy;
        _positions[index_i][2] = sz;
        _weights[index_i] = w;
        index_i++;
      }
    }

    if (index_i < number_of_positions) {
      cmac_error("The file %s has fewer sources (%" PRIuFAST32
                 ") than needed (%" PRIuFAST32 ").\n",
                 filename.c_str(), index_i, number_of_positions);
    }

    file.close();

    if (log) {
      log->write_status("AsciiFilePhotonSourceDistribution with ",
                        number_of_positions, ".");
    }
  }

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the file (default: sinks.txt)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
AsciiFilePhotonSourceDistribution::AsciiFilePhotonSourceDistribution(
    ParameterFile &params, Log *log)
      : AsciiFilePhotonSourceDistribution(
            params.get_value< std::string >(
                "PhotonSourceDistribution:filename", "sinks.txt"),
            log) {}

/**
 * @brief Get the number of sources in the ASCII file.
 *
 * @return Number of sources.
 */
photonsourcenumber_t
AsciiFilePhotonSourceDistribution::get_number_of_sources() const {
  return _positions.size();
}

/**
 * @brief Get the position of one of the sources.
 *
 * Note that this function can alter the internal state of the
 * PhotonSourceDistribution, as for some implementations the positions are
 * decided randomly based on a RandomGenerator.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Position of the given source (in m).
 */
CoordinateVector<> AsciiFilePhotonSourceDistribution::get_position(
    photonsourcenumber_t index) {
  return _positions[index];
}

/**
 * @brief Get the weight of one of the sources.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Weight of the given source.
 */
double AsciiFilePhotonSourceDistribution::get_weight(
    photonsourcenumber_t index) const {
  return _weights[index];
}

/**
 * @brief Get the total luminosity of all sources in the ASCII file.
 *
 * @return Total luminosity (in s^-1).
 */
double AsciiFilePhotonSourceDistribution::get_total_luminosity() const {
  return _total_luminosity;
}
