/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Maya Petkova (map32@st-andrews.ac.uk)
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
 * @file SPHVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution implementation with generating sites
 * imported from a file.
 *
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#ifndef SPHVORONOIGENERATORDISTRIBUTION_HPP
#define SPHVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "VoronoiGeneratorDistribution.hpp"

#include <cinttypes>
#include <fstream>

/**
 * @brief VoronoiGeneratorDistribution implementation with generating sites
 * imported from a file.
 */
class SPHVoronoiGeneratorDistribution : public VoronoiGeneratorDistribution {
private:
  /*! @brief Number of generator positions to generate. */
  const unsigned int _number_of_positions;

  /*! @brief Number of positions already generated. */
  unsigned int _current_number;

  /*! @brief Box containing the generators (in m). */
  const Box<> _box;

  /*! @brief Name of ASCII file used to import Voronoi generating sites from. */
  std::string _filename;

  /*! @brief Generator positions. */
  std::vector< CoordinateVector<> > _generator_positions;

public:
  /**
   * @brief Constructor.
   *
   * @param simulation_box Simulation box (in m).
   * @param number_of_positions Number of SPH generator positions to
   * generate.
   * @param filename Name of the ASCII text file to read.
   * @param log Log to write logging info to.
   */
  SPHVoronoiGeneratorDistribution(const Box<> &simulation_box,
                                  uint_fast32_t number_of_positions,
                                  std::string filename, Log *log = nullptr)
      : _number_of_positions(number_of_positions), _current_number(0),
        _box(simulation_box), _filename(filename) {

    _generator_positions.resize(number_of_positions);

    std::ifstream file(filename);
    if (!file.is_open()) {
      cmac_error("Could not open file \"%s\"!", filename.c_str());
    }

    uint_fast32_t index_i = 0;
    std::string line;
    while (getline(file, line)) {
      if (index_i == number_of_positions) {
        cmac_warning(
            "The file %s exceeds the number of generator positions provided.\n",
            filename.c_str());
        break;
      }
      if (line[0] != '#') {
        double cx, cy, cz;
        std::stringstream linestream(line);
        // read Voronoi cell generating cite coordinates
        linestream >> cx >> cy >> cz;
        _generator_positions[index_i][0] = cx;
        _generator_positions[index_i][1] = cy;
        _generator_positions[index_i][2] = cz;
        index_i++;
      }
    }

    if (index_i < number_of_positions) {
      cmac_error("The file %s has fewer generator positions (%" PRIuFAST32
                 ") than needed (%" PRIuFAST32 ").\n",
                 filename.c_str(), index_i, number_of_positions);
    }

    file.close();

    if (log) {
      log->write_status("SPHVoronoiGeneratorDistribution with ",
                        _number_of_positions, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - number of positions: Number of positions in the file (default: 1000)
   *  - filename: Name of the file (default: SPH.txt)
   *
   * @param simulation_box Simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SPHVoronoiGeneratorDistribution(const Box<> &simulation_box,
                                  ParameterFile &params, Log *log = nullptr)
      : SPHVoronoiGeneratorDistribution(
            simulation_box,
            params.get_value< uint_fast32_t >(
                "DensityGrid:VoronoiGeneratorDistribution:number of positions",
                1000),
            params.get_value< std::string >(
                "DensityGrid:VoronoiGeneratorDistribution:filename", "SPH.txt"),
            log) {}

  /**
   * @brief Get the number of positions that this distribution generates.
   *
   * @return Number of generated positions.
   */
  virtual unsigned int get_number_of_positions() const {
    return _number_of_positions;
  }

  /**
   * @brief Get SPH generator position.
   *
   * @return SPH generator position (in m).
   */
  virtual CoordinateVector<> get_position() {
    cmac_assert(_current_number < _generator_positions.size());
    CoordinateVector<> pos = _generator_positions[_current_number];
    ++_current_number;
    return pos;
  }
};

#endif // SPHVORONOIGENERATORDISTRIBUTION_HPP
