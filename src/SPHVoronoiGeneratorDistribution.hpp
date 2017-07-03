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
 * @file SPHVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution implementation with generating sites 
 * imported from a file.
 *
 * @author Maya A. Petkova
 */
#ifndef SPHVORONOIGENERATORDISTRIBUTION_HPP
#define SPHVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "VoronoiGeneratorDistribution.hpp"
#include "CoordinateVector.hpp"

#include <fstream>

/**
 * @brief VoronoiGeneratorDistribution implementation with generating sites 
 * imported from a file.
 */
class SPHVoronoiGeneratorDistribution
    : public VoronoiGeneratorDistribution {
private:
  /*! @brief Number of random generator positions to generate. */
  const unsigned int _number_of_positions;

  /*! @brief Number of positions already generated. */
  unsigned int _current_number;

  /*! @brief Box containing the generators (in m). */
  const Box<> _box;

  /*! @brief Name of ASCII file used to import Voronoi generating sites from. */
  std::string _filename;

  /*! @brief RandomGenerator used to generate the positions. */
  std::vector< CoordinateVector<> > _generator_positions;

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the generators (in m).
   * @param number_of_positions Number of SPH generator positions to
   * generate.
   * @param filename Name of the ASCII text file to read.
   * @param log Log to write logging info to.
   */
  SPHVoronoiGeneratorDistribution(Box<> box,
				  unsigned int number_of_positions,
				  std::string filename,
				  Log *log = nullptr)
      : _number_of_positions(number_of_positions), _current_number(0),
        _box(box), _filename(filename) {

    _generator_positions.resize(number_of_positions);

  std::ifstream file(filename);
  if (!file.is_open()) {
    cmac_error("Could not open file \"%s\"!", filename.c_str());
  }

  unsigned int index_i=0;
  std::string line;
  while (getline(file, line)) {
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
   
  file.close();


    if (log) {
      log->write_status(
          "SPHVoronoiGeneratorDistribution with ",
          _number_of_positions, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SPHVoronoiGeneratorDistribution(ParameterFile &params,
                                            Log *log = nullptr)
      : SPHVoronoiGeneratorDistribution(
					Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_anchor", "[0. m, 0. m, 0. m]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_sides", "[1. m, 1. m, 1. m]")),
            params.get_value< unsigned int >("densitygrid:voronoi_generator_"
                                             "distribution:number_of_positions",
                                             1000),
            params.get_value< std::string >("densitygrid:voronoi_generator_"
                                             "distribution:file_name",
                                             "SPH.txt"),
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
