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
 * @file UniformRandomVoronoiGeneratorDistribution.hpp
 *
 * @brief Uniform random VoronoiGeneratorDistribution implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UNIFORMRANDOMVORONOIGENERATORDISTRIBUTION_HPP
#define UNIFORMRANDOMVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "VoronoiGeneratorDistribution.hpp"

/**
 * @brief Uniform random VoronoiGeneratorDistribution implementation.
 */
class UniformRandomVoronoiGeneratorDistribution
    : public VoronoiGeneratorDistribution {
private:
  /*! @brief Number of random generator positions to generate. */
  const unsigned int _number_of_positions;

  /*! @brief Number of positions already generated. */
  unsigned int _current_number;

  /*! @brief Box containing the generators (in m). */
  const Box _box;

  /*! @brief RandomGenerator used to generate the positions. */
  RandomGenerator _random_generator;

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing hte generators (in m).
   * @param number_of_positions Number of random generator positions to
   * generate.
   * @param random_seed Seed for the random number generator.
   * @param log Log to write logging info to.
   */
  UniformRandomVoronoiGeneratorDistribution(Box box,
                                            unsigned int number_of_positions,
                                            int random_seed, Log *log = nullptr)
      : _number_of_positions(number_of_positions), _current_number(0),
        _box(box), _random_generator(random_seed) {
    if (log) {
      log->write_status(
          "Created UniformRandomVoronoiGeneratorDistribution with ",
          _number_of_positions, " positions and random seed ", random_seed,
          ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  UniformRandomVoronoiGeneratorDistribution(ParameterFile &params,
                                            Log *log = nullptr)
      : UniformRandomVoronoiGeneratorDistribution(
            Box(params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_anchor", "[0. m, 0. m, 0. m]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_sides", "[1. m, 1. m, 1. m]")),
            params.get_value< unsigned int >("densitygrid:voronoi_generator_"
                                             "distribution:number_of_positions",
                                             100),
            params.get_value< int >(
                "densitygrid:voronoi_generator_distribution:random_seed", 42),
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
   * @brief Get a uniform random generator position.
   *
   * @return Uniform random generator position (in m).
   */
  virtual CoordinateVector<> get_position() {
    ++_current_number;
    cmac_assert(_current_number <= _number_of_positions);
    CoordinateVector<> pos;
    pos[0] =
        _box.get_anchor().x() +
        _random_generator.get_uniform_random_double() * _box.get_sides().x();
    pos[1] =
        _box.get_anchor().y() +
        _random_generator.get_uniform_random_double() * _box.get_sides().y();
    pos[2] =
        _box.get_anchor().z() +
        _random_generator.get_uniform_random_double() * _box.get_sides().z();
    return pos;
  }
};

#endif // UNIFORMRANDOMVORONOIGENERATORDISTRIBUTION_HPP
