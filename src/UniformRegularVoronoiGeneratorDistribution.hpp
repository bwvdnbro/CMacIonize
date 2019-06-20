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
 * @file UniformRegularVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution that sets up generators on a regular
 * Cartesian grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UNIFORMREGULARVORONOIGENERATORDISTRIBUTION_HPP
#define UNIFORMREGULARVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "VoronoiGeneratorDistribution.hpp"

/**
 * @brief VoronoiGeneratorDistribution that sets up generators on a regular
 * Cartesian grid.
 */
class UniformRegularVoronoiGeneratorDistribution
    : public VoronoiGeneratorDistribution {
private:
  /*! @brief Box containing the generators (in m). */
  const Box<> _box;

  /*! @brief Resolution of the generator grid. */
  CoordinateVector< uint_fast32_t > _resolution;

  /*! @brief Side lengths of a single cube of the regular grid (in m). */
  CoordinateVector<> _sidelength;

  /*! @brief Indices of the next generator to return. */
  CoordinateVector< uint_fast32_t > _next_index;

public:
  /**
   * @brief Constructor.
   *
   * @param simulation_box Simulation box (in m).
   * @param resolution Resolution of the generator grid.
   * @param log Log to write logging info to.
   */
  UniformRegularVoronoiGeneratorDistribution(
      const Box<> &simulation_box, CoordinateVector< uint_fast32_t > resolution,
      Log *log = nullptr)
      : _box(simulation_box), _resolution(resolution) {

    _sidelength[0] = _box.get_sides().x() / _resolution.x();
    _sidelength[1] = _box.get_sides().y() / _resolution.y();
    _sidelength[2] = _box.get_sides().z() / _resolution.z();

    if (log) {
      log->write_status(
          "Constructed "
          "UniformRegularVoronoiGeneratorDistribution with a "
          "resolution of ",
          _resolution.x(), "x", _resolution.y(), "x", _resolution.z(),
          " generators, and a regular grid cell size of ", _sidelength.x(),
          " m x ", _sidelength.y(), " m x ", _sidelength.z(), " m.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - resolution: Resolution of the regular generator grid (default: [32, 32,
   *    32])
   *
   * @param simulation_box Simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  UniformRegularVoronoiGeneratorDistribution(const Box<> &simulation_box,
                                             ParameterFile &params,
                                             Log *log = nullptr)
      : UniformRegularVoronoiGeneratorDistribution(
            simulation_box,
            params.get_value< CoordinateVector< uint_fast32_t > >(
                "DensityGrid:VoronoiGeneratorDistribution:resolution",
                CoordinateVector< uint_fast32_t >(32)),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~UniformRegularVoronoiGeneratorDistribution() {}

  /**
   * @brief Get the number of generators in the grid.
   *
   * @return Number of generators in the grid.
   */
  virtual generatornumber_t get_number_of_positions() const {
    return _resolution.x() * _resolution.y() * _resolution.z();
  }

  /**
   * @brief Get the next generator position.
   *
   * @return Next generator position (in m).
   */
  virtual CoordinateVector<> get_position() {
    if (_next_index.x() < _resolution.x()) {
      CoordinateVector<> pos;
      pos[0] =
          (_next_index.x() + 0.5) * _sidelength.x() + _box.get_anchor().x();
      pos[1] =
          (_next_index.y() + 0.5) * _sidelength.y() + _box.get_anchor().y();
      pos[2] =
          (_next_index.z() + 0.5) * _sidelength.z() + _box.get_anchor().z();

      ++_next_index[2];
      if (_next_index.z() == _resolution.z()) {
        _next_index[2] = 0;
        ++_next_index[1];
        if (_next_index.y() == _resolution.y()) {
          _next_index[1] = 0;
          ++_next_index[0];
        }
      }

      return pos;
    } else {
      cmac_error("More generators requested than contained in the "
                 "UniformRegularVoronoiGeneratorDistribution!");
      return CoordinateVector<>();
    }
  }
};

#endif // UNIFORMREGULARVORONOIGENERATORDISTRIBUTION_HPP
