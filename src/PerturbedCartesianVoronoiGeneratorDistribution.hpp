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
 * @file PerturbedCartesianVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution implementation that generates a Cartesian
 * grid and then adds small perturbations to the generated positions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PERTURBEDCARTESIANVORONOIGENERATORDISTRIBUTION_HPP
#define PERTURBEDCARTESIANVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "VoronoiGeneratorDistribution.hpp"

#include <vector>

/**
 * @brief VoronoiGeneratorDistribution implementation that generates a Cartesian
 * grid and then adds small perturbations to the generated positions.
 */
class PerturbedCartesianVoronoiGeneratorDistribution
    : public VoronoiGeneratorDistribution {
private:
  /*! @brief Generator positions (in m). */
  std::vector< CoordinateVector<> > _generator_positions;

  /*! @brief Index of the next position that will be returned. */
  unsigned int _current_number;

public:
  /**
   * @brief Constructor.
   *
   * @param simulation_box Simulation box (in m).
   * @param resolution Number of cells in each dimensions for the base Cartesian
   * grid.
   * @param random_seed Seed for the random number generator.
   * @param amplitude Amplitude of the random perturbations (in m).
   * @param log Log to write logging info to.
   */
  PerturbedCartesianVoronoiGeneratorDistribution(
      const Box<> &simulation_box,
      const CoordinateVector< unsigned int > resolution, double random_seed,
      double amplitude, Log *log = nullptr)
      : _current_number(0) {
    RandomGenerator rg(random_seed);
    _generator_positions.resize(resolution.x() * resolution.y() *
                                resolution.z());
    for (unsigned int ix = 0; ix < resolution.x(); ++ix) {
      const unsigned int index_x = ix * resolution.y() * resolution.z();
      const double cx =
          simulation_box.get_anchor().x() +
          (ix + 0.5) * simulation_box.get_sides().x() / resolution.x();
      for (unsigned int iy = 0; iy < resolution.y(); ++iy) {
        const unsigned int index_y = index_x + iy * resolution.z();
        const double cy =
            simulation_box.get_anchor().y() +
            (iy + 0.5) * simulation_box.get_sides().y() / resolution.y();
        for (unsigned int iz = 0; iz < resolution.z(); ++iz) {
          const unsigned int index_z = index_y + iz;
          const double cz =
              simulation_box.get_anchor().z() +
              (iz + 0.5) * simulation_box.get_sides().z() / resolution.z();
          _generator_positions[index_z][0] =
              cx + amplitude * (2. * rg.get_uniform_random_double() - 1.);
          _generator_positions[index_z][1] =
              cy + amplitude * (2. * rg.get_uniform_random_double() - 1.);
          _generator_positions[index_z][2] =
              cz + amplitude * (2. * rg.get_uniform_random_double() - 1.);
        }
      }
    }

    if (log) {
      log->write_status(
          "Created PerturbedCartesianVoronoiGeneratorDistribution with ",
          resolution.x(), "x", resolution.y(), "x", resolution.z(),
          " cells with a random perturbation amplitude of ", amplitude,
          " m (random seed is ", random_seed, ").");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - number of cells: Number of cells for the basic unperturbed Cartesian
   *    grid (default: [10, 10, 10])
   *  - random seed: Random seed used to initialize the random generator that
   *    generates the random perturbations (default: 42)
   *  - amplitude: Amplitude of the random displacements (default: 0.01 m)
   *
   * @param simulation_box Simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  PerturbedCartesianVoronoiGeneratorDistribution(const Box<> &simulation_box,
                                                 ParameterFile &params,
                                                 Log *log = nullptr)
      : PerturbedCartesianVoronoiGeneratorDistribution(
            simulation_box,
            params.get_value< CoordinateVector< unsigned int > >(
                "DensityGrid:VoronoiGeneratorDistribution:number of cells",
                CoordinateVector< unsigned int >(10, 10, 10)),
            params.get_value< int >(
                "DensityGrid:VoronoiGeneratorDistribution:random seed", 42),
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityGrid:VoronoiGeneratorDistribution:amplitude", "0.01 m"),
            log) {}

  /**
   * @brief Get the number of positions that this distribution generates.
   *
   * @return Number of generated positions.
   */
  virtual unsigned int get_number_of_positions() const {
    return _generator_positions.size();
  }

  /**
   * @brief Get the next position.
   *
   * @return Next generator position (in m).
   */
  virtual CoordinateVector<> get_position() {
    cmac_assert(_current_number < _generator_positions.size());
    CoordinateVector<> pos = _generator_positions[_current_number];
    ++_current_number;
    return pos;
  }
};

#endif // PERTURBEDCARTESIANVORONOIGENERATORDISTRIBUTION_HPP
