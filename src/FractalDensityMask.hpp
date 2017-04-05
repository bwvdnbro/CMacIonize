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
 * @file FractalDensityMask.hpp
 *
 * @brief Fractal density mask that redistributes the density in the given grid
 * according to a fractal distribution.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FRACTALDENSITYMASK_HPP
#define FRACTALDENSITYMASK_HPP

#include "Box.hpp"
#include "DensityGrid.hpp"
#include "RandomGenerator.hpp"

#include <iostream>
#include <vector>

/**
 * @brief Fractal density mask that redistributes the density in the given grid
 * according to a fractal distribution.
 *
 * The algorithm used to generate the fractal density field is based on the
 * algorithm described in Elmegreen, B., 1997, ApJ, 477, 196.
 */
class FractalDensityMask {
private:
  /*! @brief Box containing the fractal distribution (in m). */
  Box _box;

  /*! @brief Resolution of the grid containing the distribution. */
  CoordinateVector< unsigned int > _resolution;

  /*! @brief Grid containing the distribution. */
  std::vector< std::vector< std::vector< double > > > _distribution;

  /**
   * @brief (Recursively) construct a fractal grid with the given number of
   * levels, and given number of points per level and fractal length scale.
   *
   * @param fractal_distribution std::vector containing the fractal distribution
   * (is updated).
   * @param random_generator RandomGenerator used to generate random positions.
   * @param N Number of points per level.
   * @param L Fractal length scale.
   * @param num_level Number of levels we want to create.
   * @param current_position Offset position for points on the current level
   * (in internal fractal units).
   * @param current_level Current level.
   */
  static void make_fractal_grid(
      std::vector< std::vector< std::vector< double > > > &fractal_distribution,
      RandomGenerator &random_generator, unsigned int N, double L,
      unsigned int num_level,
      CoordinateVector<> current_position = CoordinateVector<>(0.),
      unsigned int current_level = 0) {
    if (current_level < num_level) {
      for (unsigned int i = 0; i < N; ++i) {
        CoordinateVector<> x_level = current_position;
        x_level[0] += 2. *
                      (random_generator.get_uniform_random_double() - 0.5) /
                      std::pow(L, current_level + 1);
        x_level[1] += 2. *
                      (random_generator.get_uniform_random_double() - 0.5) /
                      std::pow(L, current_level + 1);
        x_level[2] += 2. *
                      (random_generator.get_uniform_random_double() - 0.5) /
                      std::pow(L, current_level + 1);
        make_fractal_grid(fractal_distribution, random_generator, N, L,
                          num_level, x_level, current_level + 1);
      }
    } else {
      // current_position now contains coordinates in the range [-1/L, 1/L]
      // map them to the range [0., 1.]
      current_position[0] *= 0.5 * L;
      current_position[1] *= 0.5 * L;
      current_position[2] *= 0.5 * L;
      current_position[0] += 0.5;
      current_position[1] += 0.5;
      current_position[2] += 0.5;
      // some coordinates may have fallen outside the range, periodically map
      // them to coordinates inside
      if (current_position.x() < 0.) {
        current_position[0] += 1.;
      }
      if (current_position.x() >= 1.) {
        current_position[0] -= 1.;
      }
      if (current_position.y() < 0.) {
        current_position[1] += 1.;
      }
      if (current_position.y() >= 1.) {
        current_position[1] -= 1.;
      }
      if (current_position.z() < 0.) {
        current_position[2] += 1.;
      }
      if (current_position.z() >= 1.) {
        current_position[2] -= 1.;
      }

      cmac_assert(current_position.x() >= 0. && current_position.x() < 1.);
      cmac_assert(current_position.y() >= 0. && current_position.y() < 1.);
      cmac_assert(current_position.z() >= 0. && current_position.z() < 1.);

      unsigned int ix = current_position.x() * fractal_distribution.size();
      unsigned int iy = current_position.y() * fractal_distribution[ix].size();
      unsigned int iz =
          current_position.z() * fractal_distribution[ix][iy].size();

      fractal_distribution[ix][iy][iz] += 1.;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box in which the particles are generated.
   * @param resolution Resolution of the grid on which the fractal distribution
   * is sampled.
   * @param numpart Minimum number of particles to generate. The actual number
   * can be higher due to round off in intermediate operations.
   * @param seed Seed for the random generator. Different seeds will lead to
   * different fractal distributions.
   * @param fractal_dimension Fractal dimension we want to sample.
   * @param num_level Number of levels of the fractal structure. According to
   * Elemegreen (1997), this should be a number in the range 3-6 for the ISM.
   */
  FractalDensityMask(Box box, CoordinateVector< unsigned int > resolution,
                     unsigned int numpart, int seed, double fractal_dimension,
                     unsigned int num_level)
      : _box(box), _resolution(resolution) {
    RandomGenerator random_generator(seed);

    // allocate the grid
    _distribution.resize(resolution.x());
    for (unsigned int ix = 0; ix < _resolution.x(); ++ix) {
      _distribution[ix].resize(resolution.y());
      for (unsigned int iy = 0; iy < _resolution.y(); ++iy) {
        _distribution[ix][iy].resize(_resolution.z(), 0.);
      }
    }

    // we will use equal values for the number of points (N) per level
    unsigned int N = std::ceil(std::pow(numpart, 1. / num_level));
    // the fractal length scale (L) is linked to the fractal dimension (D) and
    // the number of points per level by the formula D = log10(N)/log10(L)
    double L = std::pow(10., std::log10(N) / fractal_dimension);

    // recursively construct the fractal grid
    make_fractal_grid(_distribution, random_generator, N, L, num_level);
  }

  /**
   * @brief Apply the mask to the given DensityGrid, using the given fractal
   * fraction.
   *
   * @param grid DensityGrid to apply the mask to.
   * @param fractal_fraction Fraction of the cell density that is set by the
   * fractal distribution.
   */
  void apply(DensityGrid &grid, double fractal_fraction) const {
    double smooth_fraction = 1. - fractal_fraction;
    double Ntot = 0.;
    double Nsmooth = 0.;
    double Nfractal = 0.;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      double ncell = it.get_number_density();
      double Ncell = ncell * it.get_volume();
      Ntot += Ncell;
      CoordinateVector<> midpoint = it.get_cell_midpoint();
      unsigned int ix = (midpoint.x() - _box.get_anchor().x()) /
                        _box.get_sides().x() * _resolution.x();
      unsigned int iy = (midpoint.y() - _box.get_anchor().y()) /
                        _box.get_sides().y() * _resolution.y();
      unsigned int iz = (midpoint.z() - _box.get_anchor().z()) /
                        _box.get_sides().z() * _resolution.z();
      Nsmooth += smooth_fraction * Ncell;
      Nfractal += fractal_fraction * Ncell * _distribution[ix][iy][iz];
    }

    double fractal_norm = (Ntot - Nsmooth) / Nfractal;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      double ncell = it.get_number_density();
      CoordinateVector<> midpoint = it.get_cell_midpoint();
      unsigned int ix = (midpoint.x() - _box.get_anchor().x()) /
                        _box.get_sides().x() * _resolution.x();
      unsigned int iy = (midpoint.y() - _box.get_anchor().y()) /
                        _box.get_sides().y() * _resolution.y();
      unsigned int iz = (midpoint.z() - _box.get_anchor().z()) /
                        _box.get_sides().z() * _resolution.z();
      double nsmooth = smooth_fraction * ncell;
      double nfractal =
          fractal_fraction * fractal_norm * ncell * _distribution[ix][iy][iz];
      it.set_number_density(nsmooth + nfractal);
    }
  }
};

#endif // FRACTALDENSITYMASK_HPP
