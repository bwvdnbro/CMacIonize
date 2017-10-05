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
 * @file SPHArrayDensityGridWriter.cpp
 *
 * @brief SPHArrayDensityGridWriter implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "SPHArrayDensityGridWriter.hpp"
#include "DensityGrid.hpp"
#include "Octree.hpp"

/**
 * @brief Constructor.
 */
SPHArrayDensityGridWriter::SPHArrayDensityGridWriter()
    : DensityGridWriter("", nullptr), _octree(nullptr) {}

/**
 * @brief Reset the internal state of the writer.
 *
 * @param numpart Number of SPH particles in the underlying SPH particle
 * distribution.
 * @param octree Octree to use for neighbour searches.
 */
void SPHArrayDensityGridWriter::reset(const size_t numpart, Octree *octree) {
  _neutral_fractions.resize(numpart, 0.);
  _octree = octree;
}

/**
 * @brief Fill the given array with the remapped neutral fractions.
 *
 * Double precision version.
 *
 * @param nH Array to fill.
 */
void SPHArrayDensityGridWriter::fill_array(double *nH) {
  for (size_t i = 0; i < _neutral_fractions.size(); ++i) {
    nH[i] = _neutral_fractions[i];
  }
}

/**
 * @brief Fill the given array with the remapped neutral fractions.
 *
 * Single precision version.
 *
 * @param nH Array to fill.
 */
void SPHArrayDensityGridWriter::fill_array(float *nH) {
  for (size_t i = 0; i < _neutral_fractions.size(); ++i) {
    nH[i] = _neutral_fractions[i];
  }
}

/**
 * @brief Map the state of the grid back to the SPH particle distribution.
 *
 * @param grid DensityGrid to write out.
 * @param iteration Iteration number to use in the snapshot file name(s).
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void SPHArrayDensityGridWriter::write(DensityGrid &grid,
                                      uint_fast32_t iteration,
                                      ParameterFile &params, double time) {
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    const CoordinateVector<> p = it.get_cell_midpoint();
    uint_fast32_t closest = _octree->get_closest_ngb(p);
    _neutral_fractions[closest] =
        it.get_ionization_variables().get_ionic_fraction(ION_H_n);
  }
}
