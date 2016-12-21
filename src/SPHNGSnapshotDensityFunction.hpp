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
 * @file SPHNGSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction implementation that reads a density field from an
 * SPHNG snapshot file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHNGSNAPSHOTDENSITYFUNCTION_HPP
#define SPHNGSNAPSHOTDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"

class Log;
class Octree;
class ParameterFile;

/**
 * @brief DensityFunction implementation that reads a density field from an
 * SPHNG snapshot file.
 */
class SPHNGSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Positions of the SPH particles in the snapshot (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Masses of the SPH particles in the snapshot (in kg). */
  std::vector< double > _masses;

  /*! @brief Smoothing lengths of the SPH particles in the snapshot (in m). */
  std::vector< double > _smoothing_lengths;

  /*! @brief Octree used to speed up neighbour finding. */
  Octree *_octree;

  static double kernel(const double q, const double h);

public:
  SPHNGSnapshotDensityFunction(std::string filename, Log *log = nullptr);

  SPHNGSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  ~SPHNGSnapshotDensityFunction();

  CoordinateVector<> get_position(unsigned int index);
  double get_mass(unsigned int index);
  double get_smoothing_length(unsigned int index);

  DensityValues operator()(CoordinateVector<> position) const;
};

#endif // SPHNGSNAPSHOTDENSITYFUNCTION_HPP
