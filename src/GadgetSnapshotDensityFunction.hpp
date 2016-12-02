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
 * @file GadgetSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density field from a Gadget HDF5 snapshot
 * file: header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef GADGETSNAPSHOTDENSITYFUNCTION_HPP
#define GADGETSNAPSHOTDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"
#include "Octree.hpp"
#include <string>
#include <vector>

class Log;
class ParameterFile;

/**
 * @brief DensityFunction that reads a density field from a Gadget snapshot.
 */
class GadgetSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Simulation box, only initialized if the box is periodic (if the box
   *  is not periodic, the components of the Box will all be zero). */
  Box _box;

  /*! @brief Positions of the SPH particles in the snapshot (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Masses of the SPH particles in the snapshot (in kg). */
  std::vector< double > _masses;

  /*! @brief Smoothing lengths of the SPH particles in the snapshot (in m). */
  std::vector< double > _smoothing_lengths;

  /*! @brief Densities of the SPH particles in the snapshot (in kg m^-3). Note
   *  that these are not the densities we return; they are only used to do the
   *  correct kernel interpolation for the temperature. */
  std::vector< double > _densities;

  /*! @brief Temperatures of the SPH particles in the snapshot (in K). */
  std::vector< double > _temperatures;

  /*! @brief Octree used to speed up neighbour searching. */
  Octree *_octree;

  /*! @brief Log to write logging info to. */
  Log *_log;

  double cubic_spline_kernel(double u, double h);

public:
  GadgetSnapshotDensityFunction(std::string name,
                                bool fallback_periodic = false,
                                double fallback_unit_length_in_SI = 0.,
                                double fallback_unit_mass_in_SI = 0.,
                                double fallback_unit_temperature_in_SI = 0.,
                                Log *log = nullptr);

  GadgetSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual ~GadgetSnapshotDensityFunction();

  virtual DensityValues operator()(CoordinateVector<> position);

  double get_total_hydrogen_number();
};

#endif // GADGETSNAPSHOTDENSITYFUNCTION_HPP
