/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2017 Maya Petkova (map32@st-andrews.ac.uk)
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
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#ifndef SPHNGSNAPSHOTDENSITYFUNCTION_HPP
#define SPHNGSNAPSHOTDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"

#include <cinttypes>
#include <vector>

class Log;
class Octree;
class ParameterFile;

/**
 * @brief DensityFunction implementation that reads a density field from an
 * SPHNG snapshot file.
 */
class SPHNGSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Use the new mapping algorithm? */
  const bool _use_new_algorithm;

  /*! @brief Positions of the SPH particles in the snapshot (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Masses of the SPH particles in the snapshot (in kg). */
  std::vector< double > _masses;

  /*! @brief Smoothing lengths of the SPH particles in the snapshot (in m). */
  std::vector< double > _smoothing_lengths;

  /*! @brief Box containing all particles (in m). */
  Box<> _partbox;

  /*! @brief Octree used to speed up neighbour finding. */
  Octree *_octree;

  /*! @brief Initial temperature of the gas (in K). */
  const double _initial_temperature;

  /*! @brief Number of bins to use when writing particle statistics. */
  const uint_fast32_t _stats_numbin;

  /*! @brief Minimum distance to use for particle statistics. */
  const double _stats_mindist;

  /*! @brief Maximum distance to use for particle statistics. */
  const double _stats_maxdist;

  /*! @brief Name of the file with particle statistics. */
  const std::string _stats_filename;

  /*! @brief Log to write logging info to. */
  Log *_log;

  static double kernel(const double q, const double h);

  static double full_integral(double phi, double r0, double R_0, double h);

  static double mass_contribution(const Cell &cell,
                                  const CoordinateVector<> particle,
                                  const double h);

public:
  SPHNGSnapshotDensityFunction(
      const std::string filename, const double initial_temperature,
      const bool write_stats, const uint_fast32_t stats_numbin,
      const double stats_mindist, const double stats_maxdist,
      const std::string stats_filename, const bool use_new_algorithm = false,
      const bool binary_dump = false, const std::string binary_dump_name = "",
      Log *log = nullptr);

  SPHNGSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  ~SPHNGSnapshotDensityFunction();

  virtual void initialize();

  CoordinateVector<> get_position(uint_fast32_t index);
  double get_mass(uint_fast32_t index);
  double get_smoothing_length(uint_fast32_t index);

  virtual DensityValues operator()(const Cell &cell) const;
};

#endif // SPHNGSNAPSHOTDENSITYFUNCTION_HPP
