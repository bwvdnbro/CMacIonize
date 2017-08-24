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
 * @file SPHArrayDensityGridWriter.hpp
 *
 * @brief DensityGridWriter implementation that maps back to an SPH particle
 * data array.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHARRAYDENSITYGRIDWRITER_HPP
#define SPHARRAYDENSITYGRIDWRITER_HPP

#include "DensityGridWriter.hpp"

class Octree;

/**
 * @brief DensityGridWriter implementation that maps back to an SPH particle
 * data array.
 */
class SPHArrayDensityGridWriter : public DensityGridWriter {
private:
  /*! @brief Neutral fractions on the positions of the SPH particles. */
  std::vector< double > _neutral_fractions;

  /*! @brief Octree that contains the SPH particles. */
  Octree *_octree;

  static double cubic_spline_kernel(double u, double h);

public:
  SPHArrayDensityGridWriter();

  void reset(const size_t numpart, Octree *octree);

  void fill_array(double *nH);
  void fill_array(float *nH);

  virtual void write(DensityGrid &grid, unsigned int iteration,
                     ParameterFile &params, double time);
};

#endif // SPHARRAYDENSITYGRIDWRITER_HPP
