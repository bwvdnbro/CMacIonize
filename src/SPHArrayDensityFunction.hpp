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
 * @file SPHArrayDensityFunction.hpp
 *
 * @brief DensityFunction that is based on an SPH density field that is provided
 * from external data arrays.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHARRAYDENSITYFUNCTION_HPP
#define SPHARRAYDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"
#include "Octree.hpp"

/**
 * @brief DensityFunction that is based on an SPH density field that is provided
 * from external data arrays.
 */
class SPHArrayDensityFunction : public DensityFunction {
private:
  /*! @brief Length unit used in the input arrays (in m). */
  const double _unit_length_in_SI;

  /*! @brief Mass unit used in the input arrays (in kg). */
  const double _unit_mass_in_SI;

  /*! @brief Periodicity flag. */
  const bool _periodic;

  /*! @brief Simulation box, only initialized if the box is periodic (if the box
   *  is not periodic, the components of the Box will all be zero). */
  Box<> _box;

  /*! @brief Positions of the SPH particles in the snapshot (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Masses of the SPH particles in the snapshot (in kg). */
  std::vector< double > _masses;

  /*! @brief Smoothing lengths of the SPH particles in the snapshot (in m). */
  std::vector< double > _smoothing_lengths;

  /*! @brief Octree used to speed up neighbour searching. */
  Octree *_octree;

  static double cubic_spline_kernel(double u, double h);

public:
  SPHArrayDensityFunction(const double *x, const double *y, const double *z,
                          const double *h, const double *m, const size_t npart,
                          const double unit_length_in_SI,
                          const double unit_mass_in_SI, const bool periodic,
                          const double *box_anchor, const double *box_sides);
  SPHArrayDensityFunction(const double *x, const double *y, const double *z,
                          const float *h, const float *m, const size_t npart,
                          const double unit_length_in_SI,
                          const double unit_mass_in_SI, const bool periodic,
                          const double *box_anchor, const double *box_sides);
  SPHArrayDensityFunction(const float *x, const float *y, const float *z,
                          const float *h, const float *m, const size_t npart,
                          const double unit_length_in_SI,
                          const double unit_mass_in_SI, const bool periodic,
                          const float *box_anchor, const float *box_sides);
  ~SPHArrayDensityFunction();

  void reset(const double *x, const double *y, const double *z, const double *h,
             const double *m, const size_t npart);
  void reset(const double *x, const double *y, const double *z, const float *h,
             const float *m, const size_t npart);
  void reset(const float *x, const float *y, const float *z, const float *h,
             const float *m, const size_t npart);

  virtual void initialize();
  virtual DensityValues operator()(const Cell &cell) const;
};

#endif // SPHARRAYDENSITYFUNCTION_HPP
