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
 * @file SPHArrayInterface.hpp
 *
 * @brief DensityFunction and DensityGridWriter implementations that are coupled
 * to an underlying SPH simulation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPHARRAYINTERFACE_HPP
#define SPHARRAYINTERFACE_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"
#include "DensityGridWriter.hpp"
#include "Octree.hpp"
#include "TimeLogger.hpp"

/**
 * @brief DensityFunction and DensityGridWriter implementations that are coupled
 * to an underlying SPH simulation.
 */
class SPHArrayInterface : public DensityFunction, public DensityGridWriter {
private:
  /*! @brief Length unit used in the input arrays (in m). */
  const double _unit_length_in_SI;

  /*! @brief Mass unit used in the input arrays (in kg). */
  const double _unit_mass_in_SI;

  /*! @brief Periodicity flag. */
  const bool _is_periodic;

  /*! @brief Simulation box, only initialized if the box is periodic (if the box
   *  is not periodic, the components of the Box will all be zero). */
  Box<> _box;

  /*! @brief Positions of the SPH particles in the snapshot (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Masses of the SPH particles in the snapshot (in kg). */
  std::vector< double > _masses;

  /*! @brief Smoothing lengths of the SPH particles in the snapshot (in m). */
  std::vector< double > _smoothing_lengths;

  /*! @brief Neutral fractions on the positions of the SPH particles. */
  std::vector< double > _neutral_fractions;

  /*! @brief Grid of pre-computed cell densities. */
  std::vector< std::vector< std::vector< double > > > _density_values;

  /*! @brief Octree used to speed up neighbour searching. */
  Octree *_octree;

  /*! @brief Octree used to speed up neighbour searching. */
  TimeLogger time_log;

public:
  SPHArrayInterface(const double unit_length_in_SI,
                    const double unit_mass_in_SI);
  SPHArrayInterface(const double unit_length_in_SI,
                    const double unit_mass_in_SI, const double *box_anchor,
                    const double *box_sides);
  SPHArrayInterface(const double unit_length_in_SI,
                    const double unit_mass_in_SI, const float *box_anchor,
                    const float *box_sides);
  ~SPHArrayInterface();

  // DensityFunction functionality

  void reset(const double *x, const double *y, const double *z, const double *h,
             const double *m, const size_t npart);
  void reset(const double *x, const double *y, const double *z, const float *h,
             const float *m, const size_t npart);
  void reset(const float *x, const float *y, const float *z, const float *h,
             const float *m, const size_t npart);

  Octree *get_octree();

  // DensityMapping get_dens_map(){return _dens_map;}

  virtual void initialize();
  virtual DensityValues operator()(const Cell &cell) const;

  void gridding();

  double gridded_integral(double phi, double r0_old, double R_0_old,
                          double h_old) const;

  static double full_integral(double phi, double r0, double R_0, double h);

  double mass_contribution(const Cell &cell, const CoordinateVector<> particle,
                           const double h) const;

  // DensityGridWriter functionality

  void fill_array(double *nH);
  void fill_array(float *nH);

  virtual void write(DensityGrid &grid, uint_fast32_t iteration,
                     ParameterFile &params, double time,
                     const InternalHydroUnits *hydro_units = nullptr);

  virtual double get_gridded_density_value(int i, int j, int k) const;

  virtual void write(DensitySubGridCreator< DensitySubGrid > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.);

  /**
   * @brief Write a snapshot for a split grid with hydro.
   *
   * @param grid_creator Grid.
   * @param counter Counter value to add to the snapshot file name.
   * @param params ParameterFile containing the run parameters that should be
   * written to the file.
   * @param time Simulation time (in s).
   */
  virtual void write(DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.) {}
};

#endif // SPHARRAYINTERFACE_HPP
