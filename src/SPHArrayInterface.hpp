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
#include "CubicSplineKernel.hpp"
#include "DensityFunction.hpp"
#include "DensityGridWriter.hpp"
#include "Octree.hpp"
#include "TimeLogger.hpp"

/**
 * @brief Types of mapping that can be used to map SPH particles to grid cells.
 */
enum SPHArrayMappingType {
  /*! @brief Densities are assigned based on the reciprocal volume of the
   *  corresponding Voronoi cell, no actual mapping is performed. */
  SPHARRAY_MAPPING_M_OVER_V = 0,
  /*! @brief Densities are assigned based on the SPH density at the location
   *  of the centroid of the Voronoi cells. */
  SPHARRAY_MAPPING_CENTROID,
  /*! @brief Densities are assigned using the Petkova et al. (2018) algorithm
   *  that explicitly conserves mass. */
  SPHARRAY_MAPPING_PETKOVA
};

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

  /*! @brief Type of mapping to use to convert SPH densities to cell
   *  densities. */
  const SPHArrayMappingType _mapping_type;

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

  /*! @brief Time log used to register time consumption in various parts of the
   *  algorithm. */
  TimeLogger _time_log;

  /**
   * @brief Get the SPHArrayMappingType corresponding to the given type name
   * string.
   *
   * @param type_name Type name string (valid options are: M_over_V, centroid,
   * Petkova).
   * @return Corresponding SPHArrayMappingType.
   */
  inline static SPHArrayMappingType
  get_mapping_type(const std::string type_name) {

    if (type_name == "M_over_V") {
      return SPHARRAY_MAPPING_M_OVER_V;
    } else if (type_name == "centroid") {
      return SPHARRAY_MAPPING_CENTROID;
    } else if (type_name == "Petkova") {
      return SPHARRAY_MAPPING_PETKOVA;
    } else {
      cmac_error("Unknown SPHArrayMappingType: \"%s\"!", type_name.c_str());
      return SPHARRAY_MAPPING_PETKOVA;
    }
  }

  /**
   * @brief Functor for the inverse mapping.
   */
  class InverseMappingFunction {
  private:
    /*! @brief Reference to the SPHArrayInterface. */
    SPHArrayInterface &_array_interface;

    /*! @brief Locks protecting the neutral fraction vector entries. */
    std::vector< Lock > &_locks;

  public:
    /**
     * @brief Constructor.
     *
     * @param array_interface Reference to the SPHArrayInterface.
     * @param locks Locks protecting the neutral fraction vector entries.
     */
    InverseMappingFunction(SPHArrayInterface &array_interface,
                           std::vector< Lock > &locks)
        : _array_interface(array_interface), _locks(locks) {}

    /**
     * @brief Do the inverse mapping for a single cell.
     *
     * @param cell DensityGrid::iterator pointing to a single cell in the grid.
     */
    void operator()(DensityGrid::iterator cell) {

      if (_array_interface._mapping_type == SPHARRAY_MAPPING_M_OVER_V) {
        const CoordinateVector<> p = cell.get_cell_midpoint();
        uint_fast32_t closest = _array_interface._octree->get_closest_ngb(p);
        _locks[closest].lock();
        _array_interface._neutral_fractions[closest] =
            cell.get_ionization_variables().get_ionic_fraction(ION_H_n);
        _locks[closest].unlock();
      } else if (_array_interface._mapping_type == SPHARRAY_MAPPING_CENTROID) {
        const CoordinateVector<> p = cell.get_cell_midpoint();
        const std::vector< uint_fast32_t > ngbs =
            _array_interface._octree->get_ngbs(p);
        const unsigned int numngbs = ngbs.size();
        double cell_mass = 0.;
        for (unsigned int i = 0; i < numngbs; ++i) {
          const unsigned int index = ngbs[i];
          double r;
          if (!_array_interface._box.get_sides().x()) {
            r = (p - _array_interface._positions[index]).norm();
          } else {
            r = _array_interface._box
                    .periodic_distance(p, _array_interface._positions[index])
                    .norm();
          }
          const double h = _array_interface._smoothing_lengths[index];
          const double u = r / h;
          const double m = _array_interface._masses[index];
          const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
          cell_mass += splineval;
        }

        for (unsigned int i = 0; i < numngbs; ++i) {
          const unsigned int index = ngbs[i];
          double r;
          if (!_array_interface._box.get_sides().x()) {
            r = (p - _array_interface._positions[index]).norm();
          } else {
            r = _array_interface._box
                    .periodic_distance(p, _array_interface._positions[index])
                    .norm();
          }
          const double h = _array_interface._smoothing_lengths[index];
          const double u = r / h;
          const double m = _array_interface._masses[index];
          const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
          _locks[index].lock();
          _array_interface._neutral_fractions[index] -=
              splineval / cell_mass *
              (1. -
               cell.get_ionization_variables().get_ionic_fraction(ION_H_n));
          _locks[index].unlock();
        }
      } else if (_array_interface._mapping_type == SPHARRAY_MAPPING_PETKOVA) {
        const CoordinateVector<> position = cell.get_cell_midpoint();
        // Find the vertex that is furthest away from the cell midpoint.
        std::vector< Face > face_vector = cell.get_faces();
        double radius = 0.0;
        for (unsigned int i = 0; i < face_vector.size(); i++) {
          for (Face::Vertices j = face_vector[i].first_vertex();
               j != face_vector[i].last_vertex(); ++j) {
            double distance = (j.get_position() - position).norm();
            if (distance > radius)
              radius = distance;
          }
        }

        // Find the neighbours that are contained inside of a sphere of centre
        // the cell midpoint and radius given by the distance to the furthest
        // vertex.
        std::vector< uint_fast32_t > ngbs =
            _array_interface._octree->get_ngbs_sphere(position, radius);
        const unsigned int numngbs = ngbs.size();

        double cell_mass = 0.;
        std::vector< double > denseval;

        // Loop over all the neighbouring particles and calculate their mass
        // contributions.
        for (unsigned int i = 0; i < numngbs; i++) {
          const unsigned int index = ngbs[i];
          const double h = _array_interface._smoothing_lengths[index] / 2.0;
          const CoordinateVector<> particle =
              _array_interface._positions[index];
          denseval.push_back(
              _array_interface.mass_contribution(cell, particle, h) *
              _array_interface._masses[index]);
          // denseval.push_back(CubicSplineKernel::kernel_evaluate(0.5, h) *
          // _masses[index]);
          cell_mass += denseval[i];
        }

        for (unsigned int i = 0; i < numngbs; i++) {
          const unsigned int index = ngbs[i];
          _locks[index].lock();
          _array_interface._neutral_fractions[index] -=
              denseval[i] / cell_mass *
              (1. -
               cell.get_ionization_variables().get_ionic_fraction(ION_H_n));
          _locks[index].unlock();
        }
      }
    }
  };

public:
  SPHArrayInterface(const double unit_length_in_SI,
                    const double unit_mass_in_SI,
                    const std::string mapping_type);
  SPHArrayInterface(const double unit_length_in_SI,
                    const double unit_mass_in_SI, const double *box_anchor,
                    const double *box_sides, const std::string mapping_type);
  SPHArrayInterface(const double unit_length_in_SI,
                    const double unit_mass_in_SI, const float *box_anchor,
                    const float *box_sides, const std::string mapping_type);
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
  virtual DensityValues operator()(const Cell &cell);

  void gridding();

  double gridded_integral(const double phi, const double cosphi,
                          const double r0_old, const double R_0_old,
                          const double h_old) const;

  static double full_integral(const double phi, const double cosphi,
                              const double r0, const double R_0,
                              const double h);

  double mass_contribution(const Cell &cell, const CoordinateVector<> particle,
                           const double h) const;

  // DensityGridWriter functionality

  void fill_array(double *nH);
  void fill_array(float *nH);

  virtual void write(DensityGrid &grid, uint_fast32_t iteration,
                     ParameterFile &params, double time,
                     const InternalHydroUnits *hydro_units = nullptr);

  /**
   * @brief Get the gridded density value for the given indices.
   *
   * @param i First index.
   * @param j Second index.
   * @param k Third index.
   * @return Corresponding gridded density value.
   */
  inline double get_gridded_density_value(const uint_fast32_t i,
                                          const uint_fast32_t j,
                                          const uint_fast32_t k) const {
    return _density_values[i][j][k];
  }

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
