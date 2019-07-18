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
 * @file SPHArrayInterface.cpp
 *
 * @brief SPHArrayInterface implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "SPHArrayInterface.hpp"
#include "CubicSplineKernel.hpp"
#include "DensityGrid.hpp"
#include "DensityGridTraversalJobMarket.hpp"
#include "DensityMapping.hpp"
#include "Lock.hpp"
#include <cfloat>

/**
 * @brief Constructor.
 *
 * Non-periodic version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 */
SPHArrayInterface::SPHArrayInterface(const double unit_length_in_SI,
                                     const double unit_mass_in_SI)
    : DensityGridWriter("", nullptr), _unit_length_in_SI(unit_length_in_SI),
      _unit_mass_in_SI(unit_mass_in_SI), _is_periodic(false), _octree(nullptr) {

  _dens_map.gridding();
}

/**
 * @brief Constructor.
 *
 * Periodic double precision version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 * @param box_anchor Coordinates of the left front bottom corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 */
SPHArrayInterface::SPHArrayInterface(const double unit_length_in_SI,
                                     const double unit_mass_in_SI,
                                     const double *box_anchor,
                                     const double *box_sides)
    : DensityGridWriter("", nullptr), _unit_length_in_SI(unit_length_in_SI),
      _unit_mass_in_SI(unit_mass_in_SI), _is_periodic(true), _octree(nullptr) {

  _box.get_anchor()[0] = box_anchor[0] * _unit_length_in_SI;
  _box.get_anchor()[1] = box_anchor[1] * _unit_length_in_SI;
  _box.get_anchor()[2] = box_anchor[2] * _unit_length_in_SI;
  _box.get_sides()[0] = box_sides[0] * _unit_length_in_SI;
  _box.get_sides()[1] = box_sides[1] * _unit_length_in_SI;
  _box.get_sides()[2] = box_sides[2] * _unit_length_in_SI;

  _dens_map.gridding();
}

/**
 * @brief Constructor.
 *
 * Periodic single precision version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 * @param box_anchor Coordinates of the left front bottom corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 */
SPHArrayInterface::SPHArrayInterface(const double unit_length_in_SI,
                                     const double unit_mass_in_SI,
                                     const float *box_anchor,
                                     const float *box_sides)
    : DensityGridWriter("", nullptr), _unit_length_in_SI(unit_length_in_SI),
      _unit_mass_in_SI(unit_mass_in_SI), _is_periodic(true), _octree(nullptr) {

  _box.get_anchor()[0] = box_anchor[0] * _unit_length_in_SI;
  _box.get_anchor()[1] = box_anchor[1] * _unit_length_in_SI;
  _box.get_anchor()[2] = box_anchor[2] * _unit_length_in_SI;
  _box.get_sides()[0] = box_sides[0] * _unit_length_in_SI;
  _box.get_sides()[1] = box_sides[1] * _unit_length_in_SI;
  _box.get_sides()[2] = box_sides[2] * _unit_length_in_SI;

  _dens_map.gridding();
}

/**
 * @brief Destructor.
 *
 * Frees up memory used by the internal Octree.
 */
SPHArrayInterface::~SPHArrayInterface() { delete _octree; }

/**
 * @brief Reset the internal data values.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayInterface::reset(const double *x, const double *y, const double *z,
                              const double *h, const double *m,
                              const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  _neutral_fractions.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Reset the internal data values.
 *
 * This version uses single precision smoothing length and masses.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayInterface::reset(const double *x, const double *y, const double *z,
                              const float *h, const float *m,
                              const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  _neutral_fractions.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Reset the internal data values.
 *
 * This version uses single precision coordinates, smoothing lengths and masses.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayInterface::reset(const float *x, const float *y, const float *z,
                              const float *h, const float *m,
                              const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  _neutral_fractions.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Get a pointer to the internal Octree.
 *
 * @return Pointer to the internal Octree.
 */
Octree *SPHArrayInterface::get_octree() { return _octree; }

/**
 * @brief Initialize the internal Octree.
 */
void SPHArrayInterface::initialize() {
  _octree = new Octree(_positions, _box, _is_periodic);
  _octree->set_auxiliaries(_smoothing_lengths, Octree::max< double >);
  //_dens_map = new DensityMapping();
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues SPHArrayInterface::operator()(const Cell &cell) const {

  const bool _is_petkova_density_mapping = true;
  DensityValues values;

  const CoordinateVector<> position = cell.get_cell_midpoint();
  double density = 0.;

  if(_is_petkova_density_mapping) {
    CoordinateVector<> position = cell.get_cell_midpoint();

    // Find the vertex that is furthest away from the cell midpoint.
    std::vector< Face > face_vector = cell.get_faces();
    double radius = 0.0;
    for (unsigned int i = 0; i < face_vector.size(); i++) {
      for (Face::Vertices j = face_vector[i].first_vertex();
           j != face_vector[i].last_vertex(); ++j) {
        double distance = (j.get_position()-position).norm();
        if (distance > radius)
          radius = distance;
      }
    }

    // Find the neighbours that are contained inside of a sphere of centre the
    // cell midpoint
    // and radius given by the distance to the furthest vertex.
    std::vector< uint_fast32_t > ngbs =
        _octree->get_ngbs_sphere(position, radius);
    const unsigned int numngbs = ngbs.size();

    //printf("ngbs: %u, faces: %lu, radius: %g\n", numngbs, face_vector.size(), radius);

    // Loop over all the neighbouring particles and calculate their mass
    // contributions.
    for (unsigned int i = 0; i < numngbs; i++) {
      const unsigned int index = ngbs[i];
      const double h = _smoothing_lengths[index] / 2.0;
      const CoordinateVector<> particle = _positions[index];
      if(h < 0) printf("h < 0: %g, %u", h, index);
      density += _dens_map.mass_contribution(cell, particle, h) * _masses[index];
    }

    // Divide the cell mass by the cell volume to get density.
    density = density / cell.get_volume();

  } else {
    const std::vector< uint_fast32_t > ngbs = _octree->get_ngbs(position);
    const size_t numngbs = ngbs.size();
    for (size_t i = 0; i < numngbs; ++i) {
      const size_t index = ngbs[i];
      double r;
      if (!_box.get_sides().x()) {
        r = (position - _positions[index]).norm();
      } else {
        r = _box.periodic_distance(position, _positions[index]).norm();
      }
      const double h = _smoothing_lengths[index];
      const double u = r / h;
      const double m = _masses[index];
      const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
      density += splineval;
    }
  }

  values.set_number_density(density / 1.6737236e-27);
  values.set_temperature(8000.);
  values.set_ionic_fraction(ION_H_n, 1.e-6);
  values.set_ionic_fraction(ION_He_n, 1.e-6);
  return values;
}

/**
 * @brief Fill the given array with the remapped neutral fractions.
 *
 * Double precision version.
 *
 * @param nH Array to fill.
 */
void SPHArrayInterface::fill_array(double *nH) {
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
void SPHArrayInterface::fill_array(float *nH) {
  for (size_t i = 0; i < _neutral_fractions.size(); ++i) {
    nH[i] = _neutral_fractions[i];
  }
}

/**
 * @brief Functor for the inverse mapping.
 */
class InverseMappingFunction {
private:
  /*! @brief Reference to the Octree. */
  Octree &_octree;

  /*! @brief Reference to the neutral fraction vector. */
  std::vector< double > &_neutral_fractions;

  /*! @brief Reference to the positions vector. */
  std::vector< CoordinateVector<> > &_positions;

  /*! @brief Reference to the smoothing lengths vector. */
  std::vector< double > &_smoothing_lengths;

  /*! @brief Reference to the masses vector. */
  std::vector< double > &_masses;

  Box<> &_box;

  DensityMapping &_dens_map;

  /*! @brief Locks protecting the neutral fraction vector entries. */
  std::vector< Lock > &_locks;

public:
  /**
   * @brief Constructor.
   *
   * @param octree Reference to the Octree.
   * @param neutral_fractions Reference to the neutral fraction vector.
   * @param positions Reference to the positions vector.
   * @param smoothing_lengths Reference to the smoothing_lengths vector.
   * @param masses Reference to the masses vector.
   * @param box Reference to the simulation box.
   * @param locks Locks protecting the neutral fraction vector entries.
   */
  InverseMappingFunction(Octree &octree,
                         std::vector< double > &neutral_fractions,
			 std::vector< CoordinateVector<> > &positions,
			 std::vector< double > &smoothing_lengths,
			 std::vector< double > &masses,
			 Box<> &box,
			 DensityMapping dens_map,
                         std::vector< Lock > &locks)
      : _octree(octree), _neutral_fractions(neutral_fractions), _positions(positions), _smoothing_lengths(smoothing_lengths), _masses(masses), _box(box), 
	_dens_map(dens_map),_locks(locks) {}
  /**
   * @brief Do the inverse mapping for a single cell.
   *
   * @param cell DensityGrid::iterator pointing to a single cell in the grid.
   */
  inline void operator()(DensityGrid::iterator cell) {

    const bool _is_old_density_mapping = false;
    const bool _is_petkova_density_mapping = true;

    if(_is_old_density_mapping) {
      const CoordinateVector<> p = cell.get_cell_midpoint();
      uint_fast32_t closest = _octree.get_closest_ngb(p);
      _locks[closest].lock();
      _neutral_fractions[closest] =
        cell.get_ionization_variables().get_ionic_fraction(ION_H_n);
      _locks[closest].unlock();
    } else {
      if(_is_petkova_density_mapping) {
	const CoordinateVector<> position = cell.get_cell_midpoint();
      	// Find the vertex that is furthest away from the cell midpoint.
      	std::vector< Face > face_vector = cell.get_faces();
      	double radius = 0.0;
      	for (unsigned int i = 0; i < face_vector.size(); i++) {
          for (Face::Vertices j = face_vector[i].first_vertex();
            	j != face_vector[i].last_vertex(); ++j) {
            double distance = (j.get_position()-position).norm();
            if (distance > radius)
            	radius = distance;
           }
       	 }

      	 // Find the neighbours that are contained inside of a sphere of centre the
      	 // cell midpoint
      	 // and radius given by the distance to the furthest vertex.
      	 std::vector< uint_fast32_t > ngbs =
          	_octree.get_ngbs_sphere(position, radius);
      	 const unsigned int numngbs = ngbs.size();

      	 double cell_mass = 0.;
      	 std::vector< double > denseval;

      	 // Loop over all the neighbouring particles and calculate their mass
      	 // contributions.
      	 for (unsigned int i = 0; i < numngbs; i++) {
           const unsigned int index = ngbs[i];
           const double h = _smoothing_lengths[index] / 2.0;
           const CoordinateVector<> particle = _positions[index];
           denseval.push_back(_dens_map.mass_contribution(cell, particle, h) * _masses[index]);
           cell_mass += denseval[i];
      	 }

      	 for (unsigned int i = 0; i < numngbs; i++) {
           const unsigned int index = ngbs[i];
	   _locks[index].lock();
           _neutral_fractions[index] -= denseval[i]/cell_mass * (1. - cell.get_ionization_variables().get_ionic_fraction(ION_H_n));
	   _locks[index].unlock();
      	 }
	} else {
      	  const CoordinateVector<> p = cell.get_cell_midpoint();
      	  const std::vector< uint_fast32_t > ngbs = _octree.get_ngbs(p);
      	  const unsigned int numngbs = ngbs.size();
      	  double cell_mass = 0.;
      	  for (unsigned int i = 0; i < numngbs; ++i) {
            const unsigned int index = ngbs[i];
	    double r;
            if (!_box.get_sides().x()) {
              r = (p - _positions[index]).norm();
            } else {
              r = _box.periodic_distance(p, _positions[index]).norm();
            }
            const double h = _smoothing_lengths[index];
            const double u = r / h;
            const double m = _masses[index];
            const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
            cell_mass += splineval;
      	  }

      	  for (unsigned int i = 0; i < numngbs; ++i) {
            const unsigned int index = ngbs[i];
            double r;
            if (!_box.get_sides().x()) {
              r = (p - _positions[index]).norm();
            } else {
              r = _box.periodic_distance(p, _positions[index]).norm();
            }
            const double h = _smoothing_lengths[index];
            const double u = r / h;
            const double m = _masses[index];
            const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
	    _locks[index].lock();
            _neutral_fractions[index] -= splineval/cell_mass * (1. - cell.get_ionization_variables().get_ionic_fraction(ION_H_n));
	    _locks[index].unlock();
      	  }

	}
    }
  }
};

/**
 * @brief Map the state of the grid back to the SPH particle distribution.
 *
 * @param grid DensityGrid to write out.
 * @param iteration Iteration number to use in the snapshot file name(s).
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
/*void SPHArrayInterface::write(DensityGrid &grid, uint_fast32_t iteration,
                              ParameterFile &params, double time) {

  const bool _is_new_density_mapping = false;

  for (unsigned int i =0; i < _neutral_fractions.size(); ++i) {
    _neutral_fractions[i] = 1.0;
  }

  if(_is_new_density_mapping) {
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const CoordinateVector<> position = it.get_cell_midpoint();
      //const std::vector< unsigned int > ngbs = _octree->get_ngbs(p);
      //const unsigned int numngbs = ngbs.size();

      //const Cell &cell=it;
      // Find the vertex that is furthest away from the cell midpoint.
      std::vector< Face > face_vector = it.get_faces();
      double radius = 0.0;
      for (unsigned int i = 0; i < face_vector.size(); i++) {
        for (Face::Vertices j = face_vector[i].first_vertex();
          j != face_vector[i].last_vertex(); ++j) {
          double distance = (j.get_position()-position).norm();
          if (distance > radius)
            radius = distance;
         }
       }

      // Find the neighbours that are contained inside of a sphere of centre the
      // cell midpoint
      // and radius given by the distance to the furthest vertex.
      std::vector< uint_fast32_t > ngbs =
          _octree->get_ngbs_sphere(position, radius);
      const unsigned int numngbs = ngbs.size();

      double cell_mass = 0.;
      std::vector< double > denseval;

      // Loop over all the neighbouring particles and calculate their mass
      // contributions.
      for (unsigned int i = 0; i < numngbs; i++) {
        const unsigned int index = ngbs[i];
        const double h = _smoothing_lengths[index] / 2.0;
        const CoordinateVector<> particle = _positions[index];
        denseval.push_back(mass_contribution(it, particle, h) * _masses[index]);
        cell_mass += denseval[i];
      }

      for (unsigned int i = 0; i < numngbs; i++) {
        const unsigned int index = ngbs[i];
        //const double h = _smoothing_lengths[index] / 2.0;
        //const CoordinateVector<> particle = _positions[index];
        //const double denseval = mass_contribution(it, particle, h) * _masses[index];
        _neutral_fractions[index] -= denseval[i]/cell_mass * (1. - it.get_ionization_variables().get_ionic_fraction(ION_H_n));
      }
    }
  } else {
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const CoordinateVector<> p = it.get_cell_midpoint();
      const std::vector< uint_fast32_t > ngbs = _octree->get_ngbs(p);
      const unsigned int numngbs = ngbs.size();
      double cell_mass = 0.;
      for (unsigned int i = 0; i < numngbs; ++i) {
        const unsigned int index = ngbs[i];
        double r;
        if (!_box.get_sides().x()) {
          r = (p - _positions[index]).norm();
        } else {
          r = _box.periodic_distance(p, _positions[index]).norm();
        }
        const double h = _smoothing_lengths[index];
        const double u = r / h;
        const double m = _masses[index];
        const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
        cell_mass += splineval;
      }

      for (unsigned int i = 0; i < numngbs; ++i) {
        const unsigned int index = ngbs[i];
        double r;
        if (!_box.get_sides().x()) {
          r = (p - _positions[index]).norm();
        } else {
          r = _box.periodic_distance(p, _positions[index]).norm();
        }
        const double h = _smoothing_lengths[index];
        const double u = r / h;
        const double m = _masses[index];
        const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
        _neutral_fractions[index] -= splineval/cell_mass * (1. - it.get_ionization_variables().get_ionic_fraction(ION_H_n));
      }
    }
  }
}*/

void SPHArrayInterface::write(DensityGrid &grid, uint_fast32_t iteration,
                              ParameterFile &params, double time) {
			      
  for (unsigned int i =0; i < _neutral_fractions.size(); ++i) {
    _neutral_fractions[i] = 1.0;
  }

  std::vector< Lock > locks(_neutral_fractions.size());

  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  WorkDistributor< DensityGridTraversalJobMarket< InverseMappingFunction >,
                   DensityGridTraversalJob< InverseMappingFunction > >
      workers;
  //InverseMappingFunction do_calculation(*_octree, _neutral_fractions, locks);
  InverseMappingFunction do_calculation(*_octree,_neutral_fractions,_positions,_smoothing_lengths,_masses,_box,_dens_map,locks);

  DensityGridTraversalJobMarket< InverseMappingFunction > jobs(
      grid, do_calculation, block);
  workers.do_in_parallel(jobs);
}

