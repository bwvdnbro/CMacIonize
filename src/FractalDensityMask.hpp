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

#include "Atomic.hpp"
#include "Box.hpp"
#include "DensityGrid.hpp"
#include "DensityMask.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "WorkDistributor.hpp"

#include <iostream>
#include <vector>

/**
 * @brief Fractal density mask that redistributes the density in the given grid
 * according to a fractal distribution.
 *
 * The algorithm used to generate the fractal density field is based on the
 * algorithm described in Elmegreen, B., 1997, ApJ, 477, 196.
 */
class FractalDensityMask : public DensityMask {
private:
  /*! @brief Box containing the fractal distribution (in m). */
  const Box<> _box;

  /*! @brief Resolution of the grid containing the distribution. */
  const CoordinateVector< uint_fast32_t > _resolution;

  /*! @brief Number of particles per level. */
  const uint_fast32_t _N;

  /*! @brief Fractal length scale. */
  const double _L;

  /*! @brief Total number of levels. */
  const uint_fast8_t _num_level;

  /*! @brief First level random number seeds (to guarantee the same fractal
   *  structure for a given seed, independent of the number of threads used to
   *  construct it). */
  std::vector< int_fast32_t > _first_level_seeds;

  /*! @brief Grid containing the distribution. */
  std::vector< std::vector< std::vector< uint_least64_t > > > _distribution;

  /*! @brief Fractal fraction: maximal fraction of the number density in a cell
   *  that is affected by the mask. */
  const double _fractal_fraction;

  /**
   * @brief (Recursively) construct a fractal grid with the given number of
   * levels, and given number of points per level and fractal length scale.
   *
   * @param random_generator RandomGenerator used to generate random positions.
   * @param current_position Offset position for points on the current level
   * (in internal fractal units).
   * @param current_level Current level.
   */
  void make_fractal_grid(
      RandomGenerator &random_generator,
      CoordinateVector<> current_position = CoordinateVector<>(0.),
      uint_fast8_t current_level = 1) {

    CoordinateVector<> x_level = current_position;
    x_level[0] += 2. * (random_generator.get_uniform_random_double() - 0.5) /
                  std::pow(_L, current_level);
    x_level[1] += 2. * (random_generator.get_uniform_random_double() - 0.5) /
                  std::pow(_L, current_level);
    x_level[2] += 2. * (random_generator.get_uniform_random_double() - 0.5) /
                  std::pow(_L, current_level);

    if (current_level < _num_level) {
      for (uint_fast32_t i = 0; i < _N; ++i) {
        make_fractal_grid(random_generator, x_level, current_level + 1);
      }
    } else {
      // current_position now contains coordinates in the range [-1/L, 1/L]
      // map them to the range [0., 1.]
      x_level[0] *= 0.5 * _L;
      x_level[1] *= 0.5 * _L;
      x_level[2] *= 0.5 * _L;
      x_level[0] += 0.5;
      x_level[1] += 0.5;
      x_level[2] += 0.5;
      // some coordinates may have fallen outside the range, periodically map
      // them to coordinates inside
      if (x_level.x() < 0.) {
        x_level[0] += 1.;
      }
      if (x_level.x() >= 1.) {
        x_level[0] -= 1.;
      }
      if (x_level.y() < 0.) {
        x_level[1] += 1.;
      }
      if (x_level.y() >= 1.) {
        x_level[1] -= 1.;
      }
      if (x_level.z() < 0.) {
        x_level[2] += 1.;
      }
      if (x_level.z() >= 1.) {
        x_level[2] -= 1.;
      }

      // make sure the coordinates are inside the range
      cmac_assert(x_level.x() >= 0. && x_level.x() < 1.);
      cmac_assert(x_level.y() >= 0. && x_level.y() < 1.);
      cmac_assert(x_level.z() >= 0. && x_level.z() < 1.);

      // map the coordinates to grid indices and add a point to the
      // corresponding cell
      const uint_fast32_t ix = x_level.x() * _distribution.size();
      const uint_fast32_t iy = x_level.y() * _distribution[ix].size();
      const uint_fast32_t iz = x_level.z() * _distribution[ix][iy].size();

      // use an atomic operation to add the point, to make this method thread
      // safe
      const uint_least64_t add = 1;
      Atomic::add(_distribution[ix][iy][iz], add);
      cmac_assert(_distribution[ix][iy][iz] < 0xffffffffffffffff);
    }
  }

  /**
   * @brief Job that constructs part of the fractal density grid.
   */
  class FractalDensityMaskConstructionJob {
  private:
    /*! @brief Reference to the FractalDensityMask on which we act. */
    FractalDensityMask &_mask;

    /*! @brief RandomGenerator used by this job. */
    RandomGenerator _random_generator;

  public:
    /**
     * @brief Constructor.
     *
     * @param mask Reference to the FractalDensityMask on which we act.
     */
    inline FractalDensityMaskConstructionJob(FractalDensityMask &mask)
        : _mask(mask) {}

    /**
     * @brief Set the index for the next first level block that will be
     * constructed by this job.
     *
     * This sets the seed of the random generator to the appropriate value for
     * this first level block.
     *
     * @param index Index of the next first level block that will be constructed
     * by this job.
     */
    inline void set_index(uint_fast32_t index) {
      _random_generator.set_seed(_mask._first_level_seeds[index]);
    }

    /**
     * @brief Should the Job be deleted by the Worker when it is finished?
     *
     * @return False, since we want to keep using the same RandomGenerator.
     */
    inline bool do_cleanup() const { return false; }

    /**
     * @brief Construct the entire fractal hierarchy for _index on the first
     * level.
     */
    inline void execute() { _mask.make_fractal_grid(_random_generator); }

    /**
     * @brief Get a name tag for this job.
     *
     * @return "fractaldensitymask_construction".
     */
    inline std::string get_tag() { return "fractaldensitymask_construction"; }
  };

  /**
   * @brief JobMarket for FractalDensityMaskConstructionJobs.
   */
  class FractalDensityMaskConstructionJobMarket {
  private:
    /*! @brief Per thread FractalDensityMaskConstructionJob. */
    std::vector< FractalDensityMaskConstructionJob * > _jobs;

    /*! @brief Current index on the first level. */
    uint_fast32_t _current_index;

    /*! @brief Total number of particles on the first level. */
    const uint_fast32_t _num_part;

    /*! @brief Lock used to ensure safe access to the internal index counter. */
    Lock _lock;

  public:
    /**
     * @brief Constructor.
     *
     * @param mask FractalDensityMask to operate on.
     * @param worksize Number of threads to use.
     */
    inline FractalDensityMaskConstructionJobMarket(FractalDensityMask &mask,
                                                   int_fast32_t worksize)
        : _current_index(0), _num_part(mask._N) {

      _jobs.reserve(worksize);
      for (int_fast32_t i = 0; i < worksize; ++i) {
        _jobs.push_back(new FractalDensityMaskConstructionJob(mask));
      }
    }

    /**
     * @brief Set the number of parallel threads that will be used to execute
     * the jobs.
     *
     * @param worksize Number of parallel threads that will be used.
     */
    inline void set_worksize(int_fast32_t worksize) {}

    /**
     * @brief Destructor.
     *
     * Clean up memory used by jobs.
     */
    inline ~FractalDensityMaskConstructionJobMarket() {
      for (size_t i = 0; i < _jobs.size(); ++i) {
        delete _jobs[i];
      }
    }

    /**
     * @brief Get a FractalDensityMaskConstructionJob.
     *
     * @param thread_id Rank of the thread that wants to get a job (in a
     * parallel context).
     * @return Pointer to a FractalDensityMaskConstructionJob, or a null pointer
     * if no more jobs are available.
     */
    inline FractalDensityMaskConstructionJob *get_job(int_fast32_t thread_id) {

      if (_current_index == _num_part) {
        return nullptr;
      }
      bool has_job = false;
      uint_fast32_t current_index = 0;
      // we can only change _current_index if the lock is locked
      // similarly, we can only take decision based on the value of
      // _current_index while the lock is locked (because otherwise another
      // thread might invalidate the decision we take, while we take it)
      _lock.lock();
      ++_current_index;
      if (_current_index < _num_part) {
        has_job = true;
        current_index = _current_index;
      }
      _lock.unlock();
      if (has_job) {
        _jobs[thread_id]->set_index(current_index);
        return _jobs[thread_id];
      } else {
        return nullptr;
      }
    }
  };

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
   * @param fractal_fraction Fractal fraction: maximal fraction of the number
   * density in each cell that is affected by the mask.
   * @param log Log to write logging info to.
   */
  FractalDensityMask(const Box<> &box,
                     const CoordinateVector< uint_fast32_t > &resolution,
                     uint_fast32_t numpart, int_fast32_t seed,
                     double fractal_dimension, uint_fast8_t num_level,
                     double fractal_fraction, Log *log = nullptr)
      : _box(box), _resolution(resolution),
        // we will use equal values for the number of points (N) per level
        _N(std::ceil(std::pow(numpart, 1. / num_level))),
        // the fractal length scale (L) is linked to the fractal dimension (D)
        // and the number of points per level by the formula
        //  D = log10(N)/log10(L)
        _L(std::pow(10., std::log10(_N) / fractal_dimension)),
        _num_level(num_level), _fractal_fraction(fractal_fraction) {

    // allocate the grid
    _distribution.resize(resolution.x());
    for (uint_fast32_t ix = 0; ix < _resolution.x(); ++ix) {
      _distribution[ix].resize(resolution.y());
      for (uint_fast32_t iy = 0; iy < _resolution.y(); ++iy) {
        _distribution[ix][iy].resize(_resolution.z(), 0);
      }
    }

    // set the seeds for the random generator on the first level
    // we do this to guarantee the same fractal structure, independent of which
    // thread is used to construct which block, and in which order
    _first_level_seeds.resize(_N, 0);
    RandomGenerator random_generator(seed);
    for (uint_fast32_t i = 0; i < _N; ++i) {
      _first_level_seeds[i] = random_generator.get_random_integer();
    }

    if (log) {
      log->write_status(
          "Created FractalDensityMask with ", _num_level, " levels, ", _N,
          " points per level, and a fractal length scale of ", _L, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - box anchor: Anchor of the region where the mask is applied (default:
   *    [-5. pc, -5. pc, -5. pc])
   *  - box sides: Side lengths of the region where the mask is applied
   *    (default: [10. pc, 10. pc, 10. pc])
   *  - resolution: Resolution of the grid on which the fractal distribution is
   *    computed (default: [20, 20, 20])
   *  - number of particles: Number of particles used to sample the fractal
   *    distribution (default: 1e6)
   *  - random seed: Random seed used for the random number generator (default:
   *    42)
   *  - fractal dimension: Fractal dimension of the distribution (default: 2.6)
   *  - number of levels: Number of levels to sample (default: 4)
   *  - fractal fraction: Fraction of the density distribution that should be
   *    weighted according to the fractal distribution (default: 1.)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  FractalDensityMask(ParameterFile &params, Log *log = nullptr)
      : FractalDensityMask(
            Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                      "DensityMask:box anchor", "[-5. pc, -5. pc, -5. pc]"),
                  params.get_physical_vector< QUANTITY_LENGTH >(
                      "DensityMask:box sides", "[10. pc, 10. pc, 10. pc]")),
            params.get_value< CoordinateVector< uint_fast32_t > >(
                "DensityMask:resolution",
                CoordinateVector< uint_fast32_t >(20)),
            params.get_value< uint_fast32_t >("DensityMask:number of particles",
                                              1e6),
            params.get_value< int_fast32_t >("DensityMask:random seed", 42),
            params.get_value< double >("DensityMask:fractal dimension", 2.6),
            params.get_value< uint_fast8_t >("DensityMask:number of levels", 4),
            params.get_value< double >("DensityMask:fractal fraction", 1.),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~FractalDensityMask() {}

  /**
   * @brief Initialize the mask (in parallel) using the given number of threads.
   *
   * @param worksize Number of threads to use.
   */
  virtual void initialize(int_fast32_t worksize = -1) {

    WorkDistributor< FractalDensityMaskConstructionJobMarket,
                     FractalDensityMaskConstructionJob >
        workers(worksize);
    worksize = workers.get_worksize();
    FractalDensityMaskConstructionJobMarket jobs(*this, worksize);
    workers.do_in_parallel(jobs);
  }

  /**
   * @brief Apply the mask to the given DensityGrid, using the given fractal
   * fraction.
   *
   * Note that the mask is only applied to those cells that have their centers
   * inside the mask region.
   *
   * Note also that this operation preserves the total number of hydrogen atoms
   * in the grid (up to machine precision). So not material is removed; it is
   * only moved around to generate a fractal distribution.
   *
   * @param grid DensityGrid to apply the mask to.
   */
  virtual void apply(DensityGrid &grid) const {

    const double smooth_fraction = 1. - _fractal_fraction;
    double Ntot = 0.;
    double Nsmooth = 0.;
    double Nfractal = 0.;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      CoordinateVector<> midpoint = it.get_cell_midpoint();
      if (_box.inside(midpoint)) {
        const double ncell = it.get_ionization_variables().get_number_density();
        const double Ncell = ncell * it.get_volume();
        Ntot += Ncell;
        const uint_fast32_t ix = (midpoint.x() - _box.get_anchor().x()) /
                                 _box.get_sides().x() * _resolution.x();
        const uint_fast32_t iy = (midpoint.y() - _box.get_anchor().y()) /
                                 _box.get_sides().y() * _resolution.y();
        const uint_fast32_t iz = (midpoint.z() - _box.get_anchor().z()) /
                                 _box.get_sides().z() * _resolution.z();
        Nsmooth += smooth_fraction * Ncell;
        Nfractal += _fractal_fraction * Ncell * _distribution[ix][iy][iz];
      }
    }

    cmac_assert(Nfractal > 0.);

    const double fractal_norm = (Ntot - Nsmooth) / Nfractal;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const CoordinateVector<> midpoint = it.get_cell_midpoint();
      if (_box.inside(midpoint)) {
        const uint_fast32_t ix = (midpoint.x() - _box.get_anchor().x()) /
                                 _box.get_sides().x() * _resolution.x();
        const uint_fast32_t iy = (midpoint.y() - _box.get_anchor().y()) /
                                 _box.get_sides().y() * _resolution.y();
        const uint_fast32_t iz = (midpoint.z() - _box.get_anchor().z()) /
                                 _box.get_sides().z() * _resolution.z();
        const double ncell = it.get_ionization_variables().get_number_density();
        const double nsmooth = smooth_fraction * ncell;
        const double nfractal = _fractal_fraction * fractal_norm * ncell *
                                _distribution[ix][iy][iz];
        it.get_ionization_variables().set_number_density(nsmooth + nfractal);
      }
    }
  }
};

#endif // FRACTALDENSITYMASK_HPP
