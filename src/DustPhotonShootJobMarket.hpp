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
 * @file DustPhotonShootJobMarket.hpp
 *
 * @brief JobMarket implementation that shoots photons through a dusty
 * DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DUSTPHOTONSHOOTJOBMARKET_HPP
#define DUSTPHOTONSHOOTJOBMARKET_HPP

#include "Configuration.hpp"
#include "DustPhotonShootJob.hpp"
#include "Lock.hpp"

class PhotonSource;
class RandomGenerator;
class DensityGrid;

/**
 * @brief JobMarket implementation that shoots photons through a dusty
 * DensityGrid.
 */
class DustPhotonShootJobMarket {
private:
  /*! @brief Per thread DustPhotonShootJob. */
  DustPhotonShootJob *_jobs[MAX_NUM_THREADS];

  /*! @brief Number of threads used in the calculation. */
  int _worksize;

  /*! @brief Total number of photons to propagate through the grid. */
  unsigned int _numphoton;

  /*! @brief Number of photons to shoot during a single DustPhotonShootJob. */
  unsigned int _jobsize;

  /*! @brief Lock used to ensure safe access to the internal photon number
   *  counters. */
  Lock _lock;

public:
  /**
   * @brief Constructor.
   *
   * @param photon_source PhotonSource that emits photons.
   * @param random_seed Seed for the RandomGenerator.
   * @param density_grid DensityGrid through which photons are propagated.
   * @param numphoton Total number of photons to propagate through the grid.
   * @param jobsize Number of photons to shoot during a single
   * DustPhotonShootJob.
   * @param worksize Number of threads used in the calculation.
   */
  inline DustPhotonShootJobMarket(PhotonSource &photon_source, int random_seed,
                                  DensityGrid &density_grid,
                                  unsigned int numphoton, unsigned int jobsize,
                                  int worksize)
      : _worksize(worksize), _numphoton(numphoton), _jobsize(jobsize) {
    // create a separate RandomGenerator for each thread.
    // create a single PhotonShootJob for each thread.
    for (int i = 0; i < _worksize; ++i) {
      _jobs[i] =
          new DustPhotonShootJob(photon_source, random_seed + i, density_grid);
    }
  }

  /**
   * @brief Destructor.
   *
   * Deletes the internal job array.
   */
  inline ~DustPhotonShootJobMarket() {
    for (int i = 0; i < _worksize; ++i) {
      delete _jobs[i];
    }
  }

  /**
   * @brief Set the number of parallel threads that will be used to execute
   * the jobs.
   *
   * @param worksize Number of parallel threads that will be used.
   */
  inline void set_worksize(int worksize) {}

  /**
   * @brief Set the number of photons.
   *
   * This routine can be used to reset a DustPhotonShootJobMarket that was used
   * before.
   *
   * @param numphoton New number of photons.
   */
  inline void set_numphoton(unsigned int numphoton) { _numphoton = numphoton; }

  /**
   * @brief Update the given weight counters.
   *
   * @param totweight Total weight of all photons.
   * @param typecount Total weights per photon type.
   */
  inline void update_counters(double &totweight, double *typecount) {
    for (int i = 0; i < _worksize; ++i) {
      _jobs[i]->update_counters(totweight, typecount);
    }
  }

  /**
   * @brief Get a DustPhotonShootJob.
   *
   * @param thread_id Rank of the thread that wants to get a job (in a parallel
   * context).
   * @return DustPhotonShootJob.
   */
  inline DustPhotonShootJob *get_job(int thread_id) {
    unsigned int jobsize = std::max(_numphoton / (10 * _worksize), _jobsize);
    _lock.lock();
    if (jobsize >= _numphoton) {
      jobsize = _numphoton;
    }
    _numphoton -= jobsize;
    _lock.unlock();
    if (jobsize > 0) {
      _jobs[thread_id]->set_numphoton(jobsize);
      return _jobs[thread_id];
    } else {
      return nullptr;
    }
  }
};

#endif // DUSTPHOTONSHOOTJOBMARKET_HPP
