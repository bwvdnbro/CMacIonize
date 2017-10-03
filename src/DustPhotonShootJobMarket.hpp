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

class CCDImage;
class DensityGrid;
class DustScattering;
class PhotonSource;
class RandomGenerator;

/**
 * @brief JobMarket implementation that shoots photons through a dusty
 * DensityGrid.
 */
class DustPhotonShootJobMarket {
private:
  /*! @brief Per thread DustPhotonShootJob. */
  DustPhotonShootJob *_jobs[MAX_NUM_THREADS];

  /*! @brief Number of threads used in the calculation. */
  const int_fast32_t _worksize;

  /*! @brief Total number of photons to propagate through the grid. */
  uint_fast64_t _numphoton;

  /*! @brief Number of photons to shoot during a single DustPhotonShootJob. */
  const uint_fast64_t _jobsize;

  /*! @brief Lock used to ensure safe access to the internal photon number
   *  counters. */
  Lock _lock;

public:
  /**
   * @brief Constructor.
   *
   * @param photon_source PhotonSource that emits photons.
   * @param dust_scattering DustScattering object used to compute scattering off
   * dust.
   * @param random_seed Seed for the RandomGenerator.
   * @param density_grid DensityGrid through which photons are propagated.
   * @param numphoton Total number of photons to propagate through the grid.
   * @param image CCDImage to construct (threads will update a copy of this
   * image).
   * @param jobsize Number of photons to shoot during a single
   * DustPhotonShootJob.
   * @param worksize Number of threads used in the calculation.
   */
  inline DustPhotonShootJobMarket(PhotonSource &photon_source,
                                  const DustScattering &dust_scattering,
                                  int_fast32_t random_seed,
                                  DensityGrid &density_grid,
                                  uint_fast64_t numphoton,
                                  const CCDImage &image, uint_fast64_t jobsize,
                                  int_fast32_t worksize)
      : _worksize(worksize), _numphoton(numphoton), _jobsize(jobsize) {

    // create a separate RandomGenerator for each thread.
    // create a single PhotonShootJob for each thread.
    for (int_fast32_t i = 0; i < _worksize; ++i) {
      _jobs[i] = new DustPhotonShootJob(photon_source, dust_scattering,
                                        random_seed + i, density_grid, image);
    }
  }

  /**
   * @brief Destructor.
   *
   * Deletes the internal job array.
   */
  inline ~DustPhotonShootJobMarket() {
    for (int_fast32_t i = 0; i < _worksize; ++i) {
      delete _jobs[i];
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
   * @brief Set the number of photons.
   *
   * This routine can be used to reset a DustPhotonShootJobMarket that was used
   * before.
   *
   * @param numphoton New number of photons.
   */
  inline void set_numphoton(uint_fast64_t numphoton) { _numphoton = numphoton; }

  /**
   * @brief Update the given CCDImage.
   *
   * @param image CCDImage to update.
   */
  inline void update_image(CCDImage &image) {
    for (int_fast32_t i = 0; i < _worksize; ++i) {
      _jobs[i]->update_image(image);
    }
  }

  /**
   * @brief Get a DustPhotonShootJob.
   *
   * @param thread_id Rank of the thread that wants to get a job (in a parallel
   * context).
   * @return DustPhotonShootJob.
   */
  inline DustPhotonShootJob *get_job(int_fast32_t thread_id) {
    uint_fast64_t jobsize = std::max(_numphoton / (10 * _worksize), _jobsize);
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
