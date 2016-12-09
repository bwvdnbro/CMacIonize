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
 * @file PhotonShootJobMarket.hpp
 *
 * @brief JobMarket implementation that shoots photons through a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSHOOTJOBMARKET_HPP
#define PHOTONSHOOTJOBMARKET_HPP

#include "JobMarket.hpp"
#include "Lock.hpp"
#include "PhotonShootJob.hpp"

class PhotonSource;
class RandomGenerator;
class DensityGrid;

/**
 * @brief JobMarket implementation that shoots photons through a DensityGrid.
 */
class PhotonShootJobMarket : public JobMarket {
private:
  /*! @brief PhotonSource that emits photons. */
  PhotonSource &_photon_source;

  /*! @brief RandomGenerator used to generate random uniform numbers. */
  RandomGenerator &_random_generator;

  /*! @brief DensityGrid through which photons are propagated. */
  DensityGrid &_density_grid;

  /*! @brief Total number of photons to propagate through the grid. */
  unsigned int _numphoton;

  /*! @brief Number of photons to shoot during a single PhotonShootJob. */
  unsigned int _jobsize;

  /*! @brief Lock used to ensure safe access to the internal photon number
   *  counters. */
  Lock _lock;

public:
  /**
   * @brief Constructor.
   *
   * @param photon_source PhotonSource that emits photons.
   * @param random_generator RandomGenerator used to generate random uniform
   * numbers.
   * @param density_grid DensityGrid through which photons are propagated.
   * @param numphoton Total number of photons to propagate through the grid.
   * @param jobsize Number of photons to shoot during a single PhotonShootJob.
   */
  PhotonShootJobMarket(PhotonSource &photon_source,
                       RandomGenerator &random_generator,
                       DensityGrid &density_grid, unsigned int numphoton,
                       unsigned int jobsize)
      : _photon_source(photon_source), _random_generator(random_generator),
        _density_grid(density_grid), _numphoton(numphoton), _jobsize(jobsize) {}

  /**
   * @brief Get a PhotonShootJob.
   *
   * @return PhotonShootJob.
   */
  virtual Job *get_job() {
    unsigned int jobsize = _jobsize;
    _lock.lock();
    if (jobsize >= _numphoton) {
      jobsize = _numphoton;
    }
    _numphoton -= jobsize;
    _lock.unlock();
    if (jobsize > 0) {
      return new PhotonShootJob(_photon_source, _random_generator,
                                _density_grid, jobsize);
    } else {
      return nullptr;
    }
  }
};

#endif // PHOTONSHOOTJOBMARKET_HPP
