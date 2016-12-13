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
 * @file PhotonShootJob.hpp
 *
 * @brief Job implementation that shoots photons through a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSHOOTJOB_HPP
#define PHOTONSHOOTJOB_HPP

#include "DensityGrid.hpp"
#include "Photon.hpp"
#include "PhotonSource.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Job implementation that shoots photons through a DensityGrid.
 */
class PhotonShootJob {
private:
  /*! @brief PhotonSource that emits photons. */
  PhotonSource &_photon_source;

  /*! @brief RandomGenerator used to generate random uniform numbers. */
  RandomGenerator _random_generator;

  /*! @brief DensityGrid through which photons are propagated. */
  DensityGrid &_density_grid;

  /*! @brief Total weight of all photons. */
  double _totweight;

  /*! @brief Total weights per photon type. */
  double _typecount[PHOTONTYPE_NUMBER];

  /*! @brief Number of photons to propagate through the DensityGrid. */
  unsigned int _numphoton;

public:
  /**
   * @brief Constructor.
   *
   * @param photon_source PhotonSource that emits photons.
   * @param random_seed Seed for the RandomGenerator used by this specific
   * thread.
   * @param density_grid DensityGrid through which photons are propagated.
   */
  inline PhotonShootJob(PhotonSource &photon_source, int random_seed,
                        DensityGrid &density_grid)
      : _photon_source(photon_source), _random_generator(random_seed),
        _density_grid(density_grid), _totweight(0.), _typecount{0.},
        _numphoton(0) {}

  /**
   * @brief Set the number of photons for the next execution of the job.
   *
   * @param numphoton New number of photons.
   */
  inline void set_numphoton(unsigned int numphoton) { _numphoton = numphoton; }

  /**
   * @brief Update the given weight counters and reset the internal counters.
   *
   * @param totweight Total weight of all photons.
   * @param typecount Total weights per photon type.
   */
  inline void update_counters(double &totweight, double *typecount) {
    totweight += _totweight;
    _totweight = 0.;
    for (int i = 0; i < PHOTONTYPE_NUMBER; ++i) {
      typecount[i] += _typecount[i];
      _typecount[i] = 0.;
    }
  }

  /**
   * @brief Should the Job be deleted by the Worker when it is finished?
   *
   * @return False, since the Job is reused and managed by PhotonShootJobMarket.
   */
  inline bool do_cleanup() const { return false; }

  /**
   * @brief Shoot _numphoton photons from _photon_source through _density_grid.
   */
  inline void execute() {
    for (unsigned int i = 0; i < _numphoton; ++i) {
      Photon photon = _photon_source.get_random_photon(_random_generator);
      double tau = -std::log(_random_generator.get_uniform_random_double());
      while (_density_grid.interact(photon, tau)) {
        unsigned long new_index =
            _density_grid.get_cell_index(photon.get_position());
        if (!_photon_source.reemit(photon,
                                   _density_grid.get_cell_values(new_index),
                                   _random_generator)) {
          break;
        }
        tau = -std::log(_random_generator.get_uniform_random_double());
      }
      _totweight += photon.get_weight();
      _typecount[photon.get_type()] += photon.get_weight();
    }
  }
};

#endif // PHOTONSHOOTJOB_HPP
