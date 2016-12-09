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
#include "Job.hpp"
#include "Photon.hpp"
#include "PhotonSource.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Job implementation that shoots photons through a DensityGrid.
 */
class PhotonShootJob : public Job {
private:
  /*! @brief PhotonSource that emits photons. */
  PhotonSource &_photon_source;

  /*! @brief RandomGenerator used to generate random uniform numbers. */
  RandomGenerator &_random_generator;

  /*! @brief DensityGrid through which photons are propagated. */
  DensityGrid &_density_grid;

  /*! @brief Number of photons to propagate through the DensityGrid. */
  unsigned int _numphoton;

public:
  /**
   * @brief Constructor.
   *
   * @param photon_source PhotonSource that emits photons.
   * @param random_generator RandomGenerator used to generate random uniform
   * numbers.
   * @param density_grid DensityGrid through which photons are propagated.
   * @param numphoton Number of photons to propagate through the DensityGrid.
   */
  PhotonShootJob(PhotonSource &photon_source, RandomGenerator &random_generator,
                 DensityGrid &density_grid, unsigned int numphoton)
      : _photon_source(photon_source), _random_generator(random_generator),
        _density_grid(density_grid), _numphoton(numphoton) {}

  /**
   * @brief Shoot _numphoton photons from _photon_source through _density_grid.
   */
  virtual void execute() {
    for (unsigned int i = 0; i < _numphoton; ++i) {
      Photon photon = _photon_source.get_random_photon();
      double tau = -std::log(_random_generator.get_uniform_random_double());
      while (_density_grid.interact(photon, tau)) {
        unsigned long new_index =
            _density_grid.get_cell_index(photon.get_position());
        if (!_photon_source.reemit(photon,
                                   _density_grid.get_cell_values(new_index))) {
          break;
        }
        tau = -std::log(_random_generator.get_uniform_random_double());
      }
      _photon_source.decommission_photon(photon);
    }
  }
};

#endif // PHOTONSHOOTJOB_HPP
