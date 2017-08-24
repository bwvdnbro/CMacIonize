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
 * @file IonizationPhotonShootJob.hpp
 *
 * @brief Job implementation that shoots ionizing photons through a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONPHOTONSHOOTJOB_HPP
#define IONIZATIONPHOTONSHOOTJOB_HPP

#include "DensityGrid.hpp"
#include "Photon.hpp"
#include "PhotonSource.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Job implementation that shoots ionizing photons through a DensityGrid.
 */
class IonizationPhotonShootJob {
private:
  /*! @brief PhotonSource that emits photons. */
  const PhotonSource &_photon_source;

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
  inline IonizationPhotonShootJob(const PhotonSource &photon_source,
                                  int random_seed, DensityGrid &density_grid)
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
   * @return False, since the Job is reused and managed by the
   * IonizationPhotonShootJobMarket.
   */
  inline bool do_cleanup() const { return false; }

  /**
   * @brief Shoot _numphoton photons from _photon_source through _density_grid.
   *
   * For each photon, we first randomly generate the photon. We then draw a
   * random optical depth using the distribution function of Steinacker, J.,
   * Baes, M. & Gordon, K. D. 2013, Annu. Rev. Astro. Astrophys., 51, 63
   * (http://adsabs.harvard.edu/abs/2013ARA%26A..51...63S), equation (25). We
   * propagate the photon through the grid until that optical depth is reached,
   * or until the photon leaves the system (adding ionizing luminosity to all
   * cells covered by the path of the photon through the grid). If the photon
   * is still inside the simulation box, we randomly decide whether we need to
   * reemit it or not. If so, we again draw a random optical depth and repeat
   * the whole procedure until the photon is absorbed or leaves the system.
   */
  inline void execute() {
    for (unsigned int i = 0; i < _numphoton; ++i) {
      Photon photon = _photon_source.get_random_photon(_random_generator);
      // if a fraction of light alpha is absorbed when the light traverses a
      // small path with length dl in the material, then the spatial change of
      // the number of photons is given by
      //  dN_p / dl = alpha * N_p
      // we can rewrite this as
      //  dN_p / N_p = alpha * dl
      // which can be integrated to
      //  N_p = N_p0 * exp(-tau)
      // where tau = \int_0^L alpha dl
      // This shows that the optical depth distribution for a photon is
      // exponentially distributed. To draw a random optical depth, we can start
      // from a uniform random number ksi, and invert the cumulative
      // distribution:
      //  tau = -ln(ksi)
      // See Steinacker, Baes & Gordon (2013), equation (25)
      double tau = -std::log(_random_generator.get_uniform_random_double());
      DensityGrid::iterator it = _density_grid.interact(photon, tau);
      while (it != _density_grid.end() &&
             _photon_source.reemit(photon, it.get_ionization_variables(),
                                   _random_generator)) {
        tau = -std::log(_random_generator.get_uniform_random_double());
        it = _density_grid.interact(photon, tau);
      }
      _totweight += photon.get_weight();
      _typecount[photon.get_type()] += photon.get_weight();
    }
  }

  /**
   * @brief Get a name tag for this job.
   *
   * @return "ionizationphotonshootjob".
   */
  inline std::string get_tag() const { return "ionizationphotonshootjob"; }
};

#endif // IONIZATIONPHOTONSHOOTJOB_HPP
