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
 * @file DustPhotonShootJob.hpp
 *
 * @brief Job implementation that shoots photons through a dusty DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DUSTPHOTONSHOOTJOB_HPP
#define DUSTPHOTONSHOOTJOB_HPP

#include "CCDImage.hpp"
#include "DensityGrid.hpp"
#include "Photon.hpp"
#include "PhotonSource.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Job implementation that shoots photons through a dusty DensityGrid.
 */
class DustPhotonShootJob {
private:
  /*! @brief PhotonSource that emits photons. */
  PhotonSource &_photon_source;

  /*! @brief RandomGenerator used to generate random uniform numbers. */
  RandomGenerator _random_generator;

  /*! @brief DensityGrid through which photons are propagated. */
  DensityGrid &_density_grid;

  /*! @brief CCDImage computed by this thread. */
  CCDImage _image;

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
   * @param image CCDImage to construct.
   */
  inline DustPhotonShootJob(PhotonSource &photon_source, int random_seed,
                            DensityGrid &density_grid, const CCDImage &image)
      : _photon_source(photon_source), _random_generator(random_seed),
        _density_grid(density_grid), _image(image), _numphoton(0) {}

  /**
   * @brief Set the number of photons for the next execution of the job.
   *
   * @param numphoton New number of photons.
   */
  inline void set_numphoton(unsigned int numphoton) { _numphoton = numphoton; }

  /**
   * @brief Update the given CCDImage.
   *
   * @param image CCDImage to update.
   */
  inline void update_image(CCDImage &image) {
    image += _image;
    _image.reset();
  }

  /**
   * @brief Should the Job be deleted by the Worker when it is finished?
   *
   * @return False, since the Job is reused and managed by the
   * DustPhotonShootJobMarket.
   */
  inline bool do_cleanup() const { return false; }

  /**
   * @brief Shoot _numphoton photons from _photon_source through _density_grid.
   */
  inline void execute() {
    // parameter
    const double band_albedo = 0.54; // or 0.21 for band == 1

    for (unsigned int i = 0; i < _numphoton; ++i) {
      Photon photon = _photon_source.get_random_photon(_random_generator);
      // overwrite direction: we need the direction components to speed things
      // up in other parts of the algorithm
      double cost = 2. * _random_generator.get_uniform_random_double() - 1.;
      double sint = std::sqrt(std::max(1. - cost * cost, 0.));
      double phi = 2. * M_PI * _random_generator.get_uniform_random_double();
      double cosp = std::cos(phi);
      double sinp = std::sin(phi);
      const CoordinateVector<> direction(sint * cosp, sint * sinp, cost);
      photon.set_direction(direction);
      photon.set_direction_parameters(sint, cost, phi, sinp, cosp);
      // overwrite cross section
      photon.set_cross_section(ION_H_n, 1.);

      Photon old_photon(photon);
      old_photon.set_direction(
          _image.get_direction(sint, cost, phi, sinp, cosp));
      const double tau_old = _density_grid.integrate_optical_depth(old_photon);
      _image.add_photon(old_photon.get_position(),
                        0.25 * std::exp(-tau_old) / M_PI, 0., 0.);

      double albedo = 1.;
      // make sure the photon scatters at least once by forcing a first
      // interaction
      const double tau_max = _density_grid.integrate_optical_depth(photon);
      const double weight = (1. - std::exp(-tau_max));
      double tau = -std::log(
          1. - _random_generator.get_uniform_random_double() * weight);
      DensityGrid::iterator it = _density_grid.interact(photon, tau);
      while (it != _density_grid.end()) {

        // peel off a photon towards the observer
        Photon new_photon(photon);
        const CoordinateVector<> direction_new =
            _image.get_direction(sint, cost, phi, sinp, cosp);
        const double hgfac = _photon_source.scatter_towards(
            new_photon, direction_new, sint, cost, phi, sinp, cosp);
        const double tau_new =
            _density_grid.integrate_optical_depth(new_photon);
        double fi, fq, fu, fv;
        new_photon.get_stokes_parameters(fi, fq, fu, fv);
        // after every scattering event, the accumulated albedo is reduced
        albedo *= band_albedo;
        const double weight_new = weight * hgfac * albedo * std::exp(-tau_new);
        _image.add_photon(new_photon.get_position(), weight_new * fi,
                          weight_new * fq, weight_new * fu);

        _photon_source.scatter(photon, _random_generator);
        tau = -std::log(_random_generator.get_uniform_random_double());
        it = _density_grid.interact(photon, tau);
      }
    }
  }

  /**
   * @brief Get a name tag for this job.
   *
   * @return "dustphotonshootjob".
   */
  inline std::string get_tag() const { return "dustphotonshootjob"; }
};

#endif // DUSTPHOTONSHOOTJOB_HPP
