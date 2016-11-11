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
 * @file PhotonSource.hpp
 *
 * @brief Photon source: header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCE_HPP
#define PHOTONSOURCE_HPP

#include "CoordinateVector.hpp"
#include "HeliumLymanContinuumSpectrum.hpp"
#include "HeliumTwoPhotonContinuumSpectrum.hpp"
#include "HydrogenLymanContinuumSpectrum.hpp"
#include "Photon.hpp"
#include "RandomGenerator.hpp"
#include "Utilities.hpp"

#include <cmath>
#include <vector>

class CrossSections;
class DensityValues;
class Log;
class PhotonSourceDistribution;
class PhotonSourceSpectrum;

/**
 * @brief Photon source, that contains the actual sources that are used during
 * the radiative transfer loop.
 *
 * Note that for some PhotonSourceDistribution implementations, the information
 * contained within this class and the PhotonSourceDistribution will be very
 * similar. However, the PhotonSourceDistribution could also be a smooth
 * analytic distribution, in which case the PhotonSource will hold a discrete
 * version of this distribution.
 */
class PhotonSource {
private:
  /// discrete sources

  /*! @brief Total number of photons emitted by all discrete sources. */
  unsigned int _discrete_number_of_photons;
  /*! @brief Number of photons emitted by the currently active source. */
  unsigned int _discrete_active_number_of_photons;
  /*! @brief Currently emitted photon index. */
  unsigned int _discrete_active_photon_index;
  /*! @brief Currently active photon source index. */
  unsigned int _discrete_active_source_index;

  /*! @brief Positions of the discrete photon sources (in m). */
  std::vector< CoordinateVector<> > _discrete_positions;
  /*! @brief Weights of the discrete photon sources. */
  std::vector< double > _discrete_weights;

  /*! @brief Total luminosity of all discrete sources together (in s^-1). */
  double _discrete_total_luminosity;

  /*! @brief Spectrum of the discrete photon sources. */
  PhotonSourceSpectrum &_discrete_spectrum;

  ///

  /// continuous sources

  /*! @brief Total number of photons emitted by the continuous sources. */
  unsigned int _continuous_number_of_photons;

  /*! @brief Number of photons emitted by the continuous sources. */
  unsigned int _continuous_active_number_of_photons;

  /*! @brief Spectrum of the continuous sources. */
  PhotonSourceSpectrum *_continous_spectrum;

  ///

  /*! @brief Cross sections for photoionization. */
  CrossSections &_cross_sections;

  /*! @brief RandomGenerator used to generate random numbers. */
  RandomGenerator &_random_generator;

  /*! @brief Hydrogen Lyman continuum spectrum, used for re-emission. */
  HydrogenLymanContinuumSpectrum _HLyc_spectrum;

  /*! @brief Helium Lyman continuum spectrum, used for re-emission. */
  HeliumLymanContinuumSpectrum _HeLyc_spectrum;

  /*! @brief Helium 2-photon continuum spectrum, used for re-emission. */
  HeliumTwoPhotonContinuumSpectrum _He2pc_spectrum;

  /*! @brief Log to write logging info to. */
  Log *_log;

public:
  PhotonSource(PhotonSourceDistribution &distribution,
               PhotonSourceSpectrum &spectrum, CrossSections &cross_sections,
               RandomGenerator &random_generator, Log *log = nullptr);

  unsigned int set_number_of_photons(unsigned int number_of_photons);

  /**
   * @brief Get a random direction.
   *
   * @return CoordinateVector containing the components of a random isotropic
   * direction.
   */
  inline CoordinateVector<> get_random_direction() {
    double cost = 2. * _random_generator.get_uniform_random_double() - 1.;
    double sint = 1. - cost * cost;
    sint = std::sqrt(std::max(sint, 0.));
    double phi = 2. * M_PI * _random_generator.get_uniform_random_double();
    double cosp = std::cos(phi);
    double sinp = std::sin(phi);
    return CoordinateVector<>(sint * cosp, sint * sinp, cost);
  }

  Photon get_random_photon();

  double get_total_luminosity();

  bool reemit(Photon &photon, DensityValues &cell);
};

#endif // PHOTONSOURCE_HPP
