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

class Abundances;
class CrossSections;
class DensityValues;
class IsotropicContinuousPhotonSource;
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

  /*! @brief Spectrum of the discrete photon sources. */
  PhotonSourceSpectrum *_discrete_spectrum;

  ///

  /// continuous sources

  /*! @brief Total number of photons emitted by the continuous source. */
  unsigned int _continuous_number_of_photons;

  /*! @brief Number of photons already emitted by the continuous source. */
  unsigned int _continuous_active_number_of_photons;

  /*! @brief IsotropicContinuousPhotonSource instance used. */
  IsotropicContinuousPhotonSource *_continuous_source;

  /*! @brief Spectrum of the continuous sources. */
  PhotonSourceSpectrum *_continuous_spectrum;

  ///

  /*! @brief Fraction of photons that is emitted by discrete sources. */
  double _discrete_fraction;

  /*! @brief Total luminosity of all sources (discrete + continuous) (in s^-1).
   */
  double _total_luminosity;

  /*! @brief Abundances of the elements in the ISM. */
  Abundances &_abundances;

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

  /**
   * @brief Update the internal counters that distribute the photon sampling
   * over the various discrete and continuous sources.
   */
  inline void update_indices() {
    if (_continuous_active_number_of_photons < _continuous_number_of_photons) {
      ++_continuous_number_of_photons;
      if (_continuous_active_number_of_photons ==
          _continuous_number_of_photons) {
        // we did all continuous sources and completed the cycle, reset the
        // counters and start again with the discrete sources
        _discrete_active_source_index = 0;
        _discrete_active_photon_index = 0;
        _discrete_active_number_of_photons =
            std::round(_discrete_number_of_photons * _discrete_weights[0]);
        // in case there are no discrete sources
        if (_discrete_active_source_index == _discrete_positions.size()) {
          _continuous_active_number_of_photons = 0;
        }
      }
    } else {
      ++_discrete_active_photon_index;
      if (_discrete_active_photon_index == _discrete_active_number_of_photons) {
        _discrete_active_photon_index = 0;
        ++_discrete_active_source_index;
        if (_discrete_active_source_index < _discrete_positions.size()) {
          _discrete_active_number_of_photons =
              std::round(_discrete_number_of_photons *
                         _discrete_weights[_discrete_active_source_index]);
        } else {
          // we did all discrete sources, now do the continuous source
          _continuous_active_number_of_photons = 0;
          // in case there are no continuous sources
          if (_continuous_active_number_of_photons ==
              _continuous_number_of_photons) {
            _discrete_active_source_index = 0;
            _discrete_active_photon_index = 0;
            _discrete_active_number_of_photons =
                std::round(_discrete_number_of_photons * _discrete_weights[0]);
          }
        }
      }
    }
  }

  void set_cross_sections(Photon &photon, double energy);

public:
  PhotonSource(PhotonSourceDistribution *distribution,
               PhotonSourceSpectrum *discrete_spectrum,
               IsotropicContinuousPhotonSource *continuous_source,
               PhotonSourceSpectrum *continuous_spectrum,
               Abundances &abundances, CrossSections &cross_sections,
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
