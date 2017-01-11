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
class ContinuousPhotonSource;
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

  /*! @brief Positions of the discrete photon sources (in m). */
  std::vector< CoordinateVector<> > _discrete_positions;

  /*! @brief Spectrum of the discrete photon sources. */
  PhotonSourceSpectrum *_discrete_spectrum;

  /*! @brief Weight of discrete photons. */
  double _discrete_photon_weight;

  /*! @brief Probabilities that a photon is emitted by a specific discrete
   *  source. */
  std::vector< double > _discrete_probabilities;

  /// continuous sources

  /*! @brief ContinuousPhotonSource instance used. */
  ContinuousPhotonSource *_continuous_source;

  /*! @brief Spectrum of the continuous sources. */
  PhotonSourceSpectrum *_continuous_spectrum;

  /*! @brief Weight of continuous photons. */
  double _continuous_photon_weight;

  /*! @brief Probability that a photon is emitted by the continuous source. */
  double _continuous_probability;

  /// all sources

  /*! @brief Total luminosity of all sources (discrete + continuous) (in s^-1).
   */
  double _total_luminosity;

  /*! @brief Abundances of the elements in the ISM. */
  Abundances &_abundances;

  /*! @brief Cross sections for photoionization. */
  CrossSections &_cross_sections;

  /*! @brief Hydrogen Lyman continuum spectrum, used for re-emission. */
  HydrogenLymanContinuumSpectrum _HLyc_spectrum;

  /*! @brief Helium Lyman continuum spectrum, used for re-emission. */
  HeliumLymanContinuumSpectrum _HeLyc_spectrum;

  /*! @brief Helium 2-photon continuum spectrum, used for re-emission. */
  HeliumTwoPhotonContinuumSpectrum _He2pc_spectrum;

  /*! @brief Log to write logging info to. */
  Log *_log;

  void set_cross_sections(Photon &photon, double energy) const;

public:
  PhotonSource(PhotonSourceDistribution *distribution,
               PhotonSourceSpectrum *discrete_spectrum,
               ContinuousPhotonSource *continuous_source,
               PhotonSourceSpectrum *continuous_spectrum,
               Abundances &abundances, CrossSections &cross_sections,
               Log *log = nullptr);

  /**
   * @brief Get a random direction.
   *
   * @param random_generator RandomGenerator to use.
   * @return CoordinateVector containing the components of a random isotropic
   * direction.
   */
  inline static CoordinateVector<>
  get_random_direction(RandomGenerator &random_generator) {
    double cost = 2. * random_generator.get_uniform_random_double() - 1.;
    double sint = 1. - cost * cost;
    sint = std::sqrt(std::max(sint, 0.));
    double phi = 2. * M_PI * random_generator.get_uniform_random_double();
    double cosp = std::cos(phi);
    double sinp = std::sin(phi);
    return CoordinateVector<>(sint * cosp, sint * sinp, cost);
  }

  Photon get_random_photon(RandomGenerator &random_generator) const;

  double get_total_luminosity() const;

  bool reemit(Photon &photon, DensityValues &cell,
              RandomGenerator &random_generator) const;
};

#endif // PHOTONSOURCE_HPP
