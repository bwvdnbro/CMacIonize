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
 * @file SILCCPhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution used for post-processing SILCC snapshots.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SILCCPHOTONSOURCEDISTRIBUTION_HPP
#define SILCCPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"
#include "Utilities.hpp"

/**
 * @brief PhotonSourceDistribution used for post-procesing SILCC snapshots.
 *
 * We assume a fixed number of sources, uniformly distributed in a rectangular
 * disk in x and y, with a Gaussian distribution in z.
 */
class SILCCPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Number of individual sources. */
  const photonsourcenumber_t _num_sources;

  /*! @brief x component of the anchor of the rectangular disk (in m). */
  const double _anchor_x;

  /*! @brief y component of the anchor of the rectangular disk (in m). */
  const double _anchor_y;

  /*! @brief x side length of the rectangular disk (in m). */
  const double _sides_x;

  /*! @brief y side length of the rectangular disk (in m). */
  const double _sides_y;

  /*! @brief Origin of the Gaussian disk height distribution (in m). */
  const double _origin_z;

  /*! @brief Scale height of the Gaussian disk height distribution (in m). */
  const double _scaleheight_z;

  /*! @brief Luminosity of a single source (in s^-1). */
  const double _luminosity;

  /*! @brief RandomGenerator used to generate random numbers. We use a separate
   *  instance for this PhotonSourceDistribution, so that we can change the
   *  sources completely independently from the rest of the algorithm. */
  RandomGenerator _random_generator;

public:
  /**
   * @brief Constructor.
   *
   * @param num_sources Number of individual sources.
   * @param anchor_x x component of the anchor of the rectangular disk (in m).
   * @param sides_x x side length of the rectangular disk (in m).
   * @param anchor_y  y component of the anchor of the rectangular disk (in m).
   * @param sides_y y side length of the rectangular disk (in m).
   * @param origin_z Origin of the Gaussian disk height distribution (in m).
   * @param scaleheight_z Scale height of the Gaussian disk height distribution
   * (in m).
   * @param luminosity Luminosity of a single source (in s^-1).
   * @param random_seed Seed used for the random generator.
   * @param log Log to write logging info to.
   */
  SILCCPhotonSourceDistribution(
      photonsourcenumber_t num_sources, double anchor_x, double sides_x,
      double anchor_y, double sides_y, double origin_z, double scaleheight_z,
      double luminosity, int_fast32_t random_seed = 42, Log *log = nullptr)
      : _num_sources(num_sources), _anchor_x(anchor_x), _anchor_y(anchor_y),
        _sides_x(sides_x), _sides_y(sides_y), _origin_z(origin_z),
        _scaleheight_z(scaleheight_z), _luminosity(luminosity),
        _random_generator(random_seed) {

    if (log) {
      log->write_status("Constructed ", _num_sources,
                        " SILCC sources within the rectangle with anchor [",
                        _anchor_x, " m, ", _anchor_y, " m] and sides [",
                        _sides_x, " m, ", _sides_y, " m], a vertical origin ",
                        _origin_z, " m and scale height ", _scaleheight_z,
                        " m, and a total luminosity of ",
                        _num_sources * _luminosity, " s^-1.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - number of sources: Number of individual sources (default: 24)
   *  - anchor x: X position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides x: X side length of the 2D disc (default: 2. kpc)
   *  - anchor y: Y position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides y: Y side length of the 2D disc (default: 2. kpc)
   *  - origin z: Origin of the exponential disc profile in the z direction
   *    (default: 0. pc)
   *  - scaleheight z: Vertical scale height of the exponential disc profile
   *    (default: 63. pc)
   *  - luminosity: Luminosity of an individual source (default: 3.125e49 s^-1)
   *  - random seed: Random seed used to initialize the random generator that
   *    is used to sample the individual positions (default: 42)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SILCCPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : SILCCPhotonSourceDistribution(
            params.get_value< photonsourcenumber_t >(
                "PhotonSourceDistribution:number of sources", 24),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:anchor x", "-1. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:sides x", "2. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:anchor y", "-1. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:sides y", "2. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:origin z", "0. pc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:scaleheight z", "63. pc"),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "PhotonSourceDistribution:luminosity", "3.125e49 s^-1"),
            params.get_value< int_fast32_t >(
                "PhotonSourceDistribution:random seed", 42),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~SILCCPhotonSourceDistribution() {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_sources() const {
    return _num_sources;
  }

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {

    const double x =
        _anchor_x + _random_generator.get_uniform_random_double() * _sides_x;
    const double y =
        _anchor_y + _random_generator.get_uniform_random_double() * _sides_y;
    // we use the Box-Muller method to sample the Gaussian
    const double z =
        _scaleheight_z *
            std::sqrt(-2. *
                      std::log(_random_generator.get_uniform_random_double())) *
            std::cos(2. * M_PI *
                     _random_generator.get_uniform_random_double()) +
        _origin_z;
    return CoordinateVector<>(x, y, z);
  }

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Reciprocal of the number of sources, as every source has the same
   * weight.
   */
  virtual double get_weight(photonsourcenumber_t index) const {
    return 1. / _num_sources;
  }

  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */
  virtual double get_total_luminosity() const {
    return _num_sources * _luminosity;
  }
};

#endif // SILCCPHOTONSOURCEDISTRIBUTION_HPP
