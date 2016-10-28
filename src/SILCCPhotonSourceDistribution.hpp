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

#include "PhotonSourceDistribution.hpp"
#include "Utilities.hpp"

/**
 * @brief PhotonSourceDistribution used for post-procesing SILCC snapshots.
 *
 * We assume a fixed number of sources, uniformly distributed in a rectangular
 * disk in x and y, with a Gaussian distribution in z.
 */
class SILCCPhotonSourceDistribution {
private:
  /*! @brief Number of individual sources. */
  unsigned int _num_sources;

  /*! @brief x component of the anchor of the rectangular disk (in m). */
  double _anchor_x;

  /*! @brief y component of the anchor of the rectangular disk (in m). */
  double _anchor_y;

  /*! @brief x side length of the rectangular disk (in m). */
  double _sides_x;

  /*! @brief y side length of the rectangular disk (in m). */
  double _sides_y;

  /*! @brief Origin of the Gaussian disk height distribution. */
  double _origin_z;

  /*! @brief Scale height of the Gaussian disk height distribution. */
  double _scaleheight_z;

public:
  /**
   * @brief Constructor.
   *
   * @param num_sources Number of individual sources.
   * @param anchor_x x component of the anchor of the rectangular disk (in m).
   * @param sides_x x side length of the rectangular disk (in m).
   * @param anchor_y  y component of the anchor of the rectangular disk (in m).
   * @param sides_y y side length of the rectangular disk (in m).
   * @param origin_z Origin of the Gaussian disk height distribution.
   * @param scaleheight_z Scale height of the Gaussian disk height distribution.
   */
  SILCCPhotonSourceDistribution(unsigned int num_sources, double anchor_x,
                                double sides_x, double anchor_y, double sides_y,
                                double origin_z, double scaleheight_z)
      : _num_sources(num_sources), _anchor_x(anchor_x), _anchor_y(anchor_y),
        _sides_x(sides_x), _sides_y(sides_y), _origin_z(origin_z),
        _scaleheight_z(scaleheight_z) {}

  virtual ~SILCCPhotonSourceDistribution() {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * @return Number of sources.
   */
  virtual unsigned int get_number_of_sources() { return _num_sources; }

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(unsigned int index) {
    double x = _anchor_x + Utilities::random_double() * _sides_x;
    double y = _anchor_y + Utilities::random_double() * _sides_y;
    // we use the Box-Muller method to sample the Gaussian
    double z = _scaleheight_z *
                   std::sqrt(-2. * std::log(Utilities::random_double())) *
                   std::cos(2. * M_PI * Utilities::random_double()) +
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
  virtual double get_weight(unsigned int index) { return 1. / _num_sources; }
};

#endif // SILCCPHOTONSOURCEDISTRIBUTION_HPP
