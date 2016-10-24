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
 * @file SingleStarPhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution implementation containing a single stellar
 * source.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SINGLESTARPHOTONSOURCEDISTRIBUTION_HPP
#define SINGLESTARPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"

/**
 * @brief General interface for photon source distribution functors.
 */
class SingleStarPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Position of the single stellar source (in m). */
  CoordinateVector<> _position;

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  /**
   * @brief Constructor.
   *
   * @param position Position of the single stellar source (in m).
   * @param log Log to write logging information to.
   */
  SingleStarPhotonSourceDistribution(CoordinateVector<> position,
                                     Log *log = nullptr)
      : _position(position), _log(log) {
    if (_log) {
      _log->write_status(
          "Created SingleStarPhotonSourceDistribution at position [",
          _position.x(), " m,", _position.y(), " m,", _position.z(), " m].");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   */
  SingleStarPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : SingleStarPhotonSourceDistribution(
            params.get_physical_vector< QUANTITY_LENGTH >(
                "photonsourcedistribution.position", "[0.5 m, 0.5 m, 0.5 m]"),
            log) {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * @return 1, as this distribution contains a single stellar source
   */
  virtual unsigned int get_number_of_sources() { return 1; };

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of the single stellar source position (in m).
   */
  virtual CoordinateVector<> get_position(unsigned int index) {
    return _position;
  };

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the single photon source: 1.
   */
  virtual double get_weight(unsigned int index) { return 1.; }
};

#endif // PHOTONSOURCEDISTRIBUTION_HPP
