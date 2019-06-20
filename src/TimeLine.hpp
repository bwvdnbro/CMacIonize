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
 * @file TimeLine.hpp
 *
 * @brief Class that maps physical times to an integer power of two time line.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TIMELINE_HPP
#define TIMELINE_HPP

#include "Error.hpp"

#include <algorithm>
#include <cstdint>

/*! @brief Integer time line size, \f$2^{63}\f$. */
#define TIMELINE_MAX_INTEGER_TIMELINE_SIZE 0x8000000000000000

/**
 * @brief Class that maps physical times to an integer power of two time line.
 */
class TimeLine {
private:
  /*! @brief Minimum integer time step size. */
  uint64_t _minimum_timestep;

  /*! @brief Maximum integer time step size. */
  uint64_t _maximum_timestep;

  /*! @brief Conversion factors from integer time to physical time.
   *
   * If the integer time is given by \f$t_i\f$, the physical time \f$t_p\f$ is
   * \f[
   *   t_p = A t_i + B,
   * \f]
   * with \f$A\f$ and \f$B\f$ the two conversion factors.
   *
   * For an integer time interval \f$\Delta{}t_i\f$, the conversion to the
   * corresponding physical time interval \f$\Delta{}t_p\f$ is
   * \f[
   *   \Delta{}t_p = A \Delta{}t_i.
   * \f]
   */
  double _conversion_factors[2];

  /*! @brief Current integer time. */
  uint64_t _current_time;

  /**
   * @brief Convert an integer time to a physical time.
   *
   * @param integer_time Integer time to convert.
   * @return Physical time (in s).
   */
  inline double to_physical_time(uint64_t integer_time) const {
    return _conversion_factors[0] * integer_time + _conversion_factors[1];
  }

  /**
   * @brief Convert an integer time interval to a physical time interval.
   *
   * @param integer_time_interval Integer time interval to convert.
   * @return Physical time interval (in s).
   */
  inline double
  to_physical_time_interval(uint64_t integer_time_interval) const {
    return _conversion_factors[0] * integer_time_interval;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param start_time Start time of the time line (in s).
   * @param end_time End time of the time line (in s).
   * @param minimum_timestep Minimum size of the time step (in s).
   * @param maximum_timestep Maximum size of the time step (in s).
   */
  inline TimeLine(double start_time, double end_time, double minimum_timestep,
                  double maximum_timestep) {

    // compute the total physical time interval covered by the time line
    const double timeline_interval = end_time - start_time;

    // get conversion factor A by mapping 'timeline_interval' to the total range
    // of 64-bit unsigned integer numbers
    _conversion_factors[0] =
        timeline_interval / TIMELINE_MAX_INTEGER_TIMELINE_SIZE;

    // conversion factor B is simply given by 'start_time'
    _conversion_factors[1] = start_time;

    // convert 'minimum_timestep' to its integer counterpart. Round down to the
    // closest power of two
    _minimum_timestep = TIMELINE_MAX_INTEGER_TIMELINE_SIZE;
    if (minimum_timestep > 0) {
      while (to_physical_time_interval(_minimum_timestep) > minimum_timestep) {
        _minimum_timestep >>= 1;
      }

      // make sure '_minimum_timestep' is larger than zero
      const uint64_t one_uint64 = 1;
      _minimum_timestep = std::max(one_uint64, _minimum_timestep);
    } else {
      _minimum_timestep = 1;
    }

    // convert 'maximum_timestep' to its integer counterpart. If
    // 'maximum_timestep' is zero, set it to the total range of integers. If
    // not, make sure '_maximum_timestep' is at least as large as
    // '_minimum_timestep'.
    _maximum_timestep = TIMELINE_MAX_INTEGER_TIMELINE_SIZE;
    if (maximum_timestep > 0) {
      while (to_physical_time_interval(_maximum_timestep) > maximum_timestep) {
        _maximum_timestep >>= 1;
      }

      _maximum_timestep = std::max(_minimum_timestep, _maximum_timestep);
    }

    _current_time = 0;
  }

  /**
   * @brief Advance the time line forward in time.
   *
   * @param requested_timestep Maximum time step that can be taken (in s). The
   * actual time step will be the closest smaller time step that fits in the
   * integer time line.
   * @param actual_timestep Actual time step that was taken (in s). Is set by
   * this function.
   * @param current_time Current time in the physical time line at the end of
   * the step (in s). Is set by this function.
   * @return True if there are still valid time steps after this step, false if
   * the end of the time line was reached.
   */
  inline bool advance(double requested_timestep, double &actual_timestep,
                      double &current_time) {

    // round 'requested_timestep' down to the closest power of two that is a
    // divisor of the time left on the integer time line (such that we do not
    // overshoot the end time by taking N steps with this size)
    uint64_t integer_timestep = _maximum_timestep;
    while (to_physical_time_interval(integer_timestep) > requested_timestep) {
      integer_timestep >>= 1;
    }
    const uint64_t time_left =
        TIMELINE_MAX_INTEGER_TIMELINE_SIZE - _current_time;
    while (time_left % integer_timestep > 0) {
      integer_timestep >>= 1;
    }

    if (integer_timestep < _minimum_timestep) {
      cmac_error("Time step wants to be smaller than minimum time step: %g "
                 "(minimum: %g)!",
                 to_physical_time_interval(integer_timestep),
                 to_physical_time_interval(_minimum_timestep));
    }

    // advance '_current_time'
    _current_time += integer_timestep;

    // set the physical values
    actual_timestep = to_physical_time_interval(integer_timestep);
    current_time = to_physical_time(_current_time);

    return _current_time < TIMELINE_MAX_INTEGER_TIMELINE_SIZE;
  }
};

#endif // TIMELINE_HPP
