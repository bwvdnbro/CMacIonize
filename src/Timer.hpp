/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Timer.hpp
 *
 * @brief A simplified interface to the Unix system timer.
 *
 * This file was originally part of the public moving mesh code Shadowfax
 * (https://github.com/AstroUGent/shadowfax). We removed the restart routines,
 * everything else is unchanged.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TIMER_HPP
#define TIMER_HPP

#include <sys/time.h> // for timeval

/**
  * @brief A simplified interface to the Unix system timer.
  *
  * The Timer automatically registers the current system time when constructed
  * and returns the elapsed time in seconds when it is stopped.
  *
  * The Timer can also be used to time multiple intervals, by using the
  * functions Timer::start and Timer::stop. The function Timer::stop always
  * returns the total registered time, which is the sum of all individual
  * intervals measured.
  */
class Timer {
private:
  /*! @brief Starting time of the timer */
  timeval _start;

  /*! @brief Stop time of the timer */
  timeval _stop;

  /*! @brief Total time interval registered so far */
  timeval _diff;

public:
  /**
   * @brief Clear the internal timeval difference.
   */
  inline void reset() { timerclear(&_diff); }

  /**
   * @brief Constructor.
   *
   * Intialize the internal timeval difference and register the current system
   * time.
   */
  inline Timer() {
    reset();
    gettimeofday(&_start, nullptr);
  }

  /**
   * @brief Record the current system time as starting time.
   */
  inline void start() { gettimeofday(&_start, nullptr); }

  /**
   * @brief Record the current system time as stopping time and add the
   * difference between start and stop to the internal timeval difference.
   *
   * @return The current contents of the internal timeval difference in seconds
   * (with microsecond precision).
   */
  inline double stop() {
    gettimeofday(&_stop, NULL);
    timeval interval_diff;
    timersub(&_stop, &_start, &interval_diff);
    timeradd(&_diff, &interval_diff, &_diff);
    return _diff.tv_sec + 1.e-6 * _diff.tv_usec;
  }

  /**
   * @brief Get the current internal timeval difference.
   *
   * @return The current contents of the internal timeval difference in seconds
   * (with microsecond precision).
   */
  inline double value() const { return _diff.tv_sec + 1.e-6 * _diff.tv_usec; }

  /**
   * @brief Get the current value of the timer without affecting it.
   *
   * @return The time in seconds since the timer was last started.
   */
  inline double interval() {
    timeval tempstop;
    gettimeofday(&tempstop, NULL);
    timeval interval;
    timersub(&tempstop, &_start, &interval);
    return interval.tv_sec + 1.e-6 * interval.tv_usec;
  }

  /**
   * @brief Restart the timer by overwriting the start time.
   */
  inline void restart() { gettimeofday(&_start, nullptr); }
};

#endif // TIMER_HPP
