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
  * @brief A simplified interface to the Unix system timer
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
  Timer();
  ~Timer() {}

  void reset();
  void start();
  double stop();
  double value();

  double interval();
  void restart();
};

#endif // TIMER_HPP
