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
 * @file TimingTools.hpp
 *
 * @brief Convenient macros for statistical timing information.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TIMINGTOOLS_HPP
#define TIMINGTOOLS_HPP

#include "Error.hpp"
#include "Timer.hpp"

#include <cmath>
#include <vector>

/*! @brief Number of samples to use for statistical timing. */
#define TIMINGTOOLS_NUM_SAMPLE 10

/**
 * @brief Start a timing block with the given name.
 *
 * The block needs to be enclosed withing parentheses (like a for loop), and
 * needs to be ended with a corresponding call to timingtools_end_timing_block.
 *
 * Different timing blocks should not be nested, and a block should contain
 * at least one call to timingtools_start_timing and timingtools_stop_timing (in
 * that order).
 *
 * @param name Name of the timing block.
 */
#define timingtools_start_timing_block(name)                                   \
  {                                                                            \
    cmac_status("Starting timing %s...", name);                                \
    std::vector< double > timingtools_times_array(TIMINGTOOLS_NUM_SAMPLE, 0.); \
    for (unsigned char timingtools_index = 0;                                  \
         timingtools_index < TIMINGTOOLS_NUM_SAMPLE; ++timingtools_index)

/**
 * @brief End the timing block with the given name.
 *
 * See timingtools_start_timing_block for more information.
 *
 * @param name Name of the timing block.
 */
#define timingtools_end_timing_block(name)                                     \
  double timingtools_average_time = 0.;                                        \
  for (unsigned char timingtools_index = 0;                                    \
       timingtools_index < TIMINGTOOLS_NUM_SAMPLE; ++timingtools_index) {      \
    timingtools_average_time += timingtools_times_array[timingtools_index];    \
  }                                                                            \
  timingtools_average_time /= TIMINGTOOLS_NUM_SAMPLE;                          \
  double timingtools_standard_deviation = 0.;                                  \
  for (unsigned char timingtools_index = 0;                                    \
       timingtools_index < TIMINGTOOLS_NUM_SAMPLE; ++timingtools_index) {      \
    const double timingtools_time_diff =                                       \
        timingtools_times_array[timingtools_index] - timingtools_average_time; \
    timingtools_standard_deviation +=                                          \
        timingtools_time_diff * timingtools_time_diff;                         \
  }                                                                            \
  timingtools_standard_deviation /= TIMINGTOOLS_NUM_SAMPLE;                    \
  timingtools_standard_deviation = std::sqrt(timingtools_standard_deviation);  \
  cmac_status("Finished timing %s: %g +- %g s.", name,                         \
              timingtools_average_time, timingtools_standard_deviation);       \
  }

/**
 * @brief Start timing all code between this call and the consecutive call to
 * timingtools_stop_timing.
 */
#define timingtools_start_timing()                                             \
  Timer timingtools_timer;                                                     \
  timingtools_timer.start();

/**
 * @brief Stop timing the code (timing was started when timingtools_start_timing
 * was called).
 */
#define timingtools_stop_timing()                                              \
  timingtools_timer.stop();                                                    \
  timingtools_times_array[timingtools_index] += timingtools_timer.value();

#endif // TIMINGTOOLS_HPP
