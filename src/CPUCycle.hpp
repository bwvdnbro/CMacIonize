/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file CPUCycle.hpp
 *
 * @brief Utility macros to access the CPU cycle count.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CPUCYCLE_HPP
#define CPUCYCLE_HPP

/**
 * @brief Get the CPU cycle time stamp.
 *
 * @param time_variable Variable to store the result in.
 */
#define cpucycle_tick(time_variable)                                           \
  {                                                                            \
    unsigned int lo, hi;                                                       \
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));                        \
    time_variable = ((unsigned long)hi << 32) | lo;                            \
  }

#endif // CPUCYCLE_HPP
