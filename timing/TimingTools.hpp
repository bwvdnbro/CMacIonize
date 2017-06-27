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

#include "CommandLineParser.hpp"
#include "CompilerInfo.hpp"
#include "Timer.hpp"
#include "Utilities.hpp"
#include "WorkEnvironment.hpp"

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

/**
 * @brief Wrapper around printf.
 *
 * @param s String to print out.
 */
#define timingtools_print(s, ...) printf(s "\n", ##__VA_ARGS__)

/**
 * @brief Wrapper around printf that adds a '#' banner around the printed text.
 *
 * @param s String to print out.
 */
#define timingtools_print_header(s, ...)                                       \
  printf("%s\n", std::string(80, '#').c_str());                                \
  printf("# " s "\n", ##__VA_ARGS__);                                          \
  printf("%s\n", std::string(80, '#').c_str());

/**
 * @brief Initialize the timing tools based on the given command line arguments.
 *
 * All other macros in this file only work after this macro has been called.
 *
 * @param name Name of the timing test.
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 */
#define timingtools_init(name, argc, argv)                                     \
  CommandLineParser timingtools_command_line_parser(name);                     \
  timingtools_command_line_parser.add_option(                                  \
      "number_of_samples", 'n',                                                \
      "Set the number of samples to use for statistical timing.",              \
      COMMANDLINEOPTION_INTARGUMENT, "10");                                    \
  timingtools_command_line_parser.add_option(                                  \
      "number_of_threads", 't',                                                \
      "Set the maximum number of threads available on the system.",            \
      COMMANDLINEOPTION_INTARGUMENT, "1");                                     \
  timingtools_command_line_parser.parse_arguments(argc, argv);                 \
  const unsigned int timingtools_num_sample =                                  \
      timingtools_command_line_parser.get_value< int >("number_of_samples");   \
  const unsigned int timingtools_num_threads =                                 \
      timingtools_command_line_parser.get_value< int >("number_of_threads");   \
  (void)timingtools_num_sample;                                                \
  (void)timingtools_num_threads;

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
    timingtools_print("Starting timing %s...", name);                          \
    std::vector< double > timingtools_times_array(timingtools_num_sample, 0.); \
    for (unsigned char timingtools_index = 0;                                  \
         timingtools_index < timingtools_num_sample; ++timingtools_index)

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
       timingtools_index < timingtools_num_sample; ++timingtools_index) {      \
    timingtools_average_time += timingtools_times_array[timingtools_index];    \
  }                                                                            \
  timingtools_average_time /= timingtools_num_sample;                          \
  double timingtools_standard_deviation = 0.;                                  \
  for (unsigned char timingtools_index = 0;                                    \
       timingtools_index < timingtools_num_sample; ++timingtools_index) {      \
    const double timingtools_time_diff =                                       \
        timingtools_times_array[timingtools_index] - timingtools_average_time; \
    timingtools_standard_deviation +=                                          \
        timingtools_time_diff * timingtools_time_diff;                         \
  }                                                                            \
  timingtools_standard_deviation /= timingtools_num_sample;                    \
  timingtools_standard_deviation = std::sqrt(timingtools_standard_deviation);  \
  timingtools_print("Finished timing %s: %g +- %g s.", name,                   \
                    timingtools_average_time, timingtools_standard_deviation); \
  }

/**
 * @brief Start timing all code between this call and the consecutive call to
 * timingtools_stop_timing.
 *
 * This only works within a block started with either
 * timingtools_start_timing_block or timingtools_start_scaling_block.
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

/**
 * @brief Start a scaling block with the given name.
 *
 * The scaling block needs to be enclosed within parentheses, like a for loop.
 *
 * The scaling block needs to be ended with a corresponding call to
 * timingtools_end_scaling_block.
 *
 * @param name Name of the scaling block.
 */
#define timingtools_start_scaling_block(name)                                  \
  {                                                                            \
    timingtools_print("Starting scaling test for %s...", name);                \
    std::vector< double > timingtools_scaling_array(timingtools_num_threads,   \
                                                    0.);                       \
    std::vector< double > timingtools_scaling_standard_deviation(              \
        timingtools_num_threads, 0.);                                          \
    for (unsigned char timingtools_current_num_threads = 0;                    \
         timingtools_current_num_threads < timingtools_num_threads;            \
         ++timingtools_current_num_threads) {                                  \
      WorkEnvironment::set_max_num_threads(timingtools_current_num_threads +   \
                                           1);                                 \
      std::vector< double > timingtools_times_array(timingtools_num_sample,    \
                                                    0.);                       \
      for (unsigned char timingtools_index = 0;                                \
           timingtools_index < timingtools_num_sample; ++timingtools_index)

/**
 * @brief End the scaling block with the given name.
 *
 * See timingtools_start_scaling_block for more information.
 *
 * @param name Name of the scaling block.
 * @param filename Name of the file to write scaling statistics to.
 */
#define timingtools_end_scaling_block(name, filename)                          \
  for (unsigned char timingtools_index = 0;                                    \
       timingtools_index < timingtools_num_sample; ++timingtools_index) {      \
    timingtools_scaling_array[timingtools_current_num_threads] +=              \
        timingtools_times_array[timingtools_index];                            \
  }                                                                            \
  timingtools_scaling_array[timingtools_current_num_threads] /=                \
      timingtools_num_sample;                                                  \
  for (unsigned char timingtools_index = 0;                                    \
       timingtools_index < timingtools_num_sample; ++timingtools_index) {      \
    const double timingtools_time_diff =                                       \
        timingtools_times_array[timingtools_index] -                           \
        timingtools_scaling_array[timingtools_current_num_threads];            \
    timingtools_scaling_standard_deviation[timingtools_current_num_threads] += \
        timingtools_time_diff * timingtools_time_diff;                         \
  }                                                                            \
  timingtools_scaling_standard_deviation[timingtools_current_num_threads] /=   \
      timingtools_num_sample;                                                  \
  }                                                                            \
  timingtools_print("Finished scaling test for %s:", name);                    \
  timingtools_print("number of threads\ttotal time (s)\tstandard deviation");  \
  std::ofstream timingtools_ofile(filename);                                   \
  timingtools_ofile << "# File generated on " << Utilities::get_timestamp()    \
                    << "\n#\n";                                                \
  timingtools_ofile << "# System information:\n";                              \
  for (auto timingtools_compiler_iterator = CompilerInfo::begin();             \
       timingtools_compiler_iterator != CompilerInfo::end();                   \
       ++timingtools_compiler_iterator) {                                      \
    timingtools_ofile << "#   " << timingtools_compiler_iterator.get_key()     \
                      << ": " << timingtools_compiler_iterator.get_value()     \
                      << "\n";                                                 \
  }                                                                            \
  timingtools_ofile << "#\n";                                                  \
  timingtools_ofile << "# Number of samples used: " << timingtools_num_sample  \
                    << "\n#\n";                                                \
  timingtools_ofile << "# File contents:\n";                                   \
  timingtools_ofile                                                            \
      << "# number_of_threads\ttotal_time\tstandard_deviation\n";              \
  timingtools_ofile << "# dimensionless\t(s)\t(s)\n";                          \
  for (unsigned char timingtools_current_num_threads = 0;                      \
       timingtools_current_num_threads < timingtools_num_threads;              \
       ++timingtools_current_num_threads) {                                    \
    timingtools_print(                                                         \
        "%u\t%g\t%g", timingtools_current_num_threads + 1,                     \
        timingtools_scaling_array[timingtools_current_num_threads],            \
        timingtools_scaling_standard_deviation                                 \
            [timingtools_current_num_threads]);                                \
    timingtools_ofile                                                          \
        << timingtools_current_num_threads + 1 << "\t"                         \
        << timingtools_scaling_array[timingtools_current_num_threads] << "\t"  \
        << timingtools_scaling_standard_deviation                              \
               [timingtools_current_num_threads]                               \
        << "\n";                                                               \
  }                                                                            \
  }

#endif // TIMINGTOOLS_HPP
