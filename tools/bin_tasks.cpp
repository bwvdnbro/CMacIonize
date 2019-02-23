/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file bin_tasks.cpp
 *
 * @brief Command line tool to compress (bin) the information contained in a
 * task stat output file.
 *
 * This file can be compiled using
 * ```
 *  g++ -O3 -std=c++11 -o bin_tasks bin_tasks.cpp
 * ```
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Auxiliary structure used to represent a row in the task file.
 */
struct row {
  /*! @brief Rank that executed the task. */
  int_fast32_t rank;

  /*! @brief Thread that executed the task. */
  int_fast32_t thread;

  /*! @brief Start CPU cycle count for the task. */
  uint_fast64_t start;

  /*! @brief End CPU cycle count for the task. */
  uint_fast64_t stop;

  /*! @brief Type of the task. */
  int_fast32_t type;

  /*! @brief Total accumulated CPU cycle count interval spent executing the
   *  task. */
  uint_fast64_t interval;
};

/**
 * @brief Comparison function used to sort the rows.
 *
 * We do a sort on rank, thread and start CPU cycle count, in that order.
 * The result is a list of rank blocks consisting of thread blocks that are
 * sorted in chronological order.
 *
 * @param a Row a.
 * @param b Row b.
 * @return True if a should be before b in the sorted list.
 */
bool combosort(const struct row &a, const struct row &b) {
  // primary sort key
  if (a.rank < b.rank) {
    return true;
  }
  if (b.rank < a.rank) {
    return false;
  }
  // secondary sort key
  if (a.thread < b.thread) {
    return true;
  }
  if (b.thread < a.thread) {
    return false;
  }
  // tertiary sort key
  return a.start < b.start;
}

/**
 * @brief Command line tool to compress (bin) the information contained in a
 * task stat output file.
 *
 * The tool takes one required argument and two optional arguments. The
 * required argument is the name of a file containing task stat output. The
 * tool will parse this file and generate a new file (with the optional name
 * argument as file name, default is NAME_WITHOUT_EXTENSION.bin.EXTENSION) that
 * contains a compressed version of the same task information. We obtain the
 * compressed version by joining consecutive tasks of the same type that are
 * closer spaced in time than some fraction of the total time interval contained
 * in the stat file. This fraction is the second optional argument (default:
 * 0.1%). To compensate for the loss of information about gaps that are hidden
 * in this way, we create an additional output column per task containing the
 * accumulated active CPU cycle count interval.
 * The general idea is that the compressed file will still result in the same
 * task plot, and that the additional column will preserve all information that
 * is lost by the compression.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // parse command line options
  std::string input_name, output_name;
  double max_fraction = 0.001;

  static struct option long_options[] = {
      {"file", required_argument, 0, 'f'},
      {"output", required_argument, 0, 'o'},
      {"fraction", required_argument, 0, 'l'},
      {0, 0, 0, 0}};
  int option_index = 0;
  int c = getopt_long(argc, argv, "f:l:o:", long_options, &option_index);

  switch (c) {
  case 'f':
    input_name = std::string(optarg);
    break;
  case 'o':
    output_name = std::string(optarg);
    break;
  case 'l':
    max_fraction = atof(optarg);
    break;
  case '?':
    return 1;
  default:
    std::cerr << "This tool requires at least one option!" << std::endl;
    return 1;
  }

  if (input_name == "") {
    std::cerr << "No input file name was given!" << std::endl;
    return 1;
  }
  if (output_name == "") {
    const size_t iext = input_name.rfind('.');
    output_name = input_name.substr(0, iext) + ".bin" + input_name.substr(iext);
  }

  // open the input file and put its contents in a vector
  // while at it, find the minimum start CPU cycle and the maximum stop CPU
  // cycle, they determine the total interval
  std::ifstream file(input_name);
  std::string line;
  std::getline(file, line);
  std::vector< struct row > rows;
  uint_fast64_t minstart = 0xffffffffffffffffull;
  uint_fast64_t maxstop = 0;
  while (std::getline(file, line)) {
    std::istringstream linestream(line);
    struct row this_row;
    linestream >> this_row.rank >> this_row.thread >> this_row.start >>
        this_row.stop >> this_row.type;
    rows.push_back(this_row);
    minstart = std::min(minstart, this_row.start);
    maxstop = std::max(maxstop, this_row.stop);
  }
  file.close();

  uint_fast64_t total_interval = maxstop - minstart;
  total_interval *= max_fraction;

  // sort the vector
  std::sort(rows.begin(), rows.end(), combosort);

  // now do the compression:
  // take the first row and set it as active task
  struct row active_task;
  active_task.rank = rows[0].rank;
  active_task.thread = rows[0].thread;
  active_task.start = rows[0].start;
  active_task.stop = rows[0].stop;
  active_task.type = rows[0].type;
  // add the contribution of this task to the accumulated total CPU cycle count
  active_task.interval = active_task.stop - active_task.stop;
  // open the output file for writing
  std::ofstream ofile(output_name);
  ofile << "# rank\tthread\tstart\tstop\ttype\tinterval\n";
  // now loop over all other tasks
  for (size_t i = 1; i < rows.size(); ++i) {
    // check if the task can be added to the active task
    if (rows[i].rank == active_task.rank &&
        rows[i].thread == active_task.thread &&
        rows[i].type == active_task.type &&
        rows[i].start - active_task.stop < total_interval) {
      // it can: overwrite the stop CPU cycle count
      active_task.stop = rows[i].stop;
      // and add this task's contribution to the total interval
      active_task.interval += rows[i].stop - rows[i].start;
    } else {
      // it cannot: write the active task to the output file
      ofile << active_task.rank << "\t" << active_task.thread << "\t"
            << active_task.start << "\t" << active_task.stop << "\t"
            << active_task.type << "\t" << active_task.interval << "\n";
      // set the new active task to the current task
      active_task.rank = rows[i].rank;
      active_task.thread = rows[i].thread;
      active_task.start = rows[i].start;
      active_task.stop = rows[i].stop;
      active_task.type = rows[i].type;
      active_task.interval = rows[i].stop - rows[i].start;
    }
  }
  // output the final task
  ofile << active_task.rank << "\t" << active_task.thread << "\t"
        << active_task.start << "\t" << active_task.stop << "\t"
        << active_task.type << "\t" << active_task.interval << "\n";

  std::cout << "Done compressing " << input_name << " into " << output_name
            << "." << std::endl;

  return 0;
}
