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
 * @file TaskBasedIonizationSimulation.hpp
 *
 * @brief Ionization radiative transfer simulation using a task-based parallel
 * algorithm.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TASKBASEDIONIZATIONSIMULATION_HPP
#define TASKBASEDIONIZATIONSIMULATION_HPP

#include "RandomGenerator.hpp"
#include "Task.hpp"
#include "ThreadSafeVector.hpp"

#include <vector>

class DensityFunction;
class DensityGridWriter;
class DensitySubGrid;
class DensitySubGridCreator;
class MemorySpace;
class TaskQueue;

/**
 * @brief Ionization radiative transfer simulation using a task-based parallel
 * algorithm.
 */
class TaskBasedIonizationSimulation {
private:
  /*! @brief Grid parts. */
  std::vector< DensitySubGrid * > _subgrids;

  /*! @brief PhotonPacket buffers. */
  MemorySpace *_buffers;

  /*! @brief Queues per thread. */
  std::vector< TaskQueue * > _queues;

  /*! @brief General shared queue. */
  TaskQueue *_shared_queue;

  /*! @brief Task space. */
  ThreadSafeVector< Task > *_tasks;

  /*! @brief Random number generator per thread. */
  std::vector< RandomGenerator > _random_generators;

  /*! @brief Grid creator. */
  DensitySubGridCreator *_grid_creator;

public:
  TaskBasedIonizationSimulation(const int_fast32_t num_thread,
                                const std::string parameterfile_name);
  ~TaskBasedIonizationSimulation();

  void initialize(DensityFunction *density_function = nullptr);
  void run(DensityGridWriter *density_grid_writer = nullptr);
};

#endif // TASKBASEDIONIZATIONSIMULATION_HPP
