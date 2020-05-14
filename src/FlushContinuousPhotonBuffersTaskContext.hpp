/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file FlushContinuousPhotonBuffersTaskContext.hpp
 *
 * @brief Task context responsible for flushing the continuous source photon
 * buffers.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FLUSHCONTINUOUSPHOTONBUFFERSTASKCONTEXT_HPP
#define FLUSHCONTINUOUSPHOTONBUFFERSTASKCONTEXT_HPP

#include "Abundances.hpp"
#include "ContinuousPhotonSource.hpp"
#include "CrossSections.hpp"
#include "DensitySubGridCreator.hpp"
#include "MemorySpace.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"
#include "Task.hpp"
#include "TaskQueue.hpp"

/**
 * @brief Task context responsible for flushing the continuous source photon
 * buffers.
 */
class FlushContinuousPhotonBuffersTaskContext {
private:
  /*! @brief Photon buffer array. */
  MemorySpace &_buffers;

  /*! @brief Grid creator. */
  DensitySubGridCreator< DensitySubGrid > &_grid_creator;

  /*! @brief Task space. */
  ThreadSafeVector< Task > &_tasks;

  /*! @brief Buffers used to store sourced photon packets. */
  std::vector< std::vector< PhotonBuffer > > &_continuous_buffers;

  /*! @brief Queues per thread. */
  std::vector< TaskQueue * > &_queues;

public:
  /**
   * @brief Constructor.
   *
   * @param buffers Photon buffer array.
   * @param grid_creator Grid creator.
   * @param tasks Task space.
   * @param continuous_buffers Buffers used to store continuous photon packets.
   * @param queues Queues per thread.
   */
  inline FlushContinuousPhotonBuffersTaskContext(
      MemorySpace &buffers,
      DensitySubGridCreator< DensitySubGrid > &grid_creator,
      ThreadSafeVector< Task > &tasks,
      std::vector< std::vector< PhotonBuffer > > &continuous_buffers,
      std::vector< TaskQueue * > &queues)
      : _buffers(buffers), _grid_creator(grid_creator), _tasks(tasks),
        _continuous_buffers(continuous_buffers), _queues(queues) {}

  /**
   * @brief Execute a continuous photon source task.
   *
   * @param thread_id ID of the thread that executes the task.
   * @param tasks_to_add Array with indices of newly created tasks.
   * @param queues_to_add Array with target queue indices for the newly created
   * tasks.
   * @param task Task to execute.
   * @return Number of new tasks created by the task.
   */
  inline uint_fast32_t execute(const int_fast8_t thread_id,
                               uint_fast32_t *tasks_to_add,
                               int_fast32_t *queues_to_add, Task &task) {

    const uint_fast32_t source_copy = task.get_subgrid();
    for (uint_fast32_t i = 0; i < _continuous_buffers[source_copy].size();
         ++i) {
      if (_continuous_buffers[source_copy][i].size() > 0) {
        // yes: send the buffer off!
        uint_fast32_t buffer_index = _buffers.get_free_buffer();
        PhotonBuffer &input_buffer = _buffers[buffer_index];

        // set general buffer information
        input_buffer.grow(_continuous_buffers[source_copy][i].size());
        input_buffer.set_subgrid_index(i);
        input_buffer.set_direction(TRAVELDIRECTION_INSIDE);

        // copy over the photons
        for (uint_fast32_t iphoton = 0;
             iphoton < _continuous_buffers[source_copy][i].size(); ++iphoton) {
          input_buffer[iphoton] = _continuous_buffers[source_copy][i][iphoton];
        }

        // reset the active buffer
        _continuous_buffers[source_copy][i].reset();

        // add to the queue of the corresponding thread
        DensitySubGrid &subgrid = *_grid_creator.get_subgrid(i);
        const size_t task_index = _tasks.get_free_element();
        Task &new_task = _tasks[task_index];
        new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
        new_task.set_subgrid(i);
        new_task.set_buffer(buffer_index);

        // add dependency for task:
        //  - subgrid
        // (the output buffers belong to the subgrid and do not count
        // as a dependency)
        new_task.set_dependency(subgrid.get_dependency());

        _queues[subgrid.get_owning_thread()]->add_task(task_index);
      }
    }

    return 0;
  }
};

#endif // FLUSHCONTINUOUSPHOTONBUFFERSTASKCONTEXT_HPP
