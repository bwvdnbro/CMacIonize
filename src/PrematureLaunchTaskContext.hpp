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
 * @file PrematureLaunchTaskContext.hpp
 *
 * @brief Task context responsible for prematurely launching photon buffers.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PREMATURELAUNCHTASKCONTEXT_HPP
#define PREMATURELAUNCHTASKCONTEXT_HPP

#include "DensitySubGridCreator.hpp"
#include "MemorySpace.hpp"
#include "Task.hpp"
#include "TaskQueue.hpp"

/**
 * @brief Task context responsible for prematurely launching photon buffers.
 */
class PrematureLaunchTaskContext {
private:
  /*! @brief Photon buffer array. */
  MemorySpace &_buffers;

  /*! @brief Grid creator. */
  DensitySubGridCreator< DensitySubGrid > &_grid_creator;

  /*! @brief Task space. */
  ThreadSafeVector< Task > &_tasks;

  /*! @brief Queues per thread. */
  std::vector< TaskQueue * > &_queues;

  /*! @brief General shared queue. */
  TaskQueue &_shared_queue;

public:
  /**
   * @brief Constructor.
   *
   * @param buffers Photon buffer array.
   * @param grid_creator Grid creator.
   * @param tasks Task space.
   * @param queues Thread queues.
   * @param shared_queue Shared queue.
   */
  inline PrematureLaunchTaskContext(
      MemorySpace &buffers,
      DensitySubGridCreator< DensitySubGrid > &grid_creator,
      ThreadSafeVector< Task > &tasks, std::vector< TaskQueue * > &queues,
      TaskQueue &shared_queue)
      : _buffers(buffers), _grid_creator(grid_creator), _tasks(tasks),
        _queues(queues), _shared_queue(shared_queue) {}

  /**
   * @brief Execute a premature launch task.
   */
  inline void execute() {

    uint_fast32_t threshold_size = PHOTONBUFFER_SIZE;
    while (threshold_size > 0) {
      threshold_size >>= 1;
      for (auto gridit = _grid_creator.begin();
           gridit != _grid_creator.all_end(); ++gridit) {
        DensitySubGrid &this_subgrid = *gridit;
        if (this_subgrid.get_largest_buffer_size() > threshold_size &&
            this_subgrid.get_dependency()->try_lock()) {

          const uint_fast8_t largest_index =
              this_subgrid.get_largest_buffer_index();
          if (largest_index != TRAVELDIRECTION_NUMBER) {

            const uint_fast32_t non_full_index =
                this_subgrid.get_active_buffer(largest_index);
            this_subgrid.set_active_buffer(largest_index, NEIGHBOUR_OUTSIDE);

            const size_t task_index = _tasks.get_free_element();
            Task &new_task = _tasks[task_index];
            new_task.set_subgrid(_buffers[non_full_index].get_subgrid_index());
            new_task.set_buffer(non_full_index);
            if (largest_index > 0) {
              DensitySubGrid &subgrid = *_grid_creator.get_subgrid(
                  _buffers[non_full_index].get_subgrid_index());
              new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);

              // add dependency
              new_task.set_dependency(subgrid.get_dependency());

              const uint_fast32_t queue_index = subgrid.get_owning_thread();
              _queues[queue_index]->add_task(task_index);
            } else {
              new_task.set_type(TASKTYPE_PHOTON_REEMIT);
              // a reemit task has no dependencies
              _shared_queue.add_task(task_index);
            }

            // set the new largest index
            uint_fast8_t new_largest_index = TRAVELDIRECTION_NUMBER;
            uint_fast32_t new_largest_size = 0;
            for (uint_fast8_t ibuffer = 0; ibuffer < TRAVELDIRECTION_NUMBER;
                 ++ibuffer) {
              if (this_subgrid.get_active_buffer(ibuffer) !=
                      NEIGHBOUR_OUTSIDE &&
                  _buffers[this_subgrid.get_active_buffer(ibuffer)].size() >
                      new_largest_size) {
                new_largest_index = ibuffer;
                new_largest_size =
                    _buffers[this_subgrid.get_active_buffer(ibuffer)].size();
              }
            }
            this_subgrid.set_largest_buffer(new_largest_index,
                                            new_largest_size);

            // unlock the subgrid, we are done with it
            this_subgrid.get_dependency()->unlock();

            // we managed to activate a buffer, we are done
            threshold_size = 0;
            break;
          } else {
            // no semi-full buffers for this subgrid: release the lock
            // again
            this_subgrid.get_dependency()->unlock();
          }
        }
      }
    }
  }
};

#endif // PREMATURELAUNCHTASKCONTEXT_HPP
