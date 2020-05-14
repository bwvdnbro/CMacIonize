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
 * @file PhotonTraversalTaskContext.hpp
 *
 * @brief Task context responsible for propagating photon packets.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PHOTONTRAVERSALTASKCONTEXT_HPP
#define PHOTONTRAVERSALTASKCONTEXT_HPP

#include "Abundances.hpp"
#include "ContinuousPhotonSource.hpp"
#include "CrossSections.hpp"
#include "DensitySubGridCreator.hpp"
#include "DiffuseReemissionHandler.hpp"
#include "MemorySpace.hpp"
#include "PhotonPacketStatistics.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"
#include "Task.hpp"
#include "TaskQueue.hpp"

/**
 * @brief Task context responsible for propagating photon packets.
 */
class PhotonTraversalTaskContext {
private:
  /*! @brief Photon buffer array. */
  MemorySpace &_buffers;

  /*! @brief Grid creator. */
  DensitySubGridCreator< DensitySubGrid > &_grid_creator;

  /*! @brief Task space. */
  ThreadSafeVector< Task > &_tasks;

  /*! @brief Number of photon packets that has been terminated. */
  AtomicValue< uint_fast32_t > &_num_photon_done;

  /*! @brief Statistical information about photon packets. */
  PhotonPacketStatistics &_statistics;

  /*! @brief Whether or not to store absorbed photon packets for reemission. */
  const bool _do_reemission;

public:
  /**
   * @brief Constructor.
   *
   * @param buffers Photon buffer array.
   * @param grid_creator Grid creator.
   * @param tasks Task space.
   * @param num_photon_done Number of photon packets that has been terminated.
   * @param statistics Statistical information about photon packets.
   * @param do_reemission Whether or not to store absorbed photon packets for
   * reemission.
   */
  inline PhotonTraversalTaskContext(
      MemorySpace &buffers,
      DensitySubGridCreator< DensitySubGrid > &grid_creator,
      ThreadSafeVector< Task > &tasks,
      AtomicValue< uint_fast32_t > &num_photon_done,
      PhotonPacketStatistics &statistics, const bool do_reemission)
      : _buffers(buffers), _grid_creator(grid_creator), _tasks(tasks),
        _num_photon_done(num_photon_done), _statistics(statistics),
        _do_reemission(do_reemission) {}

  /**
   * @brief Execute a photon traversal task.
   *
   * @param thread_id ID of the thread that executes the task.
   * @param local_buffers Local temporary buffers.
   * @param local_buffer_flags Flags showing which local buffers are in use.
   * @param tasks_to_add Array with indices of newly created tasks.
   * @param queues_to_add Array with target queue indices for the newly created
   * tasks.
   * @param task Task to execute.
   * @return Number of new tasks created by the task.
   */
  inline uint_fast32_t execute(const int_fast8_t thread_id,
                               PhotonBuffer *local_buffers,
                               bool *local_buffer_flags,
                               uint_fast32_t *tasks_to_add,
                               int_fast32_t *queues_to_add, Task &task) {

    uint_fast64_t task_start, task_stop;
    cpucycle_tick(task_start);

    const uint_fast32_t current_buffer_index = task.get_buffer();
    PhotonBuffer &photon_buffer = _buffers[current_buffer_index];
    const uint_fast32_t igrid = photon_buffer.get_subgrid_index();
    DensitySubGrid &this_grid = *_grid_creator.get_subgrid(igrid);

    // prepare output buffers: make sure they are empty and that buffers
    // corresponding to directions outside the simulation box are
    // disabled
    for (int_fast8_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      const uint_fast32_t ngb = this_grid.get_neighbour(i);
      if (ngb != NEIGHBOUR_OUTSIDE) {
        local_buffer_flags[i] = true;
        local_buffers[i].reset();
      } else {
        local_buffer_flags[i] = false;
      }
    }

    // if reemission is disabled, disable output to the internal buffer
    if (!_do_reemission) {
      local_buffer_flags[TRAVELDIRECTION_INSIDE] = false;
    }

    // keep track of the original number of photons
    uint_fast32_t num_photon_done_now = photon_buffer.size();

    // now loop over the input buffer photons and traverse them one by
    // one
    for (uint_fast32_t i = 0; i < photon_buffer.size(); ++i) {

      // active photon
      PhotonPacket &photon = photon_buffer[i];

      // make sure the photon is moving in *a* direction
      cmac_assert_message(photon.get_direction()[0] != 0. ||
                              photon.get_direction()[1] != 0. ||
                              photon.get_direction()[2] != 0.,
                          "size: %" PRIuFAST32, photon_buffer.size());

      // traverse the photon through the active subgrid
      const int_fast32_t result =
          this_grid.interact(photon, photon_buffer.get_direction());

      // check that the photon ended up in a valid output buffer
      cmac_assert_message(result >= 0 && result < TRAVELDIRECTION_NUMBER,
                          "fail");

      // add the photon to an output buffer, if it still exists (if the
      // corresponding output buffer does not exist, this means the
      // photon left the simulation box)
      if (local_buffer_flags[result]) {
        // get the correct output buffer
        PhotonBuffer &output_buffer = local_buffers[result];

        // add the photon
        const uint_fast32_t index = output_buffer.get_next_free_photon();
        output_buffer[index] = photon;
      } else {
        if (result == 0) {
          _statistics.absorb_photon(photon);
        } else {
          _statistics.escape_photon(photon);
        }
      }
    }

    // add none empty buffers to the appropriate queues
    uint_fast8_t largest_index = TRAVELDIRECTION_NUMBER;
    uint_fast32_t largest_size = 0;
    uint_fast32_t num_tasks_to_add = 0;
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {

      // only process enabled, non-empty output buffers
      if (local_buffer_flags[i] && local_buffers[i].size() > 0) {

        // photon packets that are still present in an output buffer
        // are not done yet
        num_photon_done_now -= local_buffers[i].size();

        // move photon packets from the local temporary buffer (that is
        // guaranteed to be large enough) to the actual output buffer
        // for that direction (which might cause on overflow)
        const uint_fast32_t ngb = this_grid.get_neighbour(i);
        uint_fast32_t new_index = this_grid.get_active_buffer(i);

        if (new_index == NEIGHBOUR_OUTSIDE) {
          // buffer was not created yet: create it now
          new_index = _buffers.get_free_buffer();
          PhotonBuffer &buffer = _buffers[new_index];
          buffer.set_subgrid_index(ngb);
          buffer.set_direction(TravelDirections::output_to_input_direction(i));
          this_grid.set_active_buffer(i, new_index);
        }

        uint_fast32_t add_index =
            _buffers.add_photons(new_index, local_buffers[i]);

        // check if the original buffer is full
        if (add_index != new_index) {

          // new_buffers.add_photons already created a new empty
          // buffer, set it as the active buffer for this output
          // direction
          if (_buffers[add_index].size() == 0) {
            _buffers.free_buffer(add_index);
            this_grid.set_active_buffer(i, NEIGHBOUR_OUTSIDE);
          } else {
            this_grid.set_active_buffer(i, add_index);

            cmac_assert_message(_buffers[add_index].get_subgrid_index() == ngb,
                                "Wrong subgrid");
            cmac_assert_message(
                _buffers[add_index].get_direction() ==
                    TravelDirections::output_to_input_direction(i),
                "Wrong direction");
          }

          // YES: create a task for the buffer and add it to the queue
          // the task type depends on the buffer: photon packets in the
          // internal buffer were absorbed and could be reemitted,
          // photon packets in the other buffers left the subgrid and
          // need to be traversed in the neighbouring subgrid
          if (i > 0) {
            DensitySubGrid &subgrid = *_grid_creator.get_subgrid(
                _buffers[new_index].get_subgrid_index());
            const size_t task_index = _tasks.get_free_element();
            Task &new_task = _tasks[task_index];
            new_task.set_subgrid(_buffers[new_index].get_subgrid_index());
            new_task.set_buffer(new_index);
            new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);

            // add dependencies for task:
            //  - subgrid
            new_task.set_dependency(subgrid.get_dependency());

            // add the task to the queue of the corresponding thread
            const int_fast32_t queue_index =
                (*_grid_creator.get_subgrid(ngb)).get_owning_thread();
            queues_to_add[num_tasks_to_add] = queue_index;
            tasks_to_add[num_tasks_to_add] = task_index;
            ++num_tasks_to_add;
          } else {
            const size_t task_index = _tasks.get_free_element();
            Task &new_task = _tasks[task_index];
            new_task.set_subgrid(_buffers[new_index].get_subgrid_index());
            new_task.set_buffer(new_index);
            new_task.set_type(TASKTYPE_PHOTON_REEMIT);
            // a reemit task has no direct dependencies
            // add the task to the general queue
            queues_to_add[num_tasks_to_add] = -1;
            tasks_to_add[num_tasks_to_add] = task_index;
            ++num_tasks_to_add;
          }

        } // if (add_index != new_index)

      } // if (local_buffer_flags[i] &&
      //     local_buffers[i]._actual_size > 0)

      // we have to do this outside the other condition, as buffers to
      // which nothing was added can still be non-empty...
      if (local_buffer_flags[i]) {
        uint_fast32_t new_index = this_grid.get_active_buffer(i);
        if (new_index != NEIGHBOUR_OUTSIDE &&
            _buffers[new_index].size() > largest_size) {
          largest_index = i;
          largest_size = _buffers[new_index].size();
        }
      }

    } // for (int i = TRAVELDIRECTION_NUMBER - 1; i >= 0; --i)

    this_grid.set_largest_buffer(largest_index, largest_size);

    // add photons that were absorbed (if reemission was disabled) or
    // that left the system to the global count
    _num_photon_done.pre_add(num_photon_done_now);

    // delete the original buffer, as we are done with it
    _buffers.free_buffer(current_buffer_index);

    // log the end time of the task
    cpucycle_tick(task_stop);
    this_grid.add_computational_cost(task_stop - task_start);

    return num_tasks_to_add;
  }
};

#endif // PHOTONTRAVERSALTASKCONTEXT_HPP
