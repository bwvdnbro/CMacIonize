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
 * @file TaskBasedIonizationSimulation.cpp
 *
 * @brief TaskBasedIonizationSimulation implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "TaskBasedIonizationSimulation.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DensitySubGrid.hpp"
#include "DensitySubGridCreator.hpp"
#include "MemorySpace.hpp"
#include "ParameterFile.hpp"
#include "SimulationBox.hpp"
#include "TaskQueue.hpp"

#include <fstream>
#include <omp.h>

/**
 * @brief Write a file with the start and end times of all tasks.
 *
 * @param iloop Iteration number (added to file name).
 * @param tasks Tasks to print.
 * @param iteration_start Start CPU cycle count of the iteration on this
 * process.
 * @param iteration_end End CPU cycle count of the iteration on this process.
 */
inline void output_tasks(const uint_fast32_t iloop,
                         ThreadSafeVector< Task > &tasks,
                         const uint_fast64_t iteration_start,
                         const uint_fast64_t iteration_end) {

  {
    // compose the file name
    std::stringstream filename;
    filename << "tasks_";
    filename.fill('0');
    filename.width(2);
    filename << iloop;
    filename << ".txt";

    // now open the file
    std::ofstream ofile(filename.str(), std::ofstream::trunc);

    ofile << "# rank\tthread\tstart\tstop\ttype\n";

    // write the start and end CPU cycle count
    // this is a dummy task executed by thread 0 (so that the min or max
    // thread count is not affected), but with non-existing type -1
    ofile << "0\t0\t" << iteration_start << "\t" << iteration_end << "\t-1\n";

    // write the task info
    const size_t tsize = tasks.size();
    for (size_t i = 0; i < tsize; ++i) {
      const Task &task = tasks[i];
      cmac_assert_message(task.done(), "Task was never executed!");
      int_fast8_t type;
      int_fast32_t thread_id;
      uint_fast64_t start, end;
      task.get_timing_information(type, thread_id, start, end);
      ofile << "0\t" << thread_id << "\t" << start << "\t" << end << "\t"
            << static_cast< int_fast32_t >(type) << "\n";
    }
  }
}

/**
 * @brief Constructor.
 *
 * This method will read the following parameters from the parameter file:
 *  - number of buffers: Number of photon packets buffers to allocate in memory
 *    (default: 50000)
 *  - queue size per thread: Size of the queue for a single thread (default:
 *    10000)
 *  - shared queue size: Size of the shared queue (default: 100000)
 *  - number of tasks: Number of tasks to allocate in memory (default: 500000)
 *  - random seed: Seed used to initialize the random number generator (default:
 *    42)
 *
 * @param num_thread Number of shared memory parallel threads to use.
 * @param parameterfile_name Name of the parameter file to use.
 */
TaskBasedIonizationSimulation::TaskBasedIonizationSimulation(
    const int_fast32_t num_thread, const std::string parameterfile_name)
    : _parameter_file(parameterfile_name), _simulation_box(_parameter_file) {

  omp_set_num_threads(num_thread);

  const size_t number_of_buffers = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:number of buffers", 50000);
  _buffers = new MemorySpace(number_of_buffers);
  const size_t queue_size_per_thread = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:queue size per thread", 10000);
  _queues.resize(num_thread);
  for (int_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    _queues[ithread] = new TaskQueue(queue_size_per_thread);
  }
  const size_t shared_queue_size = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:shared queue size", 100000);
  _shared_queue = new TaskQueue(shared_queue_size);
  const size_t number_of_tasks = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:number of tasks", 500000);
  _tasks = new ThreadSafeVector< Task >(number_of_tasks);
  _random_generators.resize(num_thread);
  const int_fast32_t random_seed = _parameter_file.get_value< int_fast32_t >(
      "TaskBasedIonizationSimulation:random seed", 42);
  for (uint_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    _random_generators[ithread].set_seed(random_seed + ithread);
  }
  _grid_creator =
      new DensitySubGridCreator(_simulation_box.get_box(), _parameter_file);
  _subgrids.resize(_grid_creator->number_of_subgrids(), nullptr);

  _density_function =
      DensityFunctionFactory::generate(_parameter_file, nullptr);
  _density_grid_writer =
      DensityGridWriterFactory::generate(".", _parameter_file, false, nullptr);
}

/**
 * @brief Destructor.
 */
TaskBasedIonizationSimulation::~TaskBasedIonizationSimulation() {
  for (uint_fast32_t igrid = 0; igrid < _subgrids.size(); ++igrid) {
    delete _subgrids[igrid];
  }
  delete _buffers;
  for (uint_fast8_t ithread = 0; ithread < _queues.size(); ++ithread) {
    delete _queues[ithread];
  }
  delete _shared_queue;
  delete _tasks;
  delete _grid_creator;
  delete _density_function;
  delete _density_grid_writer;
}

/**
 * @brief Initialize a photoionization simulation.
 *
 * @param density_function DensityFunction to use to initialize the density
 * field.
 */
void TaskBasedIonizationSimulation::initialize(
    DensityFunction *density_function) {

  if (density_function == nullptr) {
    density_function = _density_function;
  }

#pragma omp parallel for default(shared)
  for (uint_fast32_t i = 0; i < _subgrids.size(); ++i) {
    _subgrids[i] = _grid_creator->create_subgrid(i);
    _subgrids[i]->set_owning_thread(omp_get_thread_num());
    for (auto it = _subgrids[i]->begin(); it != _subgrids[i]->end(); ++it) {
      DensityValues values = (*density_function)(it);
      it.set_number_density(values.get_number_density());
      it.set_neutral_fraction(values.get_ionic_fraction(ION_H_n));
    }
  }
}

/**
 * @brief Run a photoionization simulation.
 *
 * @param density_grid_writer DensityGridWriter to use to output the results.
 */
void TaskBasedIonizationSimulation::run(
    DensityGridWriter *density_grid_writer) {

  if (density_grid_writer == nullptr) {
    density_grid_writer = _density_grid_writer;
  }

  for (uint_fast8_t iloop = 0; iloop < 10; ++iloop) {

    uint_fast64_t iteration_start, iteration_end;
    cpucycle_tick(iteration_start);

    cmac_warning("Loop: %" PRIuFAST8, iloop);

    _tasks->get_free_elements(5000);
#pragma omp parallel for default(shared)
    for (uint_fast32_t i = 0; i < 5000; ++i) {
      (*_tasks)[i].set_type(TASKTYPE_SOURCE_PHOTON);
      (*_tasks)[i].set_buffer(PHOTONBUFFER_SIZE);
    }
    _shared_queue->add_tasks(0, 5000);

    bool global_run_flag = true;
    AtomicValue< uint_fast32_t > num_empty(TRAVELDIRECTION_NUMBER *
                                           _subgrids.size());
    AtomicValue< uint_fast32_t > num_active_buffers(0);
    AtomicValue< uint_fast32_t > num_photon_done(0);
    const uint_fast32_t num_empty_target =
        TRAVELDIRECTION_NUMBER * _subgrids.size();
#pragma omp parallel default(shared)
    {
      // thread initialisation
      const int_fast8_t thread_id = omp_get_thread_num();
      PhotonBuffer local_buffers[TRAVELDIRECTION_NUMBER];
      bool local_buffer_flags[TRAVELDIRECTION_NUMBER];
      for (int_fast8_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
        local_buffers[i].set_direction(
            TravelDirections::output_to_input_direction(i));
        local_buffers[i].reset();
        local_buffer_flags[i] = true;
      }

      // actual run flag
      while (global_run_flag) {
        uint_fast32_t current_index = _queues[thread_id]->get_task(*_tasks);
        if (current_index == NO_TASK) {
          current_index = _shared_queue->get_task(*_tasks);
        }

        if (current_index == NO_TASK) {
          uint_fast32_t threshold_size = PHOTONBUFFER_SIZE;
          while (threshold_size > 0) {
            threshold_size >>= 1;
            for (uint_fast32_t igrid = 0; igrid < _subgrids.size(); ++igrid) {
              if (_subgrids[igrid]->get_largest_buffer_size() >
                      threshold_size &&
                  _subgrids[igrid]->get_dependency()->try_lock()) {

                const uint_fast8_t largest_index =
                    _subgrids[igrid]->get_largest_buffer_index();
                if (largest_index != TRAVELDIRECTION_NUMBER) {

                  const uint_fast32_t non_full_index =
                      _subgrids[igrid]->get_active_buffer(largest_index);
                  const uint_fast32_t new_index = _buffers->get_free_buffer();
                  (*_buffers)[new_index].set_subgrid_index(
                      (*_buffers)[non_full_index].get_subgrid_index());
                  (*_buffers)[new_index].set_direction(
                      (*_buffers)[non_full_index].get_direction());
                  _subgrids[igrid]->set_active_buffer(largest_index, new_index);
                  // we are creating a new active photon buffer
                  num_active_buffers.pre_increment();
                  // we created a new empty buffer
                  num_empty.pre_increment();

                  const size_t task_index = _tasks->get_free_element();
                  Task &new_task = (*_tasks)[task_index];
                  new_task.set_subgrid(
                      (*_buffers)[non_full_index].get_subgrid_index());
                  new_task.set_buffer(non_full_index);
                  if (largest_index > 0) {
                    DensitySubGrid &subgrid =
                        *_subgrids[(*_buffers)[non_full_index]
                                       .get_subgrid_index()];
                    new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);

                    // add dependency
                    new_task.set_dependency(subgrid.get_dependency());

                    const uint_fast32_t queue_index =
                        _subgrids[(*_buffers)[non_full_index]
                                      .get_subgrid_index()]
                            ->get_owning_thread();
                    _queues[queue_index]->add_task(task_index);
                  } else {
                    new_task.set_type(TASKTYPE_PHOTON_REEMIT);
                    // a reemit task has no dependencies
                    _shared_queue->add_task(task_index);
                  }

                  // set the new largest index
                  uint_fast8_t new_largest_index = TRAVELDIRECTION_NUMBER;
                  uint_fast32_t new_largest_size = 0;
                  for (uint_fast8_t ibuffer = 0;
                       ibuffer < TRAVELDIRECTION_NUMBER; ++ibuffer) {
                    if (_subgrids[igrid]->get_active_buffer(ibuffer) !=
                            NEIGHBOUR_OUTSIDE &&
                        (*_buffers)[_subgrids[igrid]->get_active_buffer(
                                        ibuffer)]
                                .size() > new_largest_size) {
                      new_largest_index = ibuffer;
                      new_largest_size =
                          (*_buffers)[_subgrids[igrid]->get_active_buffer(
                                          ibuffer)]
                              .size();
                    }
                  }
                  _subgrids[igrid]->set_largest_buffer(new_largest_index,
                                                       new_largest_size);

                  // unlock the subgrid, we are done with it
                  _subgrids[igrid]->get_dependency()->unlock();

                  // we managed to activate a buffer, we are done with this
                  // function
                  break;
                } else {
                  // no semi-full buffers for this subgrid: release the lock
                  // again
                  _subgrids[igrid]->get_dependency()->unlock();
                }
              }
            }
          }
          current_index = _queues[thread_id]->get_task(*_tasks);
          if (current_index == NO_TASK) {
            current_index = _shared_queue->get_task(*_tasks);
          }
        }

        while (current_index != NO_TASK) {

          // execute task
          uint_fast32_t num_tasks_to_add = 0;
          uint_fast32_t tasks_to_add[TRAVELDIRECTION_NUMBER];
          int_fast32_t queues_to_add[TRAVELDIRECTION_NUMBER];

          Task &task = (*_tasks)[current_index];

          if (task.get_type() == TASKTYPE_SOURCE_PHOTON) {

            task.start(thread_id);
            num_active_buffers.pre_increment();
            uint_fast32_t num_photon_this_loop = task.get_buffer();

            // get a free photon buffer in the central queue
            uint_fast32_t buffer_index = (*_buffers).get_free_buffer();
            PhotonBuffer &input_buffer = (*_buffers)[buffer_index];
            // assign the buffer to a random thread that has a copy of the
            // subgrid that contains the source position. This should ensure
            // a balanced load for these threads.
            uint_fast32_t this_central_index = 4 * 8 * 8 + 4 * 8 + 4;

            // set general buffer information
            input_buffer.grow(num_photon_this_loop);
            input_buffer.set_subgrid_index(this_central_index);
            input_buffer.set_direction(TRAVELDIRECTION_INSIDE);

            // draw random photons and store them in the buffer
            for (uint_fast32_t i = 0; i < num_photon_this_loop; ++i) {

              PhotonPacket &photon = input_buffer[i];

              // initial position: we currently assume a single source at the
              // origin
              photon.set_position(0., 0., 0.);

              // draw two pseudo random numbers
              const double cost = 2. * _random_generators[thread_id]
                                           .get_uniform_random_double() -
                                  1.;
              const double phi =
                  2. * M_PI *
                  _random_generators[thread_id].get_uniform_random_double();

              // now use them to get all directional angles
              const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
              const double cosp = std::cos(phi);
              const double sinp = std::sin(phi);

              // set the direction...
              double direction[3];
              direction[0] = sint * cosp;
              direction[1] = sint * sinp;
              direction[2] = cost;

              photon.set_direction(direction[0], direction[1], direction[2]);

              // we currently assume equal weight for all photons
              photon.set_weight(1.);

              // target optical depth (exponential distribution)
              photon.set_target_optical_depth(-std::log(
                  _random_generators[thread_id].get_uniform_random_double()));

              // this is the fixed cross section we use for the moment
              photon.set_photoionization_cross_section(6.3e-22);
            }

            // add to the queue of the corresponding thread
            DensitySubGrid &subgrid = *_subgrids[this_central_index];
            const size_t task_index = _tasks->get_free_element();
            Task &new_task = (*_tasks)[task_index];
            new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
            new_task.set_subgrid(this_central_index);
            new_task.set_buffer(buffer_index);

            // add dependency for task:
            //  - subgrid
            // (the output buffers belong to the subgrid and do not count as a
            // dependency)
            new_task.set_dependency(subgrid.get_dependency());

            queues_to_add[num_tasks_to_add] = thread_id;
            tasks_to_add[num_tasks_to_add] = task_index;
            ++num_tasks_to_add;

            // log the end time of the task
            task.stop();

          } else if (task.get_type() == TASKTYPE_PHOTON_TRAVERSAL) {

            task.start(thread_id);

            const uint_fast32_t current_buffer_index = task.get_buffer();
            PhotonBuffer &photon_buffer = (*_buffers)[current_buffer_index];
            const uint_fast32_t igrid = photon_buffer.get_subgrid_index();
            DensitySubGrid &this_grid = *_subgrids[igrid];

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
            if (true) {
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
              cmac_assert_message(
                  result >= 0 && result < TRAVELDIRECTION_NUMBER, "fail");

              // add the photon to an output buffer, if it still exists (if the
              // corresponding output buffer does not exist, this means the
              // photon left the simulation box)
              if (local_buffer_flags[result]) {
                // get the correct output buffer
                PhotonBuffer &output_buffer = local_buffers[result];

                // add the photon
                const uint_fast32_t index =
                    output_buffer.get_next_free_photon();
                output_buffer[index] = photon;
              }
            }

            // add none empty buffers to the appropriate queues
            uint_fast8_t largest_index = TRAVELDIRECTION_NUMBER;
            uint_fast32_t largest_size = 0;
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
                  new_index = _buffers->get_free_buffer();
                  PhotonBuffer &buffer = (*_buffers)[new_index];
                  buffer.set_subgrid_index(ngb);
                  buffer.set_direction(
                      TravelDirections::output_to_input_direction(i));
                  this_grid.set_active_buffer(i, new_index);
                }

                if ((*_buffers)[new_index].size() == 0) {
                  // we are adding photons to an empty buffer
                  num_empty.pre_decrement();
                }
                uint_fast32_t add_index =
                    _buffers->add_photons(new_index, local_buffers[i]);

                // check if the original buffer is full
                if (add_index != new_index) {

                  // a new active buffer was created
                  num_active_buffers.pre_increment();

                  // new_buffers.add_photons already created a new empty
                  // buffer, set it as the active buffer for this output
                  // direction
                  this_grid.set_active_buffer(i, add_index);
                  if ((*_buffers)[add_index].size() == 0) {
                    // we have created a new empty buffer
                    num_empty.pre_increment();
                  }

                  // YES: create a task for the buffer and add it to the queue
                  // the task type depends on the buffer: photon packets in the
                  // internal buffer were absorbed and could be reemitted,
                  // photon packets in the other buffers left the subgrid and
                  // need to be traversed in the neighbouring subgrid
                  if (i > 0) {
                    DensitySubGrid &subgrid =
                        *_subgrids[(*_buffers)[new_index].get_subgrid_index()];
                    const size_t task_index = _tasks->get_free_element();
                    Task &new_task = (*_tasks)[task_index];
                    new_task.set_subgrid(
                        (*_buffers)[new_index].get_subgrid_index());
                    new_task.set_buffer(new_index);
                    new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);

                    // add dependencies for task:
                    //  - subgrid
                    new_task.set_dependency(subgrid.get_dependency());

                    // add the task to the queue of the corresponding thread
                    const int_fast32_t queue_index =
                        _subgrids[ngb]->get_owning_thread();
                    queues_to_add[num_tasks_to_add] = queue_index;
                    tasks_to_add[num_tasks_to_add] = task_index;
                    ++num_tasks_to_add;
                  } else {
                    const size_t task_index = _tasks->get_free_element();
                    Task &new_task = (*_tasks)[task_index];
                    new_task.set_subgrid(
                        (*_buffers)[new_index].get_subgrid_index());
                    new_task.set_buffer(new_index);
                    new_task.set_type(TASKTYPE_PHOTON_REEMIT);
                    // a reemit task has no direct dependencies
                    // add the task to the general queue
                    queues_to_add[num_tasks_to_add] = -1;
                    tasks_to_add[num_tasks_to_add] = task_index;
                    ++num_tasks_to_add;
                  }

                  cmac_assert_message(
                      (*_buffers)[add_index].get_subgrid_index() == ngb,
                      "Wrong subgrid");
                  cmac_assert_message(
                      (*_buffers)[add_index].get_direction() ==
                          TravelDirections::output_to_input_direction(i),
                      "Wrong direction");

                } // if (add_index != new_index)

              } // if (local_buffer_flags[i] &&
              //     local_buffers[i]._actual_size > 0)

              // we have to do this outside the other condition, as buffers to
              // which nothing was added can still be non-empty...
              if (local_buffer_flags[i]) {
                uint_fast32_t new_index = this_grid.get_active_buffer(i);
                if (new_index != NEIGHBOUR_OUTSIDE &&
                    (*_buffers)[new_index].size() > largest_size) {
                  largest_index = i;
                  largest_size = (*_buffers)[new_index].size();
                }
              }

            } // for (int i = TRAVELDIRECTION_NUMBER - 1; i >= 0; --i)

            this_grid.set_largest_buffer(largest_index, largest_size);

            // add photons that were absorbed (if reemission was disabled) or
            // that left the system to the global count
            num_photon_done.pre_add(num_photon_done_now);

            // delete the original buffer, as we are done with it
            _buffers->free_buffer(current_buffer_index);

            cmac_assert_message(num_active_buffers.value() > 0,
                                "Number of active buffers < 0!");
            num_active_buffers.pre_decrement();
            // log the end time of the task
            task.stop();
          }
          task.unlock_dependency();

          for (uint_fast32_t itask = 0; itask < num_tasks_to_add; ++itask) {
            if (queues_to_add[itask] < 0) {
              // general queue
              _shared_queue->add_task(tasks_to_add[itask]);
            } else {
              _queues[queues_to_add[itask]]->add_task(tasks_to_add[itask]);
            }
          }

          current_index = _queues[thread_id]->get_task(*_tasks);
          if (current_index == NO_TASK) {
            current_index = _shared_queue->get_task(*_tasks);
          }
        }

        if (num_empty.value() == num_empty_target &&
            num_active_buffers.value() == 0 && num_photon_done.value() == 1e6) {
          global_run_flag = false;
        }
      } // while(global_run_flag)
    }   // parallel region

#pragma omp parallel for default(shared)
    for (uint_fast32_t igrid = 0; igrid < _subgrids.size(); ++igrid) {
      _subgrids[igrid]->compute_neutral_fraction(4.26e49, 1e6);
    }

    cpucycle_tick(iteration_end);
    output_tasks(iloop, *_tasks, iteration_start, iteration_end);

    _tasks->clear();
  } // photoionization loop

  _density_grid_writer->write(_subgrids, _grid_creator->number_of_subgrids(),
                              _grid_creator->number_of_cells(),
                              _simulation_box.get_box(), 10, _parameter_file);
}
