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
 * @file SourceContinuousPhotonTaskContext.hpp
 *
 * @brief Task context responsible for generating new photon packets that
 * originate from continuous sources.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SOURCECONTINUOUSPHOTONTASKCONTEXT_HPP
#define SOURCECONTINUOUSPHOTONTASKCONTEXT_HPP

#include "Abundances.hpp"
#include "ContinuousPhotonSource.hpp"
#include "CrossSections.hpp"
#include "DensitySubGridCreator.hpp"
#include "MemorySpace.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"
#include "Task.hpp"
#include "TaskContext.hpp"
#include "TaskQueue.hpp"

/**
 * @brief Task context responsible for generating new photon packets that
 * originate from discrete sources.
 */
class SourceContinuousPhotonTaskContext : public TaskContext {
private:
  /*! @brief Continuous source of UV light. */
  ContinuousPhotonSource &_continuous_photon_source;

  /*! @brief Photon buffer array. */
  MemorySpace &_buffers;

  /*! @brief Per thread random generators. */
  std::vector< RandomGenerator > &_random_generators;

  /*! @brief Weight of an individual continuous photon packet. */
  const double _continuous_photon_weight;

  /*! @brief Spectrum for continuous sources. */
  const PhotonSourceSpectrum &_photon_source_spectrum;

  /*! @brief Abundances. */
  const Abundances _abundances;

  /*! @brief Cross sections for photoionization. */
  const CrossSections &_cross_sections;

  /*! @brief Grid creator. */
  DensitySubGridCreator< DensitySubGrid > &_grid_creator;

  /*! @brief Task space. */
  ThreadSafeVector< Task > &_tasks;

  /*! @brief Buffers used to store sourced photon packets. */
  std::vector< std::vector< PhotonBuffer > > &_continuous_buffers;

  /*! @brief Queues per thread. */
  std::vector< TaskQueue * > &_queues;

  /*! @brief General shared queue. */
  TaskQueue &_shared_queue;

  /*! @brief Counter for the number of continous photon packets that has already
   *  been sourced. */
  AtomicValue< uint_fast32_t > _number_of_continuous_photons;

  /*! @brief Control variable used to make sure that only one thread can flush
   *  the buffers. */
  AtomicValue< uint_fast32_t > _continuous_photons_flushed;

  /*! @brief Locks for the buffers. */
  std::vector< ThreadLock > &_continuous_source_lock;

public:
  /**
   * @brief Constructor.
   *
   * @param continuous_photon_source Continuous photon source.
   * @param buffers Photon buffer array.
   * @param random_generators Per thread random generator.
   * @param continuous_photon_weight Weight of an individual continuous photon
   * packet.
   * @param photon_source_spectrum Spectrum for continuous sources.
   * @param abundances Abundances.
   * @param cross_sections Cross sections for photoionization.
   * @param grid_creator Grid creator.
   * @param tasks Task space.
   * @param continuous_buffers Buffers used to store continuous photon packets.
   * @param queues Queues per thread.
   * @param shared_queue Shared task queue.
   * @param number_of_continuous_photons Number of continuous photon packets
   * to emit.
   * @param continuous_source_lock Locks for the buffers.
   */
  inline SourceContinuousPhotonTaskContext(
      ContinuousPhotonSource &continuous_photon_source, MemorySpace &buffers,
      std::vector< RandomGenerator > &random_generators,
      const double continuous_photon_weight,
      PhotonSourceSpectrum &photon_source_spectrum,
      const Abundances &abundances, CrossSections &cross_sections,
      DensitySubGridCreator< DensitySubGrid > &grid_creator,
      ThreadSafeVector< Task > &tasks,
      std::vector< std::vector< PhotonBuffer > > &continuous_buffers,
      std::vector< TaskQueue * > &queues, TaskQueue &shared_queue,
      const uint_fast32_t number_of_continuous_photons,
      std::vector< ThreadLock > &continuous_source_lock)
      : _continuous_photon_source(continuous_photon_source), _buffers(buffers),
        _random_generators(random_generators),
        _continuous_photon_weight(continuous_photon_weight),
        _photon_source_spectrum(photon_source_spectrum),
        _abundances(abundances), _cross_sections(cross_sections),
        _grid_creator(grid_creator), _tasks(tasks),
        _continuous_buffers(continuous_buffers), _queues(queues),
        _shared_queue(shared_queue),
        _number_of_continuous_photons(number_of_continuous_photons),
        _continuous_photons_flushed(0),
        _continuous_source_lock(continuous_source_lock) {}

  /**
   * @brief Execute a continuous photon source task.
   *
   * @param thread_id ID of the thread that executes the task.
   * @param thread_context Thread specific context for the task.
   * @param tasks_to_add Array with indices of newly created tasks.
   * @param queues_to_add Array with target queue indices for the newly created
   * tasks.
   * @param task Task to execute.
   * @return Number of new tasks created by the task.
   */
  virtual uint_fast32_t execute(const int_fast32_t thread_id,
                                ThreadContext *thread_context,
                                uint_fast32_t *tasks_to_add,
                                int_fast32_t *queues_to_add, Task &task) {

    const uint_fast32_t source_copy = task.get_subgrid();
    const size_t num_photon_this_loop = task.get_buffer();

    // draw random photons and store them in the continuous buffers
    for (uint_fast32_t i = 0; i < num_photon_this_loop; ++i) {

      auto posdir = _continuous_photon_source.get_random_incoming_direction(
          _random_generators[thread_id]);

      const size_t subgrid_index =
          _grid_creator.get_subgrid(posdir.first).get_index();

      PhotonBuffer &active_buffer =
          _continuous_buffers[source_copy][subgrid_index];
      const uint_fast32_t active_index = active_buffer.get_next_free_photon();

      PhotonPacket &photon = active_buffer[active_index];

      photon.set_type(PHOTONTYPE_PRIMARY);
      photon.set_scatter_counter(0);

      // initial position: we currently assume a single source at the
      // origin
      photon.set_position(posdir.first);
      photon.set_direction(posdir.second);

      // we currently assume equal weight for all photons
      photon.set_weight(_continuous_photon_weight);

      // target optical depth (exponential distribution)
      photon.set_target_optical_depth(
          -std::log(_random_generators[thread_id].get_uniform_random_double()));

      const double frequency = _photon_source_spectrum.get_random_frequency(
          _random_generators[thread_id]);
      photon.set_energy(frequency);
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        double sigma = _cross_sections.get_cross_section(ion, frequency);
#ifndef VARIABLE_ABUNDANCES
        if (ion != ION_H_n) {
          sigma *= _abundances.get_abundance(get_element(ion));
        }
#endif
        photon.set_photoionization_cross_section(ion, sigma);
      }

      // did the photon make the buffer overflow?
      if (active_buffer.size() == PHOTONBUFFER_SIZE) {
        // yes: send the buffer off!
        uint_fast32_t buffer_index = _buffers.get_free_buffer();
        PhotonBuffer &input_buffer = _buffers[buffer_index];

        // set general buffer information
        input_buffer.grow(PHOTONBUFFER_SIZE);
        input_buffer.set_subgrid_index(subgrid_index);
        input_buffer.set_direction(TRAVELDIRECTION_INSIDE);

        // copy over the photons
        for (uint_fast32_t iphoton = 0; iphoton < PHOTONBUFFER_SIZE;
             ++iphoton) {
          input_buffer[iphoton] = active_buffer[iphoton];
        }

        // reset the active buffer
        active_buffer.reset();

        // add to the queue of the corresponding thread
        DensitySubGrid &subgrid = *_grid_creator.get_subgrid(subgrid_index);
        const size_t task_index = _tasks.get_free_element();
        Task &new_task = _tasks[task_index];
        new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
        new_task.set_subgrid(subgrid_index);
        new_task.set_buffer(buffer_index);

        // add dependency for task:
        //  - subgrid
        // (the output buffers belong to the subgrid and do not count
        // as a dependency)
        new_task.set_dependency(subgrid.get_dependency());

        _queues[subgrid.get_owning_thread()]->add_task(task_index);
      }
    }

    cmac_assert(_number_of_continuous_photons.value() >= num_photon_this_loop);
    _number_of_continuous_photons.pre_subtract(num_photon_this_loop);

    // are we done with all photons?
    const uint_fast32_t number_of_continuous_blocks = _queues.size();
    if (_number_of_continuous_photons.value() == 0) {
      const uint_fast32_t flush = _continuous_photons_flushed.pre_increment();
      if (flush == 1) {
        // yes: add tasks to flush all buffers
        for (uint_fast32_t other_copy = 0;
             other_copy < number_of_continuous_blocks; ++other_copy) {
          const size_t task_index = _tasks.get_free_element();
          Task &new_task = _tasks[task_index];
          new_task.set_type(TASKTYPE_FLUSH_CONTINUOUS_PHOTON_BUFFERS);
          new_task.set_subgrid(other_copy);
          new_task.set_buffer(0);

          // add dependency for task:
          //  - subgrid
          // (the output buffers belong to the subgrid and do not
          // count as a dependency)
          new_task.set_dependency(&_continuous_source_lock[other_copy]);

          _shared_queue.add_task(task_index);
        }
      }
    }

    return 0;
  }
};

#endif // SOURCECONTINUOUSPHOTONTASKCONTEXT_HPP
