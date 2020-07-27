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
 * @file PhotonReemitTaskContext.hpp
 *
 * @brief Task context responsible for reemitting photon packets.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PHOTONREEMITTASKCONTEXT_HPP
#define PHOTONREEMITTASKCONTEXT_HPP

#include "Abundances.hpp"
#include "ContinuousPhotonSource.hpp"
#include "CrossSections.hpp"
#include "DensitySubGridCreator.hpp"
#include "DiffuseReemissionHandler.hpp"
#include "MemorySpace.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"
#include "Task.hpp"
#include "TaskContext.hpp"
#include "TaskQueue.hpp"

/**
 * @brief Task context responsible for reemitting photon packets.
 */
template < typename _subgrid_type_ >
class PhotonReemitTaskContext : public TaskContext {
private:
  /*! @brief Photon buffer array. */
  MemorySpace &_buffers;

  /*! @brief Per thread random generators. */
  std::vector< RandomGenerator > &_random_generators;

  /*! @brief Reemission handler. */
  DiffuseReemissionHandler &_reemission_handler;

  /*! @brief Abundances. */
  const Abundances _abundances;

  /*! @brief Cross sections for photoionization. */
  const CrossSections &_cross_sections;

  /*! @brief Grid creator. */
  DensitySubGridCreator< _subgrid_type_ > &_grid_creator;

  /*! @brief Task space. */
  ThreadSafeVector< Task > &_tasks;

  /*! @brief Number of photon packets that has been terminated. */
  AtomicValue< uint_fast32_t > &_num_photon_done;

public:
  /**
   * @brief Constructor.
   *
   * @param buffers Photon buffer array.
   * @param random_generators Per thread random number generators.
   * @param reemission_handler Reemission handler.
   * @param abundances Abundances.
   * @param cross_sections Cross sections for photoionization.
   * @param grid_creator Grid creator.
   * @param tasks Task space.
   * @param num_photon_done Number of photon packets that has been terminated.
   */
  inline PhotonReemitTaskContext(
      MemorySpace &buffers, std::vector< RandomGenerator > &random_generators,
      DiffuseReemissionHandler &reemission_handler,
      const Abundances &abundances, const CrossSections &cross_sections,
      DensitySubGridCreator< _subgrid_type_ > &grid_creator,
      ThreadSafeVector< Task > &tasks,
      AtomicValue< uint_fast32_t > &num_photon_done)
      : _buffers(buffers), _random_generators(random_generators),
        _reemission_handler(reemission_handler), _abundances(abundances),
        _cross_sections(cross_sections), _grid_creator(grid_creator),
        _tasks(tasks), _num_photon_done(num_photon_done) {}

  /**
   * @brief Execute a photon reemission task.
   *
   * @param thread_id ID of the thread that executes the task.
   * @param thread_context Thread specific context for the task.
   * @param tasks_to_add Array with indices of newly created tasks.
   * @param queues_to_add Array with target queue indices for the newly created
   * tasks.
   * @param task Task to execute.
   * @return Number of new tasks created by the task.
   */
  inline uint_fast32_t execute(const int_fast32_t thread_id,
                               ThreadContext *thread_context,
                               uint_fast32_t *tasks_to_add,
                               int_fast32_t *queues_to_add, Task &task) {

    const size_t current_buffer_index = task.get_buffer();
    PhotonBuffer &buffer = _buffers[current_buffer_index];

    uint_fast32_t num_photon_done_now = buffer.size();
    DensitySubGrid &subgrid = *_grid_creator.get_subgrid(task.get_subgrid());

    // reemission
    uint_fast32_t index = 0;
    for (uint_fast32_t iphoton = 0; iphoton < buffer.size(); ++iphoton) {
      PhotonPacket &old_photon = buffer[iphoton];
      const IonizationVariables &ionization_variables =
          subgrid.get_cell(old_photon.get_position())
              .get_ionization_variables();
#ifdef HAS_HELIUM
#ifdef VARIABLE_ABUNDANCES
      const double AHe =
          ionization_variables.get_abundances().get_abundance(ELEMENT_He);
#else
      // the helium abundance is already part of the cross section
      const double AHe = 1.;
#endif
#else
      const double AHe = 0.;
#endif
      PhotonType new_type;
      const double new_frequency =
          _reemission_handler.reemit(old_photon, AHe, ionization_variables,
                                     _random_generators[thread_id], new_type);
      if (new_frequency > 0.) {
        PhotonPacket &new_photon = buffer[index];
        new_photon.set_type(new_type);
        new_photon.set_scatter_counter(old_photon.get_scatter_counter() + 1);
        new_photon.set_position(old_photon.get_position());
        new_photon.set_weight(old_photon.get_weight());

        new_photon.set_energy(new_frequency);
        for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
          double sigma = _cross_sections.get_cross_section(ion, new_frequency);
#ifndef VARIABLE_ABUNDANCES
          if (ion != ION_H_n) {
            sigma *= _abundances.get_abundance(get_element(ion));
          }
#endif
          new_photon.set_photoionization_cross_section(ion, sigma);
        }

        // draw two pseudo random numbers
        const double cost =
            2. * _random_generators[thread_id].get_uniform_random_double() - 1.;
        const double phi =
            2. * M_PI *
            _random_generators[thread_id].get_uniform_random_double();

        // now use them to get all directional angles
        const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
        const double cosp = std::cos(phi);
        const double sinp = std::sin(phi);

        // set the direction...
        const CoordinateVector<> direction(sint * cosp, sint * sinp, cost);

        new_photon.set_direction(direction);

        // target optical depth (exponential distribution)
        new_photon.set_target_optical_depth(-std::log(
            _random_generators[thread_id].get_uniform_random_double()));

        ++index;
      }
    }
    // update the size of the buffer to account for photons that were
    // not reemitted
    buffer.grow(index);

    num_photon_done_now -= buffer.size();
    _num_photon_done.pre_add(num_photon_done_now);

    uint_fast32_t num_tasks_to_add = 0;
    if (index > 0) {
      // there are still photon packets left: generate a traversal
      // task
      const size_t task_index = _tasks.get_free_element();
      Task &new_task = _tasks[task_index];
      new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
      new_task.set_subgrid(task.get_subgrid());
      new_task.set_buffer(current_buffer_index);
      new_task.set_dependency(subgrid.get_dependency());

      queues_to_add[num_tasks_to_add] = subgrid.get_owning_thread();
      tasks_to_add[num_tasks_to_add] = task_index;
      ++num_tasks_to_add;
    } else {
      // delete the original buffer, as we are done with it
      _buffers.free_buffer(current_buffer_index);
    }

    return num_tasks_to_add;
  }
};

#endif // PHOTONREEMITTASKCONTEXT_HPP
