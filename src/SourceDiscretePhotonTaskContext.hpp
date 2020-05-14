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
 * @file SourceDiscretePhotonTaskContext.hpp
 *
 * @brief Task context responsible for generating new photon packets that
 * originate from discrete sources.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SOURCEDISCRETEPHOTONTASKCONTEXT_HPP
#define SOURCEDISCRETEPHOTONTASKCONTEXT_HPP

#include "CrossSections.hpp"
#include "DistributedPhotonSource.hpp"
#include "MemorySpace.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "Task.hpp"
#include "TaskContext.hpp"

/**
 * @brief Task context responsible for generating new photon packets that
 * originate from discrete sources.
 */
class SourceDiscretePhotonTaskContext : public TaskContext {
private:
  /*! @brief Discrete photon source. */
  DistributedPhotonSource< DensitySubGrid > &_photon_source;

  /*! @brief Photon buffer array. */
  MemorySpace &_buffers;

  /*! @brief Per thread random generators. */
  std::vector< RandomGenerator > &_random_generators;

  /*! @brief Weight of an individual discrete photon packet. */
  const double _discrete_photon_weight;

  /*! @brief Spectrum for discrete sources. */
  const PhotonSourceSpectrum &_photon_source_spectrum;

  /*! @brief Abundances. */
  const Abundances _abundances;

  /*! @brief Cross sections for photoionization. */
  const CrossSections &_cross_sections;

  /*! @brief Grid creator. */
  DensitySubGridCreator< DensitySubGrid > &_grid_creator;

  /*! @brief Task space. */
  ThreadSafeVector< Task > &_tasks;

public:
  /**
   * @brief Constructor.
   *
   * @param photon_source Discrete photon source.
   * @param buffers Photon buffer array.
   * @param random_generators Per thread random generator.
   * @param discrete_photon_weight Weight of an individual discrete photon
   * packet.
   * @param photon_source_spectrum Spectrum for discrete sources.
   * @param abundances Abundances.
   * @param cross_sections Cross sections for photoionization.
   * @param grid_creator Grid creator.
   * @param tasks Task space.
   */
  inline SourceDiscretePhotonTaskContext(
      DistributedPhotonSource< DensitySubGrid > &photon_source,
      MemorySpace &buffers, std::vector< RandomGenerator > &random_generators,
      const double discrete_photon_weight,
      PhotonSourceSpectrum &photon_source_spectrum,
      const Abundances &abundances, CrossSections &cross_sections,
      DensitySubGridCreator< DensitySubGrid > &grid_creator,
      ThreadSafeVector< Task > &tasks)
      : _photon_source(photon_source), _buffers(buffers),
        _random_generators(random_generators),
        _discrete_photon_weight(discrete_photon_weight),
        _photon_source_spectrum(photon_source_spectrum),
        _abundances(abundances), _cross_sections(cross_sections),
        _grid_creator(grid_creator), _tasks(tasks) {}

  /**
   * @brief Execute a discrete photon source task.
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

    const size_t source_index = task.get_subgrid();

    const size_t num_photon_this_loop = task.get_buffer();
    const size_t subgrid_index = _photon_source.get_subgrid(source_index);

    // get a free photon buffer in the central queue
    uint_fast32_t buffer_index = _buffers.get_free_buffer();
    PhotonBuffer &input_buffer = _buffers[buffer_index];

    // set general buffer information
    input_buffer.grow(num_photon_this_loop);
    input_buffer.set_subgrid_index(subgrid_index);
    input_buffer.set_direction(TRAVELDIRECTION_INSIDE);

    const CoordinateVector<> source_position =
        _photon_source.get_position(source_index);

    // draw random photons and store them in the buffer
    for (uint_fast32_t i = 0; i < num_photon_this_loop; ++i) {

      PhotonPacket &photon = input_buffer[i];

      photon.set_type(PHOTONTYPE_PRIMARY);
      photon.set_scatter_counter(0);

      // initial position: we currently assume a single source at the
      // origin
      photon.set_position(source_position);

      // draw two pseudo random numbers
      const double cost =
          2. * _random_generators[thread_id].get_uniform_random_double() - 1.;
      const double phi =
          2. * M_PI * _random_generators[thread_id].get_uniform_random_double();

      // now use them to get all directional angles
      const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
      const double cosp = std::cos(phi);
      const double sinp = std::sin(phi);

      // set the direction...
      const CoordinateVector<> direction(sint * cosp, sint * sinp, cost);

      photon.set_direction(direction);

      // we currently assume equal weight for all photons
      photon.set_weight(_discrete_photon_weight);

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
    }

    // add to the queue of the corresponding thread
    DensitySubGrid &subgrid = *_grid_creator.get_subgrid(subgrid_index);
    const size_t task_index = _tasks.get_free_element();
    Task &new_task = _tasks[task_index];
    new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
    new_task.set_subgrid(subgrid_index);
    new_task.set_buffer(buffer_index);

    // add dependency for task:
    //  - subgrid
    // (the output buffers belong to the subgrid and do not count as a
    // dependency)
    new_task.set_dependency(subgrid.get_dependency());

    queues_to_add[0] = subgrid.get_owning_thread();
    tasks_to_add[0] = task_index;
    return 1;
  }
};

#endif // SOURCEDISCRETEPHOTONTASKCONTEXT_HPP
