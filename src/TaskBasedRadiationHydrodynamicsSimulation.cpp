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
 * @file TaskBasedRadiationHydrodynamicsSimulation.cpp
 *
 * @brief TaskBasedRadiationHydrodynamicsSimulation implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "TaskBasedRadiationHydrodynamicsSimulation.hpp"
#include "ChargeTransferRates.hpp"
#include "CommandLineParser.hpp"
#include "ContinuousPhotonSourceFactory.hpp"
#include "CrossSectionsFactory.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "HydroDensitySubGrid.hpp"
#include "LineCoolingData.hpp"
#include "MemorySpace.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "RecombinationRatesFactory.hpp"
#include "SimulationBox.hpp"
#include "TaskQueue.hpp"
#include "TemperatureCalculator.hpp"
#include "TimeLine.hpp"

#include <omp.h>

/**
 * @brief Make the hydro tasks for the given subgrid.
 *
 * @param tasks Task vector.
 * @param igrid Index of the subgrid.
 * @param grid_creator Subgrids.
 */
inline void
make_hydro_tasks(ThreadSafeVector< Task > &tasks, const uint_fast32_t igrid,
                 DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {

  HydroDensitySubGrid &this_grid = *grid_creator.get_subgrid(igrid);
  const unsigned int ngbx = this_grid.get_neighbour(TRAVELDIRECTION_FACE_X_P);
  const unsigned int ngby = this_grid.get_neighbour(TRAVELDIRECTION_FACE_Y_P);
  const unsigned int ngbz = this_grid.get_neighbour(TRAVELDIRECTION_FACE_Z_P);

  /// gradient computation and prediction

  // internal gradient sweep
  {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_INTERNAL);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    this_grid.set_hydro_task(0, next_task);
  }
  // external gradient sweeps
  // x
  // positive: always apply
  if (ngbx == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_X_P);
    this_grid.set_hydro_task(1, next_task);
  } else {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_NEIGHBOUR);
    // avoid dining philosophers by sorting the dependencies on subgrid index
    if (igrid < ngbx) {
      task.set_dependency(this_grid.get_dependency());
      task.set_extra_dependency(
          (*grid_creator.get_subgrid(ngbx)).get_dependency());
    } else {
      task.set_dependency((*grid_creator.get_subgrid(ngbx)).get_dependency());
      task.set_extra_dependency(this_grid.get_dependency());
    }
    task.set_subgrid(igrid);
    task.set_buffer(ngbx);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_X_P);
    this_grid.set_hydro_task(1, next_task);
  }
  // negative: only apply if we have a non-periodic boundary
  if (this_grid.get_neighbour(TRAVELDIRECTION_FACE_X_N) == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_X_N);
    this_grid.set_hydro_task(2, next_task);
  } else {
    this_grid.set_hydro_task(2, NO_TASK);
  }
  // y
  if (ngby == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_P);
    this_grid.set_hydro_task(3, next_task);
  } else {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_NEIGHBOUR);
    if (igrid < ngby) {
      task.set_dependency(this_grid.get_dependency());
      task.set_extra_dependency(
          (*grid_creator.get_subgrid(ngby)).get_dependency());
    } else {
      task.set_dependency((*grid_creator.get_subgrid(ngby)).get_dependency());
      task.set_extra_dependency(this_grid.get_dependency());
    }
    task.set_subgrid(igrid);
    task.set_buffer(ngby);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_P);
    this_grid.set_hydro_task(3, next_task);
  }
  if (this_grid.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_N);
    this_grid.set_hydro_task(4, next_task);
  } else {
    this_grid.set_hydro_task(4, NO_TASK);
  }
  // z
  if (ngbz == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_P);
    this_grid.set_hydro_task(5, next_task);
  } else {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_NEIGHBOUR);
    if (igrid < ngbz) {
      task.set_dependency(this_grid.get_dependency());
      task.set_extra_dependency(
          (*grid_creator.get_subgrid(ngbz)).get_dependency());
    } else {
      task.set_dependency((*grid_creator.get_subgrid(ngbz)).get_dependency());
      task.set_extra_dependency(this_grid.get_dependency());
    }
    task.set_subgrid(igrid);
    task.set_buffer(ngbz);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_P);
    this_grid.set_hydro_task(5, next_task);
  }
  if (this_grid.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_N);
    this_grid.set_hydro_task(6, next_task);
  } else {
    this_grid.set_hydro_task(6, NO_TASK);
  }
  // slope limiter
  {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_SLOPE_LIMITER);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    this_grid.set_hydro_task(7, next_task);
  }
  // primitive variable prediction
  {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_PREDICT_PRIMITIVES);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    this_grid.set_hydro_task(8, next_task);
  }

  /// flux exchange and primitive variable update

  // internal flux sweep
  {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_INTERNAL);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    this_grid.set_hydro_task(9, next_task);
  }
  // external flux sweeps
  // x
  // positive: always do flux exchange
  if (ngbx == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_X_P);
    this_grid.set_hydro_task(10, next_task);
  } else {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_NEIGHBOUR);
    if (igrid < ngbx) {
      task.set_dependency(this_grid.get_dependency());
      task.set_extra_dependency(
          (*grid_creator.get_subgrid(ngbx)).get_dependency());
    } else {
      task.set_dependency((*grid_creator.get_subgrid(ngbx)).get_dependency());
      task.set_extra_dependency(this_grid.get_dependency());
    }
    task.set_subgrid(igrid);
    task.set_buffer(ngbx);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_X_P);
    this_grid.set_hydro_task(10, next_task);
  }
  // negative: only do flux exchange when we need to apply a non-periodic
  // boundary condition
  if (this_grid.get_neighbour(TRAVELDIRECTION_FACE_X_N) == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_X_N);
    this_grid.set_hydro_task(11, next_task);
  } else {
    this_grid.set_hydro_task(11, NO_TASK);
  }
  // y
  if (ngby == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_P);
    this_grid.set_hydro_task(12, next_task);
  } else {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_NEIGHBOUR);
    if (igrid < ngby) {
      task.set_dependency(this_grid.get_dependency());
      task.set_extra_dependency(
          (*grid_creator.get_subgrid(ngby)).get_dependency());
    } else {
      task.set_dependency((*grid_creator.get_subgrid(ngby)).get_dependency());
      task.set_extra_dependency(this_grid.get_dependency());
    }
    task.set_subgrid(igrid);
    task.set_buffer(ngby);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_P);
    this_grid.set_hydro_task(12, next_task);
  }
  if (this_grid.get_neighbour(TRAVELDIRECTION_FACE_Y_N) == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Y_N);
    this_grid.set_hydro_task(13, next_task);
  } else {
    this_grid.set_hydro_task(13, NO_TASK);
  }
  // z
  if (ngbz == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_P);
    this_grid.set_hydro_task(14, next_task);
  } else {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_NEIGHBOUR);
    if (igrid < ngbz) {
      task.set_dependency(this_grid.get_dependency());
      task.set_extra_dependency(
          (*grid_creator.get_subgrid(ngbz)).get_dependency());
    } else {
      task.set_dependency((*grid_creator.get_subgrid(ngbz)).get_dependency());
      task.set_extra_dependency(this_grid.get_dependency());
    }
    task.set_subgrid(igrid);
    task.set_buffer(ngbz);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_P);
    this_grid.set_hydro_task(14, next_task);
  }
  if (this_grid.get_neighbour(TRAVELDIRECTION_FACE_Z_N) == NEIGHBOUR_OUTSIDE) {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    task.set_interaction_direction(TRAVELDIRECTION_FACE_Z_N);
    this_grid.set_hydro_task(15, next_task);
  } else {
    this_grid.set_hydro_task(15, NO_TASK);
  }
  // conserved variable update
  {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_UPDATE_CONSERVED);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    this_grid.set_hydro_task(16, next_task);
  }
  // primitive variable update
  {
    const size_t next_task = tasks.get_free_element();
    Task &task = tasks[next_task];
    task.set_type(TASKTYPE_UPDATE_PRIMITIVES);
    task.set_dependency(this_grid.get_dependency());
    task.set_subgrid(igrid);
    this_grid.set_hydro_task(17, next_task);
  }
}

/**
 * @brief Set the task dependencies for the given subgrid.
 *
 * @param igrid Subgrid index.
 * @param grid_creator Subgrids.
 * @param tasks Tasks.
 */
inline void
set_dependencies(const uint_fast32_t igrid,
                 DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
                 ThreadSafeVector< Task > &tasks) {

  const HydroDensitySubGrid &this_grid = *grid_creator.get_subgrid(igrid);

  // gradient sweeps unlock slope limiter tasks
  const size_t igg = this_grid.get_hydro_task(0);
  const size_t igxp = this_grid.get_hydro_task(1);
  size_t igxn = this_grid.get_hydro_task(2);
  if (igxn == NO_TASK) {
    igxn = (*grid_creator.get_subgrid(
                this_grid.get_neighbour(TRAVELDIRECTION_FACE_X_N)))
               .get_hydro_task(1);
  }
  const size_t igyp = this_grid.get_hydro_task(3);
  size_t igyn = this_grid.get_hydro_task(4);
  if (igyn == NO_TASK) {
    igyn = (*grid_creator.get_subgrid(
                this_grid.get_neighbour(TRAVELDIRECTION_FACE_Y_N)))
               .get_hydro_task(3);
  }
  const size_t igzp = this_grid.get_hydro_task(5);
  size_t igzn = this_grid.get_hydro_task(6);
  if (igzn == NO_TASK) {
    igzn = (*grid_creator.get_subgrid(
                this_grid.get_neighbour(TRAVELDIRECTION_FACE_Z_N)))
               .get_hydro_task(5);
  }

  const size_t isl = this_grid.get_hydro_task(7);
  tasks[igg].add_child(isl);
  tasks[igxp].add_child(isl);
  tasks[igxn].add_child(isl);
  tasks[igyp].add_child(isl);
  tasks[igyn].add_child(isl);
  tasks[igzp].add_child(isl);
  tasks[igzn].add_child(isl);

  // the slope limiter task unlocks the gradient prediction task
  const size_t ipp = this_grid.get_hydro_task(8);
  tasks[isl].add_child(ipp);

  // the gradient prediction task unlocks the flux task for this cell and all
  // neighbouring cell pair flux tasks
  const size_t iff = this_grid.get_hydro_task(9);
  tasks[ipp].add_child(iff);
  // neighbours: the positive one is (always) stored in this subgrid
  const size_t ifxp = this_grid.get_hydro_task(10);
  tasks[ipp].add_child(ifxp);
  // the negative one is only stored in this subgrid if it is a non-periodic
  // boundary
  size_t ifxn = this_grid.get_hydro_task(11);
  if (ifxn == NO_TASK) {
    ifxn = (*grid_creator.get_subgrid(
                this_grid.get_neighbour(TRAVELDIRECTION_FACE_X_N)))
               .get_hydro_task(10);
  }
  tasks[ipp].add_child(ifxn);
  const size_t ifyp = this_grid.get_hydro_task(12);
  tasks[ipp].add_child(ifyp);
  size_t ifyn = this_grid.get_hydro_task(13);
  if (ifyn == NO_TASK) {
    ifyn = (*grid_creator.get_subgrid(
                this_grid.get_neighbour(TRAVELDIRECTION_FACE_Y_N)))
               .get_hydro_task(12);
  }
  tasks[ipp].add_child(ifyn);
  const size_t ifzp = this_grid.get_hydro_task(14);
  tasks[ipp].add_child(ifzp);
  size_t ifzn = this_grid.get_hydro_task(15);
  if (ifzn == NO_TASK) {
    ifzn = (*grid_creator.get_subgrid(
                this_grid.get_neighbour(TRAVELDIRECTION_FACE_Z_N)))
               .get_hydro_task(14);
  }
  tasks[ipp].add_child(ifzn);

  // the flux tasks unlock the conserved variable update
  const size_t icu = this_grid.get_hydro_task(16);
  tasks[iff].add_child(icu);
  tasks[ifxp].add_child(icu);
  tasks[ifxn].add_child(icu);
  tasks[ifyp].add_child(icu);
  tasks[ifyn].add_child(icu);
  tasks[ifzp].add_child(icu);
  tasks[ifzn].add_child(icu);

  // the conserved variable update unlocks the primitive variable update
  const size_t ipu = this_grid.get_hydro_task(17);
  tasks[icu].add_child(ipu);
}

/**
 * @brief Reset the hydro tasks for the given subgrid.
 *
 * @param tasks Tasks.
 * @param this_grid Subgrid.
 */
inline void reset_hydro_tasks(ThreadSafeVector< Task > &tasks,
                              HydroDensitySubGrid &this_grid) {

  // gradient sweeps
  // internal
  tasks[this_grid.get_hydro_task(0)].set_number_of_unfinished_parents(0);
  // external
  tasks[this_grid.get_hydro_task(1)].set_number_of_unfinished_parents(0);
  if (this_grid.get_hydro_task(2) != NO_TASK) {
    tasks[this_grid.get_hydro_task(2)].set_number_of_unfinished_parents(0);
  }
  tasks[this_grid.get_hydro_task(3)].set_number_of_unfinished_parents(0);
  if (this_grid.get_hydro_task(4) != NO_TASK) {
    tasks[this_grid.get_hydro_task(4)].set_number_of_unfinished_parents(0);
  }
  tasks[this_grid.get_hydro_task(5)].set_number_of_unfinished_parents(0);
  if (this_grid.get_hydro_task(6) != NO_TASK) {
    tasks[this_grid.get_hydro_task(6)].set_number_of_unfinished_parents(0);
  }

  // slope limiter
  tasks[this_grid.get_hydro_task(7)].set_number_of_unfinished_parents(7);
  // primitive variable prediction
  tasks[this_grid.get_hydro_task(8)].set_number_of_unfinished_parents(1);

  // flux sweeps
  // internal
  tasks[this_grid.get_hydro_task(9)].set_number_of_unfinished_parents(1);
  // external
  if (tasks[this_grid.get_hydro_task(10)].get_type() ==
      TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY) {
    tasks[this_grid.get_hydro_task(10)].set_number_of_unfinished_parents(1);
  } else {
    tasks[this_grid.get_hydro_task(10)].set_number_of_unfinished_parents(2);
  }
  if (this_grid.get_hydro_task(11) != NO_TASK) {
    tasks[this_grid.get_hydro_task(11)].set_number_of_unfinished_parents(1);
  }
  if (tasks[this_grid.get_hydro_task(12)].get_type() ==
      TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY) {
    tasks[this_grid.get_hydro_task(12)].set_number_of_unfinished_parents(1);
  } else {
    tasks[this_grid.get_hydro_task(12)].set_number_of_unfinished_parents(2);
  }
  if (this_grid.get_hydro_task(13) != NO_TASK) {
    tasks[this_grid.get_hydro_task(13)].set_number_of_unfinished_parents(1);
  }
  if (tasks[this_grid.get_hydro_task(14)].get_type() ==
      TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY) {
    tasks[this_grid.get_hydro_task(14)].set_number_of_unfinished_parents(1);
  } else {
    tasks[this_grid.get_hydro_task(14)].set_number_of_unfinished_parents(2);
  }
  if (this_grid.get_hydro_task(15) != NO_TASK) {
    tasks[this_grid.get_hydro_task(15)].set_number_of_unfinished_parents(1);
  }

  // conserved variable update
  tasks[this_grid.get_hydro_task(16)].set_number_of_unfinished_parents(7);
  // primitive variable update
  tasks[this_grid.get_hydro_task(17)].set_number_of_unfinished_parents(1);
}

/**
 * @brief Steal a task from another queue.
 *
 * @param thread_id Id of the active thread.
 * @param num_threads Total number of threads.
 * @param queues Thread queues.
 * @param tasks Task space.
 * @param grid_creator Subgrids.
 * @return Index of an available task, or NO_TASK if no tasks are available.
 */
inline uint_fast32_t
steal_task(const int_fast32_t thread_id, const int_fast32_t num_threads,
           std::vector< TaskQueue * > &queues, ThreadSafeVector< Task > &tasks,
           DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {

  // sort the queues by size
  std::vector< uint_fast32_t > queue_sizes(queues.size(), 0);
  for (uint_fast32_t i = 0; i < queues.size(); ++i) {
    queue_sizes[i] = queues[i]->size();
  }
  std::vector< size_t > sorti = Utilities::argsort(queue_sizes);

  // now try to steal from the largest queue first
  uint_fast32_t current_index = NO_TASK;
  uint_fast32_t i = 0;
  while (current_index == NO_TASK && i < queue_sizes.size() &&
         queue_sizes[sorti[queue_sizes.size() - i - 1]] > 0) {
    current_index =
        queues[sorti[queue_sizes.size() - i - 1]]->try_get_task(tasks);
    ++i;
  }
  if (current_index != NO_TASK) {
    // stealing means transferring ownership...
    (*grid_creator.get_subgrid(tasks[current_index].get_subgrid()))
        .set_owning_thread(thread_id);
  }

  return current_index;
}

/**
 * @brief Execute a task.
 *
 * @param itask Task index.
 * @param grid_creator Subgrids.
 * @param tasks Tasks.
 * @param timestep System time step (in s).
 * @param hydro Hydro instance to use.
 * @param boundary HydroBoundary to use.
 */
inline void
execute_task(const size_t itask,
             DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
             ThreadSafeVector< Task > &tasks, const double timestep,
             const Hydro &hydro, const HydroBoundary &boundary) {

  const Task &task = tasks[itask];
  HydroDensitySubGrid &subgrid = *grid_creator.get_subgrid(task.get_subgrid());
  switch (task.get_type()) {
  case TASKTYPE_GRADIENTSWEEP_INTERNAL:
    subgrid.inner_gradient_sweep(hydro);
    break;
  case TASKTYPE_GRADIENTSWEEP_EXTERNAL_NEIGHBOUR:
    subgrid.outer_gradient_sweep(task.get_interaction_direction(), hydro,
                                 *grid_creator.get_subgrid(task.get_buffer()));
    break;
  case TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY:
    subgrid.outer_ghost_gradient_sweep(task.get_interaction_direction(), hydro,
                                       boundary);
    break;
  case TASKTYPE_SLOPE_LIMITER:
    subgrid.apply_slope_limiter(hydro);
    break;
  case TASKTYPE_PREDICT_PRIMITIVES:
    subgrid.predict_primitive_variables(hydro, 0.5 * timestep);
    break;
  case TASKTYPE_FLUXSWEEP_INTERNAL:
    subgrid.inner_flux_sweep(hydro);
    break;
  case TASKTYPE_FLUXSWEEP_EXTERNAL_NEIGHBOUR:
    subgrid.outer_flux_sweep(task.get_interaction_direction(), hydro,
                             *grid_creator.get_subgrid(task.get_buffer()));
    break;
  case TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY:
    subgrid.outer_ghost_flux_sweep(task.get_interaction_direction(), hydro,
                                   boundary);
    break;
  case TASKTYPE_UPDATE_CONSERVED:
    subgrid.update_conserved_variables(timestep);
    break;
  case TASKTYPE_UPDATE_PRIMITIVES:
    subgrid.update_primitive_variables(hydro);
    break;
  default:
    cmac_error("Unknown hydro task: %" PRIiFAST32, task.get_type());
  }
}

/**
 * @brief Add RHD simulation specific command line options to the command line
 * parser.
 *
 * The following parameters are added:
 *  - output-time-unit (no abbreviation, optional, string argument): unit used
 *    for time values that are written to the Log (default: s).
 *  - restart (no abbreviation, optional, string argument): restart a run from
 *    the restart file stored in the given folder (default: .).
 *
 * @param parser CommandLineParser that has not yet parsed the command line
 * options.
 */
void TaskBasedRadiationHydrodynamicsSimulation::add_command_line_parameters(
    CommandLineParser &parser) {
  parser.add_option("output-time-unit", 0,
                    "Unit for time stat output to the terminal.",
                    COMMANDLINEOPTION_STRINGARGUMENT, "s");
  parser.add_option("restart", 0,
                    "Restart from the restart file stored in the given folder.",
                    COMMANDLINEOPTION_STRINGARGUMENT, ".");
}

/**
 * @brief Perform an RHD simulation.
 *
 * This method reads the following parameters from the parameter file:
 *  - minimum timestep: Smallest possible hydrodynamical time step, an error is
 *    thrown if the time step wants to be smaller than this value (default:
 *    0.01*(total time))
 *  - maximum timestep: Largest possible hydrodynamical time step. The time step
 *    will always be smaller or equal to this value, even if the integration
 *    does not require it to be this small (default: 0.1*(total time))
 *  - total time: Total simulation time (default: 1. s)
 *  - snapshot time: Time interval between consecutive snapshot dumps (default:
 *    0.1*(total time))
 *  - radiation time: Time interval between consecutive updates of the
 *    ionization structure by running the photoionization code. If a negative
 *    value is given, the radiation field is updated every time step. (default:
 *    -1. s: update every time step)
 *  - random seed: Seed for the random number generator (default: 42)
 *  - output folder: Folder where all output files will be placed (default: .)
 *  - number of iterations: Number of iterations of the photoionization
 *    algorithm (default: 10)
 *  - number of photons: Number of photons to use during each iteration of the
 *    photoionization algorithm (default: 1e5)
 *  - number of photons first loop: Number of photons to use during the first
 *    iteration of the photoionization algorithm (default: (number of photons))
 *  - maximum neutral fraction: Maximum value of the hydrogen neutral fraction
 *    that is allowed at the start of a radiation step (negative values do not
 *    impose an upper limit, default: -1)
 *
 * @param parser CommandLineParser that contains the parsed command line
 * arguments.
 * @param write_output Flag indicating whether this process writes output.
 * @param programtimer Total program timer.
 * @param log Log to write logging info to.
 * @return Exit code: 0 on success.
 */
int TaskBasedRadiationHydrodynamicsSimulation::do_simulation(
    CommandLineParser &parser, bool write_output, Timer &programtimer,
    Log *log) {

  const int_fast32_t num_thread = parser.get_value< int_fast32_t >("threads");
  omp_set_num_threads(num_thread);

  // set the unit for time stats terminal output
  std::string output_time_unit =
      parser.get_value< std::string >("output-time-unit");

  // second: initialize the parameters that are read in from static files
  // these files should be configured by CMake and put in a location that is
  // stored in a CMake configured header
  LineCoolingData line_cooling_data;

  // third: read in the parameters of the run from a parameter file. This file
  // should be read by a ParameterFileParser object that acts as a dictionary
  ParameterFile *params =
      new ParameterFile(parser.get_value< std::string >("params"));

  // fourth: construct the density grid. This should be stored in a separate
  // DensityGrid object with geometrical and physical properties
  DensityFunction *density_function =
      DensityFunctionFactory::generate(*params, log);
  CrossSections *cross_sections = CrossSectionsFactory::generate(*params, log);
  RecombinationRates *recombination_rates =
      RecombinationRatesFactory::generate(*params, log);

  // initialize the simulation box
  const SimulationBox simulation_box(*params);

  const double hydro_total_time = params->get_physical_value< QUANTITY_TIME >(
      "TaskBasedRadiationHydrodynamicsSimulation:total time", "1. s");

  double hydro_minimum_timestep = params->get_physical_value< QUANTITY_TIME >(
      "TaskBasedRadiationHydrodynamicsSimulation:minimum timestep", "-1. s");
  if (hydro_minimum_timestep < 0.) {
    hydro_minimum_timestep = 1.e-10 * hydro_total_time;
  }

  double hydro_maximum_timestep = params->get_physical_value< QUANTITY_TIME >(
      "TaskBasedRadiationHydrodynamicsSimulation:maximum timestep", "-1. s");
  if (hydro_maximum_timestep < 0.) {
    hydro_maximum_timestep = 0.1 * hydro_total_time;
  }

  double hydro_snaptime = params->get_physical_value< QUANTITY_TIME >(
      "TaskBasedRadiationHydrodynamicsSimulation:snapshot time", "-1. s");
  if (hydro_snaptime < 0.) {
    hydro_snaptime = 0.1 * hydro_total_time;
  }
  uint_fast32_t hydro_lastsnap = 1;

  const double hydro_radtime = params->get_physical_value< QUANTITY_TIME >(
      "TaskBasedRadiationHydrodynamicsSimulation:radiation time", "-1. s");
  uint_fast32_t hydro_lastrad = 0;

  const double maximum_neutral_fraction = params->get_value< double >(
      "TaskBasedRadiationHydrodynamicsSimulation:maximum neutral fraction",
      -1.);

  const size_t number_of_buffers = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:number of buffers", 50000);
  MemorySpace *buffers = new MemorySpace(number_of_buffers);
  const size_t queue_size_per_thread = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:queue size per thread", 10000);
  std::vector< TaskQueue * > queues(num_thread);
  for (int_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    std::stringstream queue_name;
    queue_name << "Queue for Thread " << static_cast< int_fast32_t >(ithread);
    queues[ithread] = new TaskQueue(queue_size_per_thread, queue_name.str());
  }
  const size_t shared_queue_size = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:shared queue size", 100000);
  TaskQueue *shared_queue = new TaskQueue(shared_queue_size, "Shared queue");
  const size_t number_of_tasks = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:number of tasks", 500000);
  ThreadSafeVector< Task > *tasks =
      new ThreadSafeVector< Task >(number_of_tasks, "Tasks");
  std::vector< RandomGenerator > random_generators(num_thread);
  const int_fast32_t random_seed = params->get_value< int_fast32_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:random seed", 42);
  for (uint_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    random_generators[ithread].set_seed(random_seed + ithread);
  }
  DensitySubGridCreator< HydroDensitySubGrid > *grid_creator =
      new DensitySubGridCreator< HydroDensitySubGrid >(simulation_box.get_box(),
                                                       *params);

  // fifth: construct the stellar sources. These should be stored in a
  // separate StellarSources object with geometrical and physical properties.
  PhotonSourceDistribution *sourcedistribution =
      PhotonSourceDistributionFactory::generate(*params, log);
  PhotonSourceSpectrum *spectrum = PhotonSourceSpectrumFactory::generate(
      "PhotonSourceSpectrum", *params, log);

  if (sourcedistribution != nullptr && spectrum == nullptr) {
    cmac_error("No spectrum provided for the discrete photon sources!");
  }
  if (sourcedistribution == nullptr && spectrum != nullptr) {
    cmac_warning("Discrete photon source spectrum provided, but no discrete "
                 "photon source distributions. The given spectrum will be "
                 "ignored.");
  }

  ContinuousPhotonSource *continuoussource =
      ContinuousPhotonSourceFactory::generate(simulation_box.get_box(), *params,
                                              log);
  PhotonSourceSpectrum *continuousspectrum =
      PhotonSourceSpectrumFactory::generate("ContinuousPhotonSourceSpectrum",
                                            *params, log);

  if (continuoussource != nullptr && continuousspectrum == nullptr) {
    cmac_error("No spectrum provided for the continuous photon sources!");
  }
  if (continuoussource == nullptr && continuousspectrum != nullptr) {
    cmac_warning("Continuous photon source spectrum provided, but no "
                 "continuous photon source. The given spectrum will be "
                 "ignored.");
  }

  Abundances abundances(*params, log);

  PhotonSource source(sourcedistribution, spectrum, continuoussource,
                      continuousspectrum, abundances, *cross_sections, *params,
                      log);

  // set up output
  std::string output_folder =
      Utilities::get_absolute_path(params->get_value< std::string >(
          "TaskBasedRadiationHydrodynamicsSimulation:output folder", "."));
  DensityGridWriter *writer =
      DensityGridWriterFactory::generate(output_folder, *params, true, log);

  uint_fast32_t nloop = params->get_value< uint_fast32_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:number of iterations", 10);

  uint_fast64_t numphoton = params->get_value< uint_fast64_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:number of photons", 1e5);
  uint_fast64_t numphoton1 = params->get_value< uint_fast64_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:number of photons first loop",
      numphoton);

  const double CFL = params->get_value< double >(
      "TaskBasedRadiationHydrodynamicsSimulation:CFL", 0.2);

  Hydro hydro(params->get_value< double >(
      "TaskBasedRadiationHydrodynamicsSimulation:polytropic index", 5. / 3.));
  InflowHydroBoundary hydro_boundary;

  ChargeTransferRates charge_transfer_rates;

  // used to calculate both the ionization state and the temperature
  TemperatureCalculator *temperature_calculator = new TemperatureCalculator(
      source.get_total_luminosity(), abundances, line_cooling_data,
      *recombination_rates, charge_transfer_rates, *params, log);

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output) {
    std::ofstream pfile(output_folder + "/parameters-usedvalues.param");
    params->print_contents(pfile);
    pfile.close();
    if (log) {
      log->write_status("Wrote used parameters to ", output_folder,
                        "/parameters-usedvalues.param.");
    }
  }

  if (parser.get_value< bool >("dry-run")) {
    if (log) {
      log->write_warning("Dry run requested. Program will now halt.");
    }
    return 0;
  }

  if (log) {
    log->write_status("Initializing DensityFunction...");
  }
  density_function->initialize();
  if (log) {
    log->write_status("Done.");
  }

  if (log) {
    log->write_status("Initializing grid...");
  }
  grid_creator->initialize(*density_function);
  double requested_timestep = DBL_MAX;
  for (auto cellit = grid_creator->begin();
       cellit != grid_creator->original_end(); ++cellit) {
    requested_timestep =
        std::min(requested_timestep,
                 (*cellit).initialize_hydrodynamic_variables(hydro, true));
    make_hydro_tasks(*tasks, cellit.get_index(), *grid_creator);
  }
  for (auto cellit = grid_creator->begin();
       cellit != grid_creator->original_end(); ++cellit) {
    set_dependencies(cellit.get_index(), *grid_creator, *tasks);
  }
  if (log) {
    log->write_status("Done.");
  }

  if (write_output) {
    writer->write(*grid_creator, 0, *params, 0.);
  }

  double maximum_timestep = hydro_maximum_timestep;
  if (hydro_radtime > 0.) {
    // make sure the system is evolved hydrodynamically in between successive
    // ionization steps
    maximum_timestep = std::min(maximum_timestep, hydro_radtime);
  }

  TimeLine *timeline = new TimeLine(
      0., hydro_total_time, hydro_minimum_timestep, maximum_timestep, log);
  uint_fast32_t num_step = 0;
  double actual_timestep, current_time;
  requested_timestep *= CFL;
  bool has_next_step =
      timeline->advance(requested_timestep, actual_timestep, current_time);
  bool stop_simulation = false;
  while (has_next_step && !stop_simulation) {

    if (log) {
      log->write_status("Starting hydro integration step ", num_step,
                        ", t = ", current_time);
    }
    ++num_step;

    // reset the hydro tasks and add them to the queue
    AtomicValue< uint_fast32_t > number_of_tasks;
    for (auto cellit = grid_creator->begin();
         cellit != grid_creator->original_end(); ++cellit) {
      reset_hydro_tasks(*tasks, *cellit);
      for (int_fast8_t i = 0; i < 18; ++i) {
        const size_t itask = (*cellit).get_hydro_task(i);
        if (itask != NO_TASK &&
            (*tasks)[itask].get_number_of_unfinished_parents() == 0) {
          queues[(*cellit).get_owning_thread()]->add_task(itask);
          number_of_tasks.pre_increment();
        }
      }
    }

#pragma omp parallel default(shared)
    {
      const int_fast32_t thread_id = omp_get_thread_num();
      while (number_of_tasks.value() > 0) {
        size_t current_task = queues[thread_id]->get_task(*tasks);
        if (current_task == NO_TASK) {
          current_task =
              steal_task(thread_id, num_thread, queues, *tasks, *grid_creator);
        }
        if (current_task != NO_TASK) {
          (*tasks)[current_task].start(thread_id);
          execute_task(current_task, *grid_creator, *tasks, actual_timestep,
                       hydro, hydro_boundary);
          (*tasks)[current_task].stop();
          (*tasks)[current_task].unlock_dependency();
          const unsigned char numchild =
              (*tasks)[current_task].get_number_of_children();
          for (uint_fast8_t i = 0; i < numchild; ++i) {
            const size_t ichild = (*tasks)[current_task].get_child(i);
            if ((*tasks)[ichild].decrement_number_of_unfinished_parents() ==
                0) {
              queues[(*grid_creator->get_subgrid(
                          (*tasks)[ichild].get_subgrid()))
                         .get_owning_thread()]
                  ->add_task(ichild);
              number_of_tasks.pre_increment();
            }
          }
          number_of_tasks.pre_decrement();
        }
      }
    }

    requested_timestep = DBL_MAX;
    for (auto gridit = grid_creator->begin();
         gridit != grid_creator->original_end(); ++gridit) {
      for (auto cellit = (*gridit).hydro_begin();
           cellit != (*gridit).hydro_end(); ++cellit) {
        requested_timestep = std::min(
            requested_timestep, hydro.get_timestep(cellit.get_hydro_variables(),
                                                   cellit.get_volume()));
      }
    }
    requested_timestep *= CFL;
    has_next_step =
        timeline->advance(requested_timestep, actual_timestep, current_time);
  }

  if (write_output) {
    writer->write(*grid_creator, hydro_lastsnap, *params, current_time);
  }

  (void)num_step;
  (void)hydro_lastsnap;
  (void)hydro_lastrad;
  (void)maximum_neutral_fraction;
  (void)nloop;
  (void)numphoton1;
  if (sourcedistribution != nullptr) {
    delete sourcedistribution;
  }
  if (continuoussource != nullptr) {
    delete continuoussource;
  }
  delete density_function;
  delete writer;
  delete temperature_calculator;
  delete continuousspectrum;
  delete spectrum;

  delete cross_sections;
  delete recombination_rates;

  delete params;

  delete buffers;
  for (int_fast32_t ithread = 0; ithread < num_thread; ++ithread) {
    delete queues[ithread];
  }
  delete shared_queue;
  delete tasks;
  delete grid_creator;

  return 0;
}
