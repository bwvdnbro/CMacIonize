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
#include "DiffuseReemissionHandler.hpp"
#include "DistributedPhotonSource.hpp"
#include "HydroBoundaryManager.hpp"
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

/*! @brief Stop the serial time timer and start the parallel time timer. */
#define start_parallel_timing_block()                                          \
  serial_timer.stop();                                                         \
  parallel_timer.start();

/*! @brief Stop the parallel time timer and start the serial time timer. */
#define stop_parallel_timing_block()                                           \
  parallel_timer.stop();                                                       \
  serial_timer.start();

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
 * @param boundary_manager HydroBoundaryManager to use.
 */
inline void
execute_task(const size_t itask,
             DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
             ThreadSafeVector< Task > &tasks, const double timestep,
             const Hydro &hydro, const HydroBoundaryManager &boundary_manager) {

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
                                       boundary_manager.get_boundary_condition(
                                           task.get_interaction_direction()));
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
                                   boundary_manager.get_boundary_condition(
                                       task.get_interaction_direction()));
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
 *  - task-plot (no abbreviation, optional, integer argument): output task
 *    information for the first N steps of the algorithm (default: 0).
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
  parser.add_option(
      "task-plot", 0,
      "Output task information for the first N steps of the algorithm",
      COMMANDLINEOPTION_INTARGUMENT, "0");
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
 *  - diffuse field: Enable diffuse reemission? (default: no)
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

  Timer total_timer;
  Timer serial_timer;
  Timer parallel_timer;
  Timer worktimer;

  total_timer.start();
  serial_timer.start();

  const int_fast32_t num_thread = parser.get_value< int_fast32_t >("threads");
  omp_set_num_threads(num_thread);

  // set the unit for time stats terminal output
  std::string output_time_unit =
      parser.get_value< std::string >("output-time-unit");

  const int_fast32_t task_plot_N =
      parser.get_value< int_fast32_t >("task-plot");

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
  DiffuseReemissionHandler *reemission_handler = nullptr;
  if (params->get_value< bool >(
          "TaskBasedRadiationHydrodynamicsSimulation:diffuse field", false)) {
    reemission_handler = new DiffuseReemissionHandler(*cross_sections);
  }

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
  int_fast32_t task_plot_i = 0;

  const double hydro_radtime = params->get_physical_value< QUANTITY_TIME >(
      "TaskBasedRadiationHydrodynamicsSimulation:radiation time", "-1. s");
  uint_fast32_t hydro_lastrad = 0;

  const double maximum_neutral_fraction = params->get_value< double >(
      "TaskBasedRadiationHydrodynamicsSimulation:maximum neutral fraction",
      -1.);

  const size_t number_of_buffers = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:number of buffers", 50000);
  const size_t queue_size_per_thread = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:queue size per thread", 10000);
  const size_t shared_queue_size = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:shared queue size", 100000);
  const size_t number_of_tasks = params->get_value< size_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:number of tasks", 500000);
  const int_fast32_t random_seed = params->get_value< int_fast32_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:random seed", 42);
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
      "TaskBasedRadiationHydrodynamicsSimulation:number of photons", 1e6);

  const double CFL = params->get_value< double >(
      "TaskBasedRadiationHydrodynamicsSimulation:CFL", 0.2);

  const double source_copy_level = params->get_value< double >(
      "TaskBasedRadiationHydrodynamicsSimulation:source copy level", 4);

  Hydro hydro(*params);
  HydroBoundaryManager hydro_boundary_manager(*params);

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
  start_parallel_timing_block();
  grid_creator->initialize(*density_function);
  stop_parallel_timing_block();
  if (log) {
    log->write_status("Done.");
  }

  density_function->free();

  if (log) {
    log->write_status("Initializing task-based structures...");
  }
  if (log) {
    log->write_status("Allocating memory space...");
  }
  MemorySpace *buffers = new MemorySpace(number_of_buffers);
  if (log) {
    log->write_status("Allocating task memory...");
  }
  ThreadSafeVector< Task > *tasks =
      new ThreadSafeVector< Task >(number_of_tasks, "Tasks");
  if (log) {
    log->write_status("Allocating shared queue...");
  }
  TaskQueue *shared_queue = new TaskQueue(shared_queue_size, "Shared queue");
  if (log) {
    log->write_status("Allocating per thread queues...");
  }
  std::vector< TaskQueue * > queues(num_thread);
  for (int_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    std::stringstream queue_name;
    queue_name << "Queue for Thread " << static_cast< int_fast32_t >(ithread);
    queues[ithread] = new TaskQueue(queue_size_per_thread, queue_name.str());
  }
  if (log) {
    log->write_status("Initializing per thread random generators...");
  }
  std::vector< RandomGenerator > random_generators(num_thread);
  for (uint_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    random_generators[ithread].set_seed(random_seed + ithread);
  }
  if (log) {
    log->write_status("Done initializing task-based structures.");
  }

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
  const size_t radiation_task_offset = tasks->size();

  if (write_output) {
    writer->write(*grid_creator, 0, *params, 0.);
  }

  double maximum_timestep = hydro_maximum_timestep;
  if (hydro_radtime > 0.) {
    // make sure the system is evolved hydrodynamically in between successive
    // ionization steps
    maximum_timestep = std::min(maximum_timestep, hydro_radtime);
  }

  /// RADIATION
  std::vector< uint_fast8_t > levels(
      grid_creator->number_of_original_subgrids(), 0);

  // set the copy level off all subgrids containing a source to the given
  // parameter value (for now)
  {
    const photonsourcenumber_t number_of_sources =
        sourcedistribution->get_number_of_sources();
    for (photonsourcenumber_t isource = 0; isource < number_of_sources;
         ++isource) {
      const CoordinateVector<> position =
          sourcedistribution->get_position(isource);
      DensitySubGridCreator< HydroDensitySubGrid >::iterator gridit =
          grid_creator->get_subgrid(position);
      levels[gridit.get_index()] = source_copy_level;
    }
  }

  // impose copy restrictions
  {
    uint_fast8_t max_level = 0;
    const size_t levelsize = levels.size();
    for (size_t i = 0; i < levelsize; ++i) {
      max_level = std::max(max_level, levels[i]);
    }

    size_t ngbs[6];
    while (max_level > 0) {
      for (size_t i = 0; i < levelsize; ++i) {
        if (levels[i] == max_level) {
          const uint_fast8_t numngbs = grid_creator->get_neighbours(i, ngbs);
          for (uint_fast8_t ingb = 0; ingb < numngbs; ++ingb) {
            const size_t ngbi = ngbs[ingb];
            if (levels[ngbi] < levels[i] - 1) {
              levels[ngbi] = levels[i] - 1;
            }
          }
        }
      }
      --max_level;
    }
  }
  grid_creator->create_copies(levels);

  /// SIMULATION
  TimeLine *timeline = new TimeLine(
      0., hydro_total_time, hydro_minimum_timestep, maximum_timestep, log);
  uint_fast32_t num_step = 0;
  double actual_timestep, current_time;
  requested_timestep *= CFL;
  bool has_next_step =
      timeline->advance(requested_timestep, actual_timestep, current_time);
  bool stop_simulation = false;
  while (has_next_step && !stop_simulation) {

    uint_fast64_t iteration_start, iteration_end;
    cpucycle_tick(iteration_start);
    std::vector< uint_fast64_t > active_time(num_thread, 0);

    ++num_step;

    if (log) {
      log->write_status("Starting hydro step ", num_step, ", t = ",
                        UnitConverter::to_unit_string< QUANTITY_TIME >(
                            current_time - actual_timestep, output_time_unit),
                        ", dt = ",
                        UnitConverter::to_unit_string< QUANTITY_TIME >(
                            actual_timestep, output_time_unit),
                        ".");
    }

    // decide whether or not to do the radiation step
    if (hydro_radtime < 0. ||
        (current_time - actual_timestep) >= hydro_lastrad * hydro_radtime) {

      if (log) {
        log->write_status("Starting radiation step...");
      }

      ++hydro_lastrad;
      // update the PhotonSource
      if (sourcedistribution->update(current_time)) {
        source.update(sourcedistribution);
        temperature_calculator->update_luminosity(
            source.get_total_luminosity());
      }

      if (source.get_total_luminosity() > 0.) {
        {
          AtomicValue< size_t > igrid(0);
          start_parallel_timing_block();
#pragma omp parallel default(shared)
          while (igrid.value() < grid_creator->number_of_original_subgrids()) {
            const size_t this_igrid = igrid.post_increment();
            if (this_igrid < grid_creator->number_of_original_subgrids()) {
              HydroDensitySubGrid &subgrid =
                  *grid_creator->get_subgrid(this_igrid);
              subgrid.update_ionization_variables(hydro,
                                                  maximum_neutral_fraction);
            }
          }
          stop_parallel_timing_block();
        }

        DistributedPhotonSource< HydroDensitySubGrid > photon_source(
            numphoton, *sourcedistribution, *grid_creator);
        {
          AtomicValue< size_t > igrid(0);
          start_parallel_timing_block();
#pragma omp parallel default(shared)
          while (igrid.value() < grid_creator->number_of_actual_subgrids()) {
            const size_t this_igrid = igrid.post_increment();
            if (this_igrid < grid_creator->number_of_actual_subgrids()) {
              HydroDensitySubGrid &subgrid =
                  *grid_creator->get_subgrid(this_igrid);
              for (int ingb = 0; ingb < TRAVELDIRECTION_NUMBER; ++ingb) {
                subgrid.set_active_buffer(ingb, NEIGHBOUR_OUTSIDE);
              }
            }
          }
          stop_parallel_timing_block();
        }

        for (uint_fast32_t iloop = 0; iloop < nloop; ++iloop) {

          if (log) {
            log->write_status("Starting loop ", iloop, ".");
          }

          worktimer.start();

          // update copies
          grid_creator->update_copy_properties();

          // reset the photon source information
          photon_source.reset();

          // reset the diffuse field variables
          if (reemission_handler != nullptr) {
            AtomicValue< size_t > igrid(0);
#pragma omp parallel default(shared)
            while (igrid.value() <
                   grid_creator->number_of_original_subgrids()) {
              const size_t this_igrid = igrid.post_increment();
              if (this_igrid < grid_creator->number_of_original_subgrids()) {
                auto gridit = grid_creator->get_subgrid(this_igrid);
                for (auto cellit = (*gridit).begin(); cellit != (*gridit).end();
                     ++cellit) {
                  IonizationVariables &vars = cellit.get_ionization_variables();
                  DiffuseReemissionHandler::set_reemission_probabilities(vars);
                }
              }
            }
          }

          for (size_t ibatch = 0; ibatch < photon_source.get_number_of_batches(
                                               0, PHOTONBUFFER_SIZE);
               ++ibatch) {
            for (size_t isrc = 0; isrc < photon_source.get_number_of_sources();
                 ++isrc) {
              const size_t new_task = tasks->get_free_element();
              (*tasks)[new_task].set_type(TASKTYPE_SOURCE_PHOTON);
              (*tasks)[new_task].set_subgrid(isrc);
              (*tasks)[new_task].set_buffer(
                  photon_source.get_photon_batch(isrc, PHOTONBUFFER_SIZE));
              shared_queue->add_task(new_task);
            }
          }

          bool global_run_flag = true;
          const uint_fast32_t num_empty_target =
              TRAVELDIRECTION_NUMBER *
              grid_creator->number_of_actual_subgrids();
          AtomicValue< uint_fast32_t > num_empty(num_empty_target);
          AtomicValue< uint_fast32_t > num_active_buffers(0);
          AtomicValue< uint_fast32_t > num_photon_done(0);
          start_parallel_timing_block();
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
            uint_fast32_t current_index = shared_queue->get_task(*tasks);
            while (global_run_flag) {

              if (current_index == NO_TASK) {
                uint_fast32_t threshold_size = PHOTONBUFFER_SIZE;
                while (threshold_size > 0) {
                  threshold_size >>= 1;
                  for (auto gridit = grid_creator->begin();
                       gridit != grid_creator->all_end(); ++gridit) {
                    DensitySubGrid &this_subgrid = *gridit;
                    if (this_subgrid.get_largest_buffer_size() >
                            threshold_size &&
                        this_subgrid.get_dependency()->try_lock()) {

                      const uint_fast8_t largest_index =
                          this_subgrid.get_largest_buffer_index();
                      if (largest_index != TRAVELDIRECTION_NUMBER) {

                        const uint_fast32_t non_full_index =
                            this_subgrid.get_active_buffer(largest_index);
                        this_subgrid.set_active_buffer(largest_index,
                                                       NEIGHBOUR_OUTSIDE);
                        // we are creating a new active photon buffer
                        num_active_buffers.pre_increment();
                        // we created a new empty buffer
                        num_empty.pre_increment();

                        const size_t task_index = tasks->get_free_element();
                        Task &new_task = (*tasks)[task_index];
                        new_task.set_subgrid(
                            (*buffers)[non_full_index].get_subgrid_index());
                        new_task.set_buffer(non_full_index);
                        if (largest_index > 0) {
                          DensitySubGrid &subgrid = *grid_creator->get_subgrid(
                              (*buffers)[non_full_index].get_subgrid_index());
                          new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);

                          // add dependency
                          new_task.set_dependency(subgrid.get_dependency());

                          const uint_fast32_t queue_index =
                              subgrid.get_owning_thread();
                          queues[queue_index]->add_task(task_index);
                        } else {
                          new_task.set_type(TASKTYPE_PHOTON_REEMIT);
                          // a reemit task has no dependencies
                          shared_queue->add_task(task_index);
                        }

                        // set the new largest index
                        uint_fast8_t new_largest_index = TRAVELDIRECTION_NUMBER;
                        uint_fast32_t new_largest_size = 0;
                        for (uint_fast8_t ibuffer = 0;
                             ibuffer < TRAVELDIRECTION_NUMBER; ++ibuffer) {
                          if (this_subgrid.get_active_buffer(ibuffer) !=
                                  NEIGHBOUR_OUTSIDE &&
                              (*buffers)[this_subgrid.get_active_buffer(
                                             ibuffer)]
                                      .size() > new_largest_size) {
                            new_largest_index = ibuffer;
                            new_largest_size =
                                (*buffers)[this_subgrid.get_active_buffer(
                                               ibuffer)]
                                    .size();
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
                        // no semi-full buffers for this subgrid: release the
                        // lock again
                        this_subgrid.get_dependency()->unlock();
                      }
                    }
                  }
                }
                current_index = queues[thread_id]->get_task(*tasks);
                if (current_index == NO_TASK) {
                  current_index = steal_task(thread_id, num_thread, queues,
                                             *tasks, *grid_creator);
                  if (current_index == NO_TASK) {
                    current_index = shared_queue->get_task(*tasks);
                  }
                }
              }

              while (current_index != NO_TASK) {

                // execute task
                uint_fast32_t num_tasks_to_add = 0;
                uint_fast32_t tasks_to_add[TRAVELDIRECTION_NUMBER];
                int_fast32_t queues_to_add[TRAVELDIRECTION_NUMBER];

                Task &task = (*tasks)[current_index];
                uint_fast64_t task_start, task_stop;
                cpucycle_tick(task_start);

                if (task.get_type() == TASKTYPE_SOURCE_PHOTON) {

                  task.start(thread_id);
                  num_active_buffers.pre_increment();
                  const size_t source_index = task.get_subgrid();

                  const size_t num_photon_this_loop = task.get_buffer();
                  const size_t subgrid_index =
                      photon_source.get_subgrid(source_index);

                  // get a free photon buffer in the central queue
                  uint_fast32_t buffer_index = (*buffers).get_free_buffer();
                  PhotonBuffer &input_buffer = (*buffers)[buffer_index];

                  // set general buffer information
                  input_buffer.grow(num_photon_this_loop);
                  input_buffer.set_subgrid_index(subgrid_index);
                  input_buffer.set_direction(TRAVELDIRECTION_INSIDE);

                  const CoordinateVector<> source_position =
                      photon_source.get_position(source_index);

                  // draw random photons and store them in the buffer
                  for (uint_fast32_t i = 0; i < num_photon_this_loop; ++i) {

                    PhotonPacket &photon = input_buffer[i];

                    // initial position: we currently assume a single source at
                    // the origin
                    photon.set_position(source_position);

                    // draw two pseudo random numbers
                    const double cost = 2. * random_generators[thread_id]
                                                 .get_uniform_random_double() -
                                        1.;
                    const double phi = 2. * M_PI *
                                       random_generators[thread_id]
                                           .get_uniform_random_double();

                    // now use them to get all directional angles
                    const double sint =
                        std::sqrt(std::max(1. - cost * cost, 0.));
                    const double cosp = std::cos(phi);
                    const double sinp = std::sin(phi);

                    // set the direction...
                    const CoordinateVector<> direction(sint * cosp, sint * sinp,
                                                       cost);

                    photon.set_direction(direction);

                    // we currently assume equal weight for all photons
                    photon.set_weight(1.);

                    // target optical depth (exponential distribution)
                    photon.set_target_optical_depth(
                        -std::log(random_generators[thread_id]
                                      .get_uniform_random_double()));

                    const double frequency = spectrum->get_random_frequency(
                        random_generators[thread_id]);
                    photon.set_energy(frequency);
                    for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES;
                         ++ion) {
                      double sigma =
                          cross_sections->get_cross_section(ion, frequency);
                      if (ion != ION_H_n) {
                        sigma *= abundances.get_abundance(get_element(ion));
                      }
                      // this is the fixed cross section we use for the moment
                      photon.set_photoionization_cross_section(ion, sigma);
                    }
                  }

                  // add to the queue of the corresponding thread
                  DensitySubGrid &subgrid =
                      *grid_creator->get_subgrid(subgrid_index);
                  const size_t task_index = tasks->get_free_element();
                  Task &new_task = (*tasks)[task_index];
                  new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
                  new_task.set_subgrid(subgrid_index);
                  new_task.set_buffer(buffer_index);

                  // add dependency for task:
                  //  - subgrid
                  // (the output buffers belong to the subgrid and do not count
                  // as a dependency)
                  new_task.set_dependency(subgrid.get_dependency());

                  queues_to_add[num_tasks_to_add] = subgrid.get_owning_thread();
                  tasks_to_add[num_tasks_to_add] = task_index;
                  ++num_tasks_to_add;

                  // log the end time of the task
                  task.stop();

                } else if (task.get_type() == TASKTYPE_PHOTON_REEMIT) {

                  task.start(thread_id);

                  const size_t current_buffer_index = task.get_buffer();
                  PhotonBuffer &buffer = (*buffers)[current_buffer_index];

                  uint_fast32_t num_photon_done_now = buffer.size();
                  DensitySubGrid &subgrid =
                      *grid_creator->get_subgrid(task.get_subgrid());

                  // reemission
                  uint_fast32_t index = 0;
                  for (uint_fast32_t iphoton = 0; iphoton < buffer.size();
                       ++iphoton) {
                    PhotonPacket &old_photon = buffer[iphoton];
                    const IonizationVariables &ionization_variables =
                        subgrid.get_cell(old_photon.get_position())
                            .get_ionization_variables();
                    const double AHe = 0.;
                    const double new_frequency = reemission_handler->reemit(
                        old_photon, AHe, ionization_variables,
                        random_generators[thread_id]);
                    if (new_frequency > 0.) {
                      PhotonPacket &new_photon = buffer[index];
                      new_photon.set_position(old_photon.get_position());
                      new_photon.set_weight(old_photon.get_weight());

                      new_photon.set_energy(new_frequency);
                      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES;
                           ++ion) {
                        double sigma = cross_sections->get_cross_section(
                            ion, new_frequency);
                        // this is the fixed cross section we use for the moment
                        new_photon.set_photoionization_cross_section(ion,
                                                                     sigma);
                      }

                      // draw two pseudo random numbers
                      const double cost =
                          2. * random_generators[thread_id]
                                   .get_uniform_random_double() -
                          1.;
                      const double phi = 2. * M_PI *
                                         random_generators[thread_id]
                                             .get_uniform_random_double();

                      // now use them to get all directional angles
                      const double sint =
                          std::sqrt(std::max(1. - cost * cost, 0.));
                      const double cosp = std::cos(phi);
                      const double sinp = std::sin(phi);

                      // set the direction...
                      const CoordinateVector<> direction(sint * cosp,
                                                         sint * sinp, cost);

                      new_photon.set_direction(direction);

                      // target optical depth (exponential distribution)
                      new_photon.set_target_optical_depth(
                          -std::log(random_generators[thread_id]
                                        .get_uniform_random_double()));

                      ++index;
                    }
                  }
                  // update the size of the buffer to account for photons that
                  // were not reemitted
                  buffer.grow(index);

                  num_photon_done_now -= buffer.size();
                  num_photon_done.pre_add(num_photon_done_now);

                  const size_t task_index = tasks->get_free_element();
                  Task &new_task = (*tasks)[task_index];
                  new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
                  new_task.set_subgrid(task.get_subgrid());
                  new_task.set_buffer(current_buffer_index);
                  new_task.set_dependency(subgrid.get_dependency());

                  queues_to_add[num_tasks_to_add] = subgrid.get_owning_thread();
                  tasks_to_add[num_tasks_to_add] = task_index;
                  ++num_tasks_to_add;

                  task.stop();

                } else if (task.get_type() == TASKTYPE_PHOTON_TRAVERSAL) {

                  task.start(thread_id);

                  const uint_fast32_t current_buffer_index = task.get_buffer();
                  PhotonBuffer &photon_buffer =
                      (*buffers)[current_buffer_index];
                  const uint_fast32_t igrid = photon_buffer.get_subgrid_index();
                  DensitySubGrid &this_grid = *grid_creator->get_subgrid(igrid);

                  // prepare output buffers: make sure they are empty and that
                  // buffers corresponding to directions outside the simulation
                  // box are disabled
                  for (int_fast8_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
                    const uint_fast32_t ngb = this_grid.get_neighbour(i);
                    if (ngb != NEIGHBOUR_OUTSIDE) {
                      local_buffer_flags[i] = true;
                      local_buffers[i].reset();
                    } else {
                      local_buffer_flags[i] = false;
                    }
                  }

                  // only enable output to internal buffer if reemission
                  // is enabled
                  if (reemission_handler == nullptr) {
                    local_buffer_flags[TRAVELDIRECTION_INSIDE] = false;
                  }

                  // keep track of the original number of photons
                  uint_fast32_t num_photon_done_now = photon_buffer.size();

                  // now loop over the input buffer photons and traverse them
                  // one by one
                  for (uint_fast32_t i = 0; i < photon_buffer.size(); ++i) {

                    // active photon
                    PhotonPacket &photon = photon_buffer[i];

                    // make sure the photon is moving in *a* direction
                    cmac_assert_message(photon.get_direction()[0] != 0. ||
                                            photon.get_direction()[1] != 0. ||
                                            photon.get_direction()[2] != 0.,
                                        "size: %" PRIuFAST32,
                                        photon_buffer.size());

                    // traverse the photon through the active subgrid
                    const int_fast32_t result = this_grid.interact(
                        photon, photon_buffer.get_direction());

                    // check that the photon ended up in a valid output buffer
                    cmac_assert_message(
                        result >= 0 && result < TRAVELDIRECTION_NUMBER, "fail");

                    // add the photon to an output buffer, if it still exists
                    // (if the corresponding output buffer does not exist, this
                    // means the photon left the simulation box)
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

                      // photon packets that are still present in an output
                      // buffer are not done yet
                      num_photon_done_now -= local_buffers[i].size();

                      // move photon packets from the local temporary buffer
                      // (that is guaranteed to be large enough) to the actual
                      // output buffer for that direction (which might cause on
                      // overflow)
                      const uint_fast32_t ngb = this_grid.get_neighbour(i);
                      uint_fast32_t new_index = this_grid.get_active_buffer(i);

                      if (new_index == NEIGHBOUR_OUTSIDE) {
                        // buffer was not created yet: create it now
                        new_index = buffers->get_free_buffer();
                        PhotonBuffer &buffer = (*buffers)[new_index];
                        buffer.set_subgrid_index(ngb);
                        buffer.set_direction(
                            TravelDirections::output_to_input_direction(i));
                        this_grid.set_active_buffer(i, new_index);
                      }

                      if ((*buffers)[new_index].size() == 0) {
                        // we are adding photons to an empty buffer
                        num_empty.pre_decrement();
                      }
                      uint_fast32_t add_index =
                          buffers->add_photons(new_index, local_buffers[i]);

                      // check if the original buffer is full
                      if (add_index != new_index) {

                        // a new active buffer was created
                        num_active_buffers.pre_increment();

                        // new_buffers.add_photons already created a new empty
                        // buffer, set it as the active buffer for this output
                        // direction
                        if ((*buffers)[add_index].size() == 0) {
                          buffers->free_buffer(add_index);
                          this_grid.set_active_buffer(i, NEIGHBOUR_OUTSIDE);
                          // we have created a new empty buffer
                          num_empty.pre_increment();
                        } else {
                          this_grid.set_active_buffer(i, add_index);

                          cmac_assert_message(
                              (*buffers)[add_index].get_subgrid_index() == ngb,
                              "Wrong subgrid");
                          cmac_assert_message(
                              (*buffers)[add_index].get_direction() ==
                                  TravelDirections::output_to_input_direction(
                                      i),
                              "Wrong direction");
                        }

                        // YES: create a task for the buffer and add it to the
                        // queue the task type depends on the buffer: photon
                        // packets in the internal buffer were absorbed and
                        // could be reemitted, photon packets in the other
                        // buffers left the subgrid and need to be traversed in
                        // the neighbouring subgrid
                        if (i > 0) {
                          DensitySubGrid &subgrid = *grid_creator->get_subgrid(
                              (*buffers)[new_index].get_subgrid_index());
                          const size_t task_index = tasks->get_free_element();
                          Task &new_task = (*tasks)[task_index];
                          new_task.set_subgrid(
                              (*buffers)[new_index].get_subgrid_index());
                          new_task.set_buffer(new_index);
                          new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);

                          // add dependencies for task:
                          //  - subgrid
                          new_task.set_dependency(subgrid.get_dependency());

                          // add the task to the queue of the corresponding
                          // thread
                          const int_fast32_t queue_index =
                              (*grid_creator->get_subgrid(ngb))
                                  .get_owning_thread();
                          queues_to_add[num_tasks_to_add] = queue_index;
                          tasks_to_add[num_tasks_to_add] = task_index;
                          ++num_tasks_to_add;
                        } else {
                          const size_t task_index = tasks->get_free_element();
                          Task &new_task = (*tasks)[task_index];
                          new_task.set_subgrid(
                              (*buffers)[new_index].get_subgrid_index());
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

                    // we have to do this outside the other condition, as
                    // buffers to which nothing was added can still be
                    // non-empty...
                    if (local_buffer_flags[i]) {
                      uint_fast32_t new_index = this_grid.get_active_buffer(i);
                      if (new_index != NEIGHBOUR_OUTSIDE &&
                          (*buffers)[new_index].size() > largest_size) {
                        largest_index = i;
                        largest_size = (*buffers)[new_index].size();
                      }
                    }

                  } // for (int i = TRAVELDIRECTION_NUMBER - 1; i >= 0; --i)

                  this_grid.set_largest_buffer(largest_index, largest_size);

                  // add photons that were absorbed (if reemission was disabled)
                  // or that left the system to the global count
                  num_photon_done.pre_add(num_photon_done_now);

                  // delete the original buffer, as we are done with it
                  buffers->free_buffer(current_buffer_index);

                  cmac_assert_message(num_active_buffers.value() > 0,
                                      "Number of active buffers < 0!");
                  num_active_buffers.pre_decrement();
                  // log the end time of the task
                  task.stop();
                }
                task.unlock_dependency();

                cpucycle_tick(task_stop);
                active_time[thread_id] += task_stop - task_start;

                if (task_plot_i >= task_plot_N) {
                  tasks->free_element(current_index);
                }

                for (uint_fast32_t itask = 0; itask < num_tasks_to_add;
                     ++itask) {
                  if (queues_to_add[itask] < 0) {
                    // general queue
                    shared_queue->add_task(tasks_to_add[itask]);
                  } else {
                    queues[queues_to_add[itask]]->add_task(tasks_to_add[itask]);
                  }
                }

                current_index = queues[thread_id]->get_task(*tasks);
                if (current_index == NO_TASK) {
                  current_index = steal_task(thread_id, num_thread, queues,
                                             *tasks, *grid_creator);
                  if (current_index == NO_TASK) {
                    current_index = shared_queue->get_task(*tasks);
                  }
                }
              }

              if (num_empty.value() == num_empty_target &&
                  num_active_buffers.value() == 0 &&
                  num_photon_done.value() == numphoton) {
                global_run_flag = false;
              } else {
                current_index = queues[thread_id]->get_task(*tasks);
                if (current_index == NO_TASK) {
                  current_index = steal_task(thread_id, num_thread, queues,
                                             *tasks, *grid_creator);
                  if (current_index == NO_TASK) {
                    current_index = shared_queue->get_task(*tasks);
                  }
                }
              }
            } // while(global_run_flag)
          }   // parallel region
          stop_parallel_timing_block();

          buffers->reset();

          // update copies
          grid_creator->update_original_counters();

          cmac_assert_message(buffers->is_empty(),
                              "Number of active buffers: %zu",
                              buffers->get_number_of_active_buffers());

          {
            AtomicValue< size_t > igrid(0);
            start_parallel_timing_block();
#pragma omp parallel default(shared)
            while (igrid.value() <
                   grid_creator->number_of_original_subgrids()) {
              const size_t this_igrid = igrid.post_increment();
              if (this_igrid < grid_creator->number_of_original_subgrids()) {
                auto gridit = grid_creator->get_subgrid(this_igrid);

                const size_t itask = tasks->get_free_element();
                Task &task = (*tasks)[itask];

                uint_fast64_t task_start, task_stop;
                cpucycle_tick(task_start);

                task.set_type(TASKTYPE_TEMPERATURE_STATE);
                task.start(omp_get_thread_num());

                // correct the intensity counters for abundance factors
                for (auto cellit = (*gridit).begin(); cellit != (*gridit).end();
                     ++cellit) {
                  IonizationVariables &vars = cellit.get_ionization_variables();
                  for (int_fast32_t ion = 1; ion < NUMBER_OF_IONNAMES; ++ion) {
                    const double abundance =
                        abundances.get_abundance(get_element(ion));
                    if (abundance > 0.) {
                      vars.set_mean_intensity(
                          ion, vars.get_mean_intensity(ion) / abundance);
                    }
                  }
                }
                temperature_calculator->calculate_temperature(iloop, numphoton,
                                                              *gridit);
                task.stop();
                cpucycle_tick(task_stop);
                active_time[omp_get_thread_num()] += task_stop - task_start;
              }
            }
            stop_parallel_timing_block();
          }

          worktimer.stop();
        }

      } else {

        if (log) {
          log->write_status("No ionising sources!");
        }

        // there are no ionising sources: skip radiation for this step
        // manually set all neutral fractions to 1
        {
          AtomicValue< size_t > igrid(0);
          start_parallel_timing_block();
#pragma omp parallel default(shared)
          while (igrid.value() < grid_creator->number_of_original_subgrids()) {
            const size_t this_igrid = igrid.post_increment();
            if (this_igrid < grid_creator->number_of_original_subgrids()) {
              auto gridit = grid_creator->get_subgrid(this_igrid);
              for (auto cellit = (*gridit).begin(); cellit != (*gridit).end();
                   ++cellit) {
                cellit.get_ionization_variables().set_ionic_fraction(ION_H_n,
                                                                     1.);
              }
            }
          }
          stop_parallel_timing_block();
        }
      }

      cmac_assert_message(buffers->is_empty(), "Number of active buffers: %zu",
                          buffers->get_number_of_active_buffers());

      {
        AtomicValue< size_t > igrid(0);
        start_parallel_timing_block();
#pragma omp parallel default(shared)
        while (igrid.value() < grid_creator->number_of_original_subgrids()) {
          const size_t this_igrid = igrid.post_increment();
          if (this_igrid < grid_creator->number_of_original_subgrids()) {
            auto gridit = grid_creator->get_subgrid(this_igrid);
            (*gridit).add_ionization_energy(hydro);
          }
        }
        stop_parallel_timing_block();
      }

      if (log) {
        log->write_status("Done with radiation step.");
      }
    }

    if (log) {
      log->write_status("Starting hydro step...");
    }

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

    start_parallel_timing_block();
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

          uint_fast64_t task_start, task_stop;
          cpucycle_tick(task_start);

          execute_task(current_task, *grid_creator, *tasks, actual_timestep,
                       hydro, hydro_boundary_manager);
          (*tasks)[current_task].stop();

          cpucycle_tick(task_stop);
          active_time[thread_id] += task_stop - task_start;

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
    stop_parallel_timing_block();

    cpucycle_tick(iteration_end);

    if (log) {
      log->write_status("Done with hydro step.");
    }

    if (log) {
      const uint_fast64_t total_interval = iteration_end - iteration_start;
      log->write_status("Thread activity level:");
      for (int_fast32_t ithread = 0; ithread < num_thread; ++ithread) {
        const double percentage =
            (100. * active_time[ithread]) / total_interval;
        log->write_status("Thread ", ithread, ": ", percentage, "%.");
      }
    }

    // write snapshot
    // we don't write if this is the last snapshot, because then it is written
    // outside the integration loop
    if (write_output && hydro_lastsnap * hydro_snaptime <= current_time &&
        has_next_step) {
      writer->write(*grid_creator, hydro_lastsnap, *params, current_time);
      ++hydro_lastsnap;
    }

    if (write_output && task_plot_i < task_plot_N) {
      if (log) {
        log->write_status("Writing task plot file...");
      }
      output_tasks(task_plot_i, *tasks, iteration_start, iteration_end);
      if (log) {
        log->write_status("Done writing task plot file.");
      }
      ++task_plot_i;
    }

    // remove radiation tasks
    tasks->clear_after(radiation_task_offset);

    requested_timestep = DBL_MAX;
    for (auto gridit = grid_creator->begin();
         gridit != grid_creator->original_end(); ++gridit) {
      for (auto cellit = (*gridit).hydro_begin();
           cellit != (*gridit).hydro_end(); ++cellit) {
        requested_timestep =
            std::min(requested_timestep,
                     hydro.get_timestep(cellit.get_hydro_variables(),
                                        cellit.get_ionization_variables(),
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

  serial_timer.stop();
  total_timer.stop();
  programtimer.stop();
  if (log) {
    log->write_status("Total serial time: ",
                      Utilities::human_readable_time(serial_timer.value()),
                      ".");
    log->write_status("Total parallel time: ",
                      Utilities::human_readable_time(parallel_timer.value()),
                      ".");
    log->write_status("Total overall time: ",
                      Utilities::human_readable_time(total_timer.value()), ".");
    log->write_status("Total program time: ",
                      Utilities::human_readable_time(programtimer.value()),
                      ".");

    const size_t memory_usage = OperatingSystem::get_peak_memory_usage();
    log->write_status("Peak memory usage: ",
                      Utilities::human_readable_bytes(memory_usage), ".");
    log->write_status("Total photon shooting time: ",
                      Utilities::human_readable_time(worktimer.value()), ".");
  }

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
  if (reemission_handler != nullptr) {
    delete reemission_handler;
  }

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
