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
#include "AlveliusTurbulenceForcing.hpp"
#include "ChargeTransferRates.hpp"
#include "CommandLineParser.hpp"
#include "ContinuousPhotonSourceFactory.hpp"
#include "CrossSectionsFactory.hpp"
#include "DeRijckeRadiativeCooling.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DiffuseReemissionHandlerFactory.hpp"
#include "DistributedPhotonSource.hpp"
#include "ExternalPotentialFactory.hpp"
#include "HydroBoundaryManager.hpp"
#include "HydroDensitySubGrid.hpp"
#include "HydroMaskFactory.hpp"
#include "LineCoolingData.hpp"
#include "LiveOutputManager.hpp"
#include "MemoryLogger.hpp"
#include "MemorySpace.hpp"
#include "OpenMP.hpp"
#include "ParameterFile.hpp"
#include "PhotonReemitTaskContext.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "PhotonTraversalTaskContext.hpp"
#include "PhotonTraversalThreadContext.hpp"
#include "PrematureLaunchTaskContext.hpp"
#include "RecombinationRatesFactory.hpp"
#include "RestartManager.hpp"
#include "Scheduler.hpp"
#include "SimulationBox.hpp"
#include "SourceDiscretePhotonTaskContext.hpp"
#include "TaskQueue.hpp"
#include "TemperatureCalculator.hpp"
#include "TimeLine.hpp"
#include "TimeLogger.hpp"

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
  std::vector< uint_fast32_t > sorti = Utilities::argsort(queue_sizes);

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
    subgrid.inner_flux_sweep(hydro, timestep);
    break;
  case TASKTYPE_FLUXSWEEP_EXTERNAL_NEIGHBOUR:
    subgrid.outer_flux_sweep(task.get_interaction_direction(), hydro,
                             *grid_creator.get_subgrid(task.get_buffer()),
                             timestep);
    break;
  case TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY:
    subgrid.outer_ghost_flux_sweep(task.get_interaction_direction(), hydro,
                                   boundary_manager.get_boundary_condition(
                                       task.get_interaction_direction()),
                                   timestep);
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
 *  - task-plot-rhd (no abbreviation, optional, integer argument): output task
 *    information for the first N steps of the algorithm (default: 0).
 *  - number-of-steps (no abbreviation, optional, integer argument): number of
 *    time steps to execute before halting the code (negative values mean no
 *    limit on the number of time steps, default: -1).
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
      "task-plot-rhd", 0,
      "Output task information for the first N steps of the algorithm",
      COMMANDLINEOPTION_INTARGUMENT, "0");
  parser.add_option("number-of-steps", 0,
                    "Number of time steps to execute before halting the code.",
                    COMMANDLINEOPTION_INTARGUMENT, "-1");
}

/**
 * @brief Apply the cooling to the given cell over the given time step.
 *
 * We need to solve the implicit equation
 * @f[
 *   E_t(t + \Delta{}t) = E_t(t) - \Delta{}t \Lambda{}(E_t(t + \Delta{}t)),
 * @f]
 * where @f$E_t(t)@f$ is the thermal energy at time @f$t@f$, @f$\Delta{}t@f$
 * is the time step, and @f$\Lambda{}(E_t)@f$ is the cooling function for a
 * given thermal energy.
 *
 * This is essentially a root finding problem that we solve using a bisection
 * method.
 *
 * @param ionization_variables Copy of the IonizationVariables to update
 * (updated with final result).
 * @param hydro_variables Copy of the HydroVariables to update (updated with
 * final result).
 * @param inverse_volume Inverse cell volume (in m^-3).
 * @param nH2V Hydrogen number density squared times cell volume (in m^-3).
 * @param total_dt Total time step (in s).
 * @param radiative_cooling Radiative cooling tables to use.
 * @param hydro Hydro instance to use.
 */
inline static void do_cooling(IonizationVariables &ionization_variables,
                              HydroVariables &hydro_variables,
                              const double inverse_volume, const double nH2V,
                              const double total_dt,
                              DeRijckeRadiativeCooling &radiative_cooling,
                              Hydro &hydro) {

  const double energy_start =
      hydro_variables.get_conserved_total_energy() -
      0.5 * CoordinateVector<>::dot_product(
                hydro_variables.get_primitives_velocity(),
                hydro_variables.get_conserved_momentum());

  if (energy_start == 0.) {
    // don't cool gas that has no thermal energy
    return;
  }

  const double temperature_start = ionization_variables.get_temperature();

  double energy_low = energy_start;
  double energy_high = energy_start;

  double temperature = temperature_start;
  double cooling = radiative_cooling.get_cooling_rate(temperature) * nH2V;

  cmac_assert(energy_high - energy_start + total_dt * cooling > 0.);
  cmac_assert(energy_low - energy_start + total_dt * cooling > 0.);

  uint_fast32_t loopcount = 0;
  while (loopcount < 1e6 &&
         energy_low - energy_start + total_dt * cooling > 0.) {

    ++loopcount;

    if (energy_low == 0.) {
      hydro.update_energy_variables(ionization_variables, hydro_variables,
                                    inverse_volume, -energy_start);
      return;
    }

    energy_high = energy_low;
    energy_low *= 0.5;

    cmac_assert(energy_high - energy_start + total_dt * cooling > 0.);

    temperature =
        temperature_start + hydro.get_temperature_difference(
                                ionization_variables, hydro_variables,
                                inverse_volume, energy_start - energy_low);
    cooling = radiative_cooling.get_cooling_rate(temperature) * nH2V;
  }
  cmac_assert_message(
      loopcount < 1e6,
      "energy_low: %g, energy_start: %g, total_dt: %g, cooling: %g", energy_low,
      energy_start, total_dt, cooling);

  double energy_next = 0.5 * (energy_low + energy_high);
  temperature =
      temperature_start + hydro.get_temperature_difference(
                              ionization_variables, hydro_variables,
                              inverse_volume, energy_start - energy_next);
  cooling = radiative_cooling.get_cooling_rate(temperature) * nH2V;
  loopcount = 0;
  while (loopcount < 1e6 &&
         std::abs(energy_low - energy_high) > 1.e-5 * energy_next) {

    ++loopcount;

    if (energy_next - energy_start + total_dt * cooling > 0.) {
      energy_high = energy_next;
    } else {
      energy_low = energy_next;
    }
    energy_next = 0.5 * (energy_low + energy_high);
    temperature =
        temperature_start + hydro.get_temperature_difference(
                                ionization_variables, hydro_variables,
                                inverse_volume, energy_start - energy_next);
    cooling = radiative_cooling.get_cooling_rate(temperature) * nH2V;
  }
  cmac_assert(loopcount < 1e6);

  double dE = energy_next - energy_start;
  dE = std::max(dE, -0.5 * energy_start);
  hydro.update_energy_variables(ionization_variables, hydro_variables,
                                inverse_volume, dE);
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
 *  - external gravity: Enable external gravity? (default: no)
 *  - use mask: Use a mask to disable hydrodynamics and radiation in part of
 *    the box? (default: no)
 *  - turbulent forcing: Enable turbulent forcing? (default: no)
 *  - first snapshot: Index of the first snapshot to write out (default: 0)
 *  - do radiation: Enable radiation? (default: yes)
 *  - do radiative cooling: Enable radiative cooling? (default: no)
 *  - do stellar feedback: Enable stellar feedback? (default: no)
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

  OperatingSystem::install_signal_handlers(false);

  RestartReader *restart_reader = nullptr;
  if (parser.was_found("restart")) {
    restart_reader = RestartManager::get_restart_reader(
        parser.get_value< std::string >("restart"), log);
    total_timer = Timer(*restart_reader);
    serial_timer = Timer(*restart_reader);
    parallel_timer = Timer(*restart_reader);
    worktimer = Timer(*restart_reader);
    // there is no point in restarting the cycle counts, as they won't make
    // any sense
  }

  uint_fast64_t program_start, program_end;
  cpucycle_tick(program_start);

  total_timer.start();
  serial_timer.start();

  MemoryLogger memory_logger;
  memory_logger.add_entry("start");
  TimeLogger time_logger;

  time_logger.start("initialization");

  time_logger.start("parameter reading");

  const int_fast32_t num_thread = parser.get_value< int_fast32_t >("threads");
  set_number_of_threads(num_thread);

  // set the unit for time stats terminal output
  std::string output_time_unit =
      parser.get_value< std::string >("output-time-unit");

  const int_fast32_t task_plot_N =
      parser.get_value< int_fast32_t >("task-plot-rhd");

  const int_fast32_t number_of_steps =
      parser.get_value< int_fast32_t >("number-of-steps");

  // second: initialize the parameters that are read in from static files
  // these files should be configured by CMake and put in a location that is
  // stored in a CMake configured header
  LineCoolingData line_cooling_data;

  // third: read in the parameters of the run from a parameter file. This file
  // should be read by a ParameterFileParser object that acts as a dictionary
  ParameterFile *params = nullptr;
  if (restart_reader == nullptr) {
    params = new ParameterFile(parser.get_value< std::string >("params"));
  } else {
    params = new ParameterFile(*restart_reader);
  }

  // fourth: construct the density grid. This should be stored in a separate
  // DensityGrid object with geometrical and physical properties
  time_logger.start("density function creation");
  DensityFunction *density_function =
      DensityFunctionFactory::generate(*params, log);
  time_logger.end("density function creation");
  CrossSections *cross_sections = CrossSectionsFactory::generate(*params, log);
  RecombinationRates *recombination_rates =
      RecombinationRatesFactory::generate(*params, log);
  DiffuseReemissionHandler *reemission_handler = nullptr;
  if (params->get_value< bool >(
          "TaskBasedRadiationHydrodynamicsSimulation:diffuse field", false)) {
    reemission_handler = DiffuseReemissionHandlerFactory::generate(
        *cross_sections, *params, log);
  }

  // initialize the simulation box
  const SimulationBox simulation_box(*params);

  ExternalPotential *external_potential = nullptr;
  if (params->get_value< bool >(
          "TaskBasedRadiationHydrodynamicsSimulation:external gravity",
          false)) {
    external_potential = ExternalPotentialFactory::generate(*params, log);
  }
  HydroMask *hydro_mask = nullptr;
  if (params->get_value< bool >(
          "TaskBasedRadiationHydrodynamicsSimulation:use mask", false)) {
    if (restart_reader == nullptr) {
      hydro_mask = HydroMaskFactory::generate(*params, log);
    } else {
      hydro_mask = HydroMaskFactory::restart(*restart_reader, log);
    }
    if (hydro_mask == nullptr) {
      cmac_error("Hydro mask requested, but no hydro mask was found!");
    }
  }

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
  const uint_fast32_t hydro_firstsnap = params->get_value< uint_fast32_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:first snapshot", 0);

  const double hydro_radtime = params->get_physical_value< QUANTITY_TIME >(
      "TaskBasedRadiationHydrodynamicsSimulation:radiation time", "-1. s");
  uint_fast32_t hydro_lastrad = 0;
  const bool do_radiation = params->get_value< bool >(
      "TaskBasedRadiationHydrodynamicsSimulation:do radiation", true);

  if (restart_reader != nullptr) {
    hydro_lastsnap = restart_reader->read< uint_fast32_t >();
    hydro_lastrad = restart_reader->read< uint_fast32_t >();
  }

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
  int_fast32_t random_seed = params->get_value< int_fast32_t >(
      "TaskBasedRadiationHydrodynamicsSimulation:random seed", 42);
  if (restart_reader != nullptr) {
    random_seed = restart_reader->read< int_fast32_t >();
  }
  time_logger.start("density grid creation");
  DensitySubGridCreator< HydroDensitySubGrid > *grid_creator = nullptr;
  if (restart_reader == nullptr) {
    grid_creator = new DensitySubGridCreator< HydroDensitySubGrid >(
        simulation_box.get_box(), *params);
  } else {
    grid_creator =
        new DensitySubGridCreator< HydroDensitySubGrid >(*restart_reader);
  }
  time_logger.end("density grid creation");

  AlveliusTurbulenceForcing *turbulence_forcing = nullptr;
  if (params->get_value< bool >(
          "TaskBasedRadiationHydrodynamicsSimulation:turbulent forcing",
          false)) {
    if (restart_reader == nullptr) {
      time_logger.start("turbulence initialization");
      turbulence_forcing =
          new AlveliusTurbulenceForcing(grid_creator->get_subgrid_layout(),
                                        grid_creator->get_subgrid_cell_layout(),
                                        simulation_box.get_box(), *params, log);
      time_logger.end("turbulence initialization");
    } else {
      turbulence_forcing = new AlveliusTurbulenceForcing(*restart_reader);
    }
  }

  DeRijckeRadiativeCooling *radiative_cooling = nullptr;
  if (params->get_value< bool >(
          "TaskBasedRadiationHydrodynamicsSimulation:do radiative cooling",
          false)) {
    radiative_cooling = new DeRijckeRadiativeCooling();
  }

  const bool do_stellar_feedback = params->get_value< bool >(
      "TaskBasedRadiationHydrodynamicsSimulation:do stellar feedback", false);

  // fifth: construct the stellar sources. These should be stored in a
  // separate StellarSources object with geometrical and physical properties.
  PhotonSourceDistribution *sourcedistribution = nullptr;
  if (restart_reader == nullptr) {
    sourcedistribution =
        PhotonSourceDistributionFactory::generate(*params, log);
  } else {
    sourcedistribution =
        PhotonSourceDistributionFactory::restart(*restart_reader, log);
  }
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
      sourcedistribution->get_total_luminosity(), abundances, line_cooling_data,
      *recombination_rates, charge_transfer_rates, *params, log);

  RestartManager restart_manager(*params);
  RandomGenerator restart_generator(random_seed);

  LiveOutputManager live_output_manager(grid_creator->get_subgrid_layout(),
                                        grid_creator->get_subgrid_cell_layout(),
                                        *params);
  if (restart_reader != nullptr) {
    live_output_manager.read_restart_info(*restart_reader);
  }

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output && restart_reader == nullptr) {
    const std::string usedvaluename =
        parser.get_value< std::string >("params") + ".used-values";
    std::ofstream pfile(usedvaluename);
    params->print_contents(pfile);
    pfile.close();
    if (log) {
      log->write_status("Wrote used parameters to ", usedvaluename, ".");
    }
  }

  time_logger.end("parameter reading");

  if (parser.get_value< bool >("dry-run")) {
    if (log) {
      log->write_warning("Dry run requested. Program will now halt.");
    }
    return 0;
  }

  memory_logger.add_entry("preinit");

  if (restart_reader == nullptr) {
    time_logger.start("grid initialization");
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
    memory_logger.add_entry("grid");
    start_parallel_timing_block();
    grid_creator->initialize(*density_function);
    stop_parallel_timing_block();

#ifdef VARIABLE_ABUNDANCES
    if (!density_function->has_abundances()) {
      for (auto gridit = grid_creator->begin();
           gridit != grid_creator->original_end(); ++gridit) {
        for (auto cellit = (*gridit).begin(); cellit != (*gridit).end();
             ++cellit) {
          cellit.get_ionization_variables().get_abundances().set_abundances(
              abundances);
        }
      }
    }
#endif

    memory_logger.finalize_entry();
    if (log) {
      log->write_status("Done.");
    }

    memory_logger.add_entry("postinit");

    density_function->free();
    time_logger.end("grid initialization");
  }

  memory_logger.add_entry("pretasks");

  time_logger.start("task initialization");

  if (log) {
    log->write_status("Initializing task-based structures...");
  }
  if (log) {
    log->write_status("Allocating memory space...");
  }
  memory_logger.add_entry("memory space");
  MemorySpace *buffers = new MemorySpace(number_of_buffers);
  memory_logger.finalize_entry();
  if (log) {
    log->write_status("Allocating task memory...");
  }
  memory_logger.add_entry("tasks");
  ThreadSafeVector< Task > *tasks =
      new ThreadSafeVector< Task >(number_of_tasks, "Tasks");
  memory_logger.finalize_entry();
  if (log) {
    log->write_status("Allocating shared queue...");
  }
  memory_logger.add_entry("shared queue");
  TaskQueue *shared_queue = new TaskQueue(shared_queue_size, "Shared queue");
  memory_logger.finalize_entry();
  if (log) {
    log->write_status("Allocating per thread queues...");
  }
  memory_logger.add_entry("per thread queues");
  std::vector< TaskQueue * > queues(num_thread);
  for (int_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    std::stringstream queue_name;
    queue_name << "Queue for Thread " << static_cast< int_fast32_t >(ithread);
    queues[ithread] = new TaskQueue(queue_size_per_thread, queue_name.str());
  }
  memory_logger.finalize_entry();
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

  memory_logger.add_entry("posttasks");

  time_logger.end("task initialization");

  time_logger.start("initial time step");
  double requested_timestep = DBL_MAX;
  if (restart_reader == nullptr) {
    for (auto cellit = grid_creator->begin();
         cellit != grid_creator->original_end(); ++cellit) {
      requested_timestep =
          std::min(requested_timestep,
                   (*cellit).initialize_hydrodynamic_variables(hydro, true));
    }
  }
  time_logger.end("initial time step");

  time_logger.start("hydro task creation");
  for (auto cellit = grid_creator->begin();
       cellit != grid_creator->original_end(); ++cellit) {
    make_hydro_tasks(*tasks, cellit.get_index(), *grid_creator);
  }
  for (auto cellit = grid_creator->begin();
       cellit != grid_creator->original_end(); ++cellit) {
    set_dependencies(cellit.get_index(), *grid_creator, *tasks);
  }
  const size_t radiation_task_offset = tasks->size();
  time_logger.end("hydro task creation");

  // initialize the mask (if applicable). MUST BE DONE IN SERIAL!
  if (hydro_mask != nullptr && restart_reader == nullptr) {
    for (auto gridit = grid_creator->begin();
         gridit != grid_creator->original_end(); ++gridit) {
      hydro_mask->initialize_mask(gridit.get_index(), *gridit);
    }
  }

  // apply the mask (if applicable)
  if (hydro_mask != nullptr && restart_reader == nullptr) {
    AtomicValue< size_t > igrid(0);
    start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    while (igrid.value() < grid_creator->number_of_original_subgrids()) {
      const size_t this_igrid = igrid.post_increment();
      if (this_igrid < grid_creator->number_of_original_subgrids()) {
        hydro_mask->apply_mask(this_igrid,
                               *grid_creator->get_subgrid(this_igrid), 0., 0.);
      }
    }
    stop_parallel_timing_block();
  }

  // update the gravitational accelerations if applicable (just to make sure
  // they are present in the first snapshot)
  if (external_potential != nullptr && restart_reader == nullptr) {
    AtomicValue< size_t > igrid(0);
    start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    while (igrid.value() < grid_creator->number_of_original_subgrids()) {
      const size_t this_igrid = igrid.post_increment();
      if (this_igrid < grid_creator->number_of_original_subgrids()) {
        HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
        for (auto it = subgrid.hydro_begin(); it != subgrid.hydro_end(); ++it) {
          const CoordinateVector<> a =
              external_potential->get_acceleration(it.get_cell_midpoint());
          it.get_hydro_variables().set_gravitational_acceleration(a);
        }
      }
    }
    stop_parallel_timing_block();
  }

  // do the initial stellar feedback
  if (restart_reader == nullptr && do_stellar_feedback &&
      sourcedistribution->do_stellar_feedback(0.)) {
    AtomicValue< size_t > igrid(0);
    start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    while (igrid.value() < grid_creator->number_of_original_subgrids()) {
      const size_t this_igrid = igrid.post_increment();
      if (this_igrid < grid_creator->number_of_original_subgrids()) {
        HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
        sourcedistribution->add_stellar_feedback(subgrid);
        auto gridit = grid_creator->get_subgrid(this_igrid);
        for (auto cellit = (*gridit).hydro_begin();
             cellit != (*gridit).hydro_end(); ++cellit) {
          const double dE = cellit.get_hydro_variables().get_energy_term();
          hydro.update_energy_variables(cellit.get_ionization_variables(),
                                        cellit.get_hydro_variables(),
                                        1. / cellit.get_volume(), dE);
          cellit.get_hydro_variables().set_energy_term(0.);
        }
      }
    }
    stop_parallel_timing_block();
    sourcedistribution->done_stellar_feedback();
  }

  if (write_output && restart_reader == nullptr && hydro_firstsnap == 0) {
    time_logger.start("snapshot");
    writer->write(*grid_creator, 0, *params, 0.);
    time_logger.end("snapshot");
  }

  double maximum_timestep = hydro_maximum_timestep;
  if (hydro_radtime > 0.) {
    // make sure the system is evolved hydrodynamically in between successive
    // ionization steps
    maximum_timestep = std::min(maximum_timestep, hydro_radtime);
  }

  if (restart_reader == nullptr) {
    time_logger.start("first time step");
    requested_timestep = DBL_MAX;
    {
      // first figure out the time step for each subgrid, then do the global
      // time step
      std::vector< double > requested_timestep_list(
          grid_creator->number_of_original_subgrids(), DBL_MAX);
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          for (auto cellit = subgrid.hydro_begin();
               cellit != subgrid.hydro_end(); ++cellit) {
            requested_timestep_list[this_igrid] =
                std::min(requested_timestep_list[this_igrid],
                         hydro.get_timestep(cellit.get_hydro_variables(),
                                            cellit.get_ionization_variables(),
                                            cellit.get_volume()));
          }
        }
      }
      stop_parallel_timing_block();
      for (uint_fast32_t i = 0; i < requested_timestep_list.size(); ++i) {
        requested_timestep =
            std::min(requested_timestep, requested_timestep_list[i]);
      }
    }
    time_logger.end("first time step");
  }

  /// RADIATION
  if (restart_reader == nullptr) {
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
    memory_logger.add_entry("subgrid copies");
    grid_creator->create_copies(levels);
    memory_logger.finalize_entry();
  }

  {
    if (log) {
      log->write_status("Outputting memory allocation stats to memory.txt.");
    }
    std::ofstream mfile("memory.txt");
    memory_logger.print(mfile, false);
  }

  /// SIMULATION
  TimeLine *timeline = nullptr;
  int_fast32_t num_step = 0;
  double actual_timestep, current_time;
  requested_timestep *= CFL;
  bool has_next_step;
  if (restart_reader == nullptr) {
    timeline = new TimeLine(0., hydro_total_time, hydro_minimum_timestep,
                            maximum_timestep, log);
    has_next_step =
        timeline->advance(requested_timestep, actual_timestep, current_time);
  } else {
    timeline = new TimeLine(*restart_reader);
    num_step = restart_reader->read< uint_fast32_t >();
    requested_timestep = restart_reader->read< double >();
    has_next_step = restart_reader->read< bool >();
    actual_timestep = restart_reader->read< double >();
    current_time = restart_reader->read< double >();
    delete restart_reader;
    restart_reader = nullptr;
  }

  time_logger.end("initialization");

  time_logger.output("time_log.txt");

  bool stop_simulation = false;
  while (has_next_step && !stop_simulation) {

    // check if we want to prematurely stop the algorithm
    if (number_of_steps > 0 && num_step == number_of_steps) {
      break;
    }

    uint_fast64_t iteration_start, iteration_end;
    cpucycle_tick(iteration_start);
    std::vector< uint_fast64_t > active_time(num_thread, 0);

    ++num_step;
    std::stringstream num_step_line;
    num_step_line << "step " << num_step;
    time_logger.start(num_step_line.str());

    if (number_of_steps > 0) {
      memory_logger.add_entry(num_step_line.str());
    }
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
    if (do_radiation &&
        (hydro_radtime < 0. ||
         (current_time - actual_timestep) >= hydro_lastrad * hydro_radtime)) {

      time_logger.start("radiation");

      if (log) {
        log->write_status("Starting radiation step...");
      }

      ++hydro_lastrad;
      // update the PhotonSource
      if (sourcedistribution->update(current_time)) {
        time_logger.start("source update");

        temperature_calculator->update_luminosity(
            sourcedistribution->get_total_luminosity());

        if (log) {
          log->write_status("Updating subgrid copy hierarchy after source "
                            "distribution change...");
        }

        // update subgrid copies
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
                const uint_fast8_t numngbs =
                    grid_creator->get_neighbours(i, ngbs);
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
        grid_creator->update_copies(levels);

        time_logger.end("source update");
      }

      if (sourcedistribution->get_total_luminosity() > 0.) {
        time_logger.start("radiation transfer");

        {
          AtomicValue< size_t > igrid(0);
          start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
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
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
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

          // reset mean intensity counters
          {
            AtomicValue< size_t > igrid(0);
            start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
            while (igrid.value() < grid_creator->number_of_actual_subgrids()) {
              const size_t this_igrid = igrid.post_increment();
              if (this_igrid < grid_creator->number_of_actual_subgrids()) {
                auto gridit = grid_creator->get_subgrid(this_igrid);
                (*gridit).reset_intensities();
              }
            }
            stop_parallel_timing_block();
          }

          // update copies
          grid_creator->update_copy_properties();

          // reset the photon source information
          photon_source.reset();

          // reset the diffuse field variables
          if (reemission_handler != nullptr) {
            AtomicValue< size_t > igrid(0);
            start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
            while (igrid.value() < grid_creator->number_of_actual_subgrids()) {
              const size_t this_igrid = igrid.post_increment();
              if (this_igrid < grid_creator->number_of_actual_subgrids()) {
                auto gridit = grid_creator->get_subgrid(this_igrid);
                for (auto cellit = (*gridit).begin(); cellit != (*gridit).end();
                     ++cellit) {
                  IonizationVariables &vars = cellit.get_ionization_variables();
                  reemission_handler->set_reemission_probabilities(vars);
                }
              }
            }
            stop_parallel_timing_block();
          }

          size_t number_of_photons_done = 0;
          while (number_of_photons_done < numphoton) {
            for (size_t isrc = 0; isrc < photon_source.get_number_of_sources();
                 ++isrc) {

              const size_t number_of_photons_this_batch =
                  photon_source.get_photon_batch(isrc, PHOTONBUFFER_SIZE);
              if (number_of_photons_this_batch > 0) {
                const size_t new_task = tasks->get_free_element();
                (*tasks)[new_task].set_type(TASKTYPE_SOURCE_DISCRETE_PHOTON);
                (*tasks)[new_task].set_subgrid(isrc);
                (*tasks)[new_task].set_buffer(number_of_photons_this_batch);
                shared_queue->add_task(new_task);
                number_of_photons_done += number_of_photons_this_batch;
              }
            }
          }
          cmac_assert(number_of_photons_done == numphoton);

          bool global_run_flag = true;
          AtomicValue< uint_fast32_t > num_photon_done(0);

          // create task contexts
          TaskContext *task_contexts[TASKTYPE_NUMBER] = {nullptr};
          task_contexts[TASKTYPE_SOURCE_DISCRETE_PHOTON] =
              new SourceDiscretePhotonTaskContext< HydroDensitySubGrid >(
                  photon_source, *buffers, random_generators, 1., *spectrum,
                  abundances, *cross_sections, *grid_creator, *tasks);
          if (reemission_handler) {
            task_contexts[TASKTYPE_PHOTON_REEMIT] =
                new PhotonReemitTaskContext< HydroDensitySubGrid >(
                    *buffers, random_generators, *reemission_handler,
                    abundances, *cross_sections, *grid_creator, *tasks,
                    num_photon_done);
          }
          task_contexts[TASKTYPE_PHOTON_TRAVERSAL] =
              new PhotonTraversalTaskContext< HydroDensitySubGrid >(
                  *buffers, *grid_creator, *tasks, num_photon_done, nullptr,
                  reemission_handler != nullptr);

          PrematureLaunchTaskContext< HydroDensitySubGrid > premature_launch(
              *buffers, *grid_creator, *tasks, queues, *shared_queue);

          Scheduler scheduler(*tasks, queues, *shared_queue);

          start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
          {
            // thread initialisation
            const int_fast8_t thread_id = get_thread_index();
            ThreadContext *thread_contexts[TASKTYPE_NUMBER] = {nullptr};
            for (int_fast32_t itask = 0; itask < TASKTYPE_NUMBER; ++itask) {
              if (task_contexts[itask]) {
                thread_contexts[itask] =
                    task_contexts[itask]->get_thread_context();
              }
            }

            // actual run flag
            uint_fast32_t current_index = shared_queue->get_task(*tasks);
            while (global_run_flag) {

              if (current_index == NO_TASK) {
                premature_launch.execute();
                current_index = scheduler.get_task(thread_id);
              }

              while (current_index != NO_TASK) {

                // execute task
                uint_fast32_t num_tasks_to_add = 0;
                uint_fast32_t tasks_to_add[TRAVELDIRECTION_NUMBER];
                int_fast32_t queues_to_add[TRAVELDIRECTION_NUMBER];

                Task &task = (*tasks)[current_index];
                uint_fast64_t task_start, task_stop;
                cpucycle_tick(task_start);

                task.start(thread_id);

                num_tasks_to_add = task_contexts[task.get_type()]->execute(
                    thread_id, thread_contexts[task.get_type()], tasks_to_add,
                    queues_to_add, task);

                // log the end time of the task
                task.stop();

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

                current_index = scheduler.get_task(thread_id);
              }

              if (buffers->is_empty() && num_photon_done.value() == numphoton) {
                global_run_flag = false;
              } else {
                current_index = scheduler.get_task(thread_id);
              }
            } // while(global_run_flag)

            for (int_fast32_t itask = 0; itask < TASKTYPE_NUMBER; ++itask) {
              delete thread_contexts[itask];
            }
          } // parallel region
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
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
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
                task.start(get_thread_index());

#ifndef VARIABLE_ABUNDANCES
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
#endif
                temperature_calculator->calculate_temperature(iloop, numphoton,
                                                              *gridit);
                task.stop();
                cpucycle_tick(task_stop);
                active_time[get_thread_index()] += task_stop - task_start;
              }
            }
            stop_parallel_timing_block();
          }

          for (int_fast32_t itask = 0; itask < TASKTYPE_NUMBER; ++itask) {
            delete task_contexts[itask];
          }

          worktimer.stop();
        }

        time_logger.end("radiation transfer");

      } else {

        if (log) {
          log->write_status("No ionising sources!");
        }

        // there are no ionising sources: skip radiation for this step
        // manually set all neutral fractions to 1
        {
          AtomicValue< size_t > igrid(0);
          start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
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
        time_logger.start("ionizing energy update");
        AtomicValue< size_t > igrid(0);
        start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
        while (igrid.value() < grid_creator->number_of_original_subgrids()) {
          const size_t this_igrid = igrid.post_increment();
          if (this_igrid < grid_creator->number_of_original_subgrids()) {
            auto gridit = grid_creator->get_subgrid(this_igrid);
            (*gridit).add_ionization_energy(hydro, actual_timestep);
          }
        }
        stop_parallel_timing_block();
        time_logger.end("ionizing energy update");
      }

      if (log) {
        log->write_status("Done with radiation step.");
      }

      time_logger.end("radiation");
    }

    if (radiative_cooling != nullptr) {
      time_logger.start("cooling");
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          auto gridit = grid_creator->get_subgrid(this_igrid);
          for (auto cellit = (*gridit).hydro_begin();
               cellit != (*gridit).hydro_end(); ++cellit) {
            IonizationVariables ionization_variables =
                cellit.get_ionization_variables();
            HydroVariables hydro_variables = cellit.get_hydro_variables();
            const double nH = ionization_variables.get_number_density();
            const double nH2 = nH * nH;
            do_cooling(ionization_variables, hydro_variables,
                       1. / cellit.get_volume(), nH2 * cellit.get_volume(),
                       actual_timestep, *radiative_cooling, hydro);
            cellit.get_ionization_variables().set_temperature(
                ionization_variables.get_temperature());
            cellit.get_hydro_variables().set_primitives_pressure(
                hydro_variables.get_primitives_pressure());
            cellit.get_hydro_variables().set_conserved_total_energy(
                hydro_variables.get_conserved_total_energy());
            if (ionization_variables.get_temperature() <
                radiative_cooling->get_minimum_temperature()) {
              hydro.set_temperature(
                  cellit.get_ionization_variables(),
                  cellit.get_hydro_variables(), cellit.get_volume(),
                  radiative_cooling->get_minimum_temperature());
            }
          }
        }
      }
      stop_parallel_timing_block();
      time_logger.end("cooling");
    }

    time_logger.start("hydro");

    if (log) {
      log->write_status("Starting hydro step...");
    }

    // update the gravitational accelerations if applicable
    if (external_potential != nullptr) {
      time_logger.start("gravity");
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          uint_fast64_t task_start, task_stop;
          cpucycle_tick(task_start);
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          for (auto it = subgrid.hydro_begin(); it != subgrid.hydro_end();
               ++it) {
            const CoordinateVector<> a =
                external_potential->get_acceleration(it.get_cell_midpoint());
            it.get_hydro_variables().set_gravitational_acceleration(a);
          }
          cpucycle_tick(task_stop);
          active_time[get_thread_index()] += task_stop - task_start;
        }
      }
      stop_parallel_timing_block();
      time_logger.end("gravity");
    }

    // apply the turbulent forcing if applicable
    if (turbulence_forcing != nullptr) {
      time_logger.start("turbulence");
      if (log) {
        log->write_status("Applying turbulence forcing...");
      }
      time_logger.start("update");
      turbulence_forcing->update_turbulence(current_time + actual_timestep);
      time_logger.end("update");
      time_logger.start("grid update");
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          uint_fast64_t task_start, task_stop;
          cpucycle_tick(task_start);
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          turbulence_forcing->add_turbulent_forcing(this_igrid, subgrid);
          cpucycle_tick(task_stop);
          active_time[get_thread_index()] += task_stop - task_start;
        }
      }
      stop_parallel_timing_block();
      time_logger.end("grid update");
      if (log) {
        log->write_status("Done applying turbulence forcing.");
      }
      time_logger.end("turbulence");
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
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    {
      const int_fast32_t thread_id = get_thread_index();
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

    // apply the mask (if applicable)
    if (hydro_mask != nullptr) {
      time_logger.start("mask");
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          hydro_mask->apply_mask(this_igrid,
                                 *grid_creator->get_subgrid(this_igrid),
                                 actual_timestep, current_time);
        }
      }
      stop_parallel_timing_block();
      time_logger.end("mask");
    }

    cpucycle_tick(iteration_end);

    if (log) {
      log->write_status("Done with hydro step.");
    }

    time_logger.end("hydro");

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
      if (hydro_firstsnap <= hydro_lastsnap) {
        time_logger.start("snapshot");
        writer->write(*grid_creator, hydro_lastsnap, *params, current_time);
        time_logger.end("snapshot");
      }
      ++hydro_lastsnap;
    }

    // check for live output
    if (live_output_manager.do_output(current_time)) {
      time_logger.start("live output");
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          live_output_manager.compute_output(
              this_igrid, *grid_creator->get_subgrid(this_igrid));
        }
      }
      stop_parallel_timing_block();
      live_output_manager.write_output(simulation_box.get_box());
      time_logger.end("live output");
    }

    if (write_output && task_plot_i < task_plot_N) {
      time_logger.start("task plot");
      if (log) {
        log->write_status("Writing task plot file...");
      }
      output_tasks(task_plot_i, *tasks, iteration_start, iteration_end);
      if (log) {
        log->write_status("Done writing task plot file.");
      }
      time_logger.end("task plot");
      ++task_plot_i;
    }

    // remove radiation tasks
    time_logger.start("task cleanup");
    tasks->clear_after(radiation_task_offset);
    time_logger.end("task cleanup");

    if (do_stellar_feedback &&
        sourcedistribution->do_stellar_feedback(current_time)) {
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          sourcedistribution->add_stellar_feedback(subgrid);
          auto gridit = grid_creator->get_subgrid(this_igrid);
          for (auto cellit = (*gridit).hydro_begin();
               cellit != (*gridit).hydro_end(); ++cellit) {
            const double dE = cellit.get_hydro_variables().get_energy_term();
            hydro.update_energy_variables(cellit.get_ionization_variables(),
                                          cellit.get_hydro_variables(),
                                          1. / cellit.get_volume(), dE);
            cellit.get_hydro_variables().set_energy_term(0.);
          }
        }
      }
      stop_parallel_timing_block();
      sourcedistribution->done_stellar_feedback();
    }

    time_logger.start("time step");
    requested_timestep = DBL_MAX;
    {
      // first figure out the time step for each subgrid, then do the global
      // time step
      std::vector< double > requested_timestep_list(
          grid_creator->number_of_original_subgrids(), DBL_MAX);
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < grid_creator->number_of_original_subgrids()) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          for (auto cellit = subgrid.hydro_begin();
               cellit != subgrid.hydro_end(); ++cellit) {
            requested_timestep_list[this_igrid] =
                std::min(requested_timestep_list[this_igrid],
                         hydro.get_timestep(cellit.get_hydro_variables(),
                                            cellit.get_ionization_variables(),
                                            cellit.get_volume()));
          }
        }
      }
      stop_parallel_timing_block();
      for (uint_fast32_t i = 0; i < requested_timestep_list.size(); ++i) {
        requested_timestep =
            std::min(requested_timestep, requested_timestep_list[i]);
      }
    }
    requested_timestep *= CFL;
    has_next_step =
        timeline->advance(requested_timestep, actual_timestep, current_time);
    time_logger.end("time step");

    random_seed = restart_generator.get_random_integer();
    stop_simulation = restart_manager.stop_simulation() || Signals::interrupt();
    if (restart_manager.write_restart_file() || stop_simulation) {
      time_logger.start("restart file");
      RestartWriter *restart_writer = restart_manager.get_restart_writer(log);

      total_timer.stop();
      total_timer.write_restart_file(*restart_writer);
      total_timer.start();

      serial_timer.stop();
      serial_timer.write_restart_file(*restart_writer);
      serial_timer.start();

      parallel_timer.write_restart_file(*restart_writer);
      worktimer.write_restart_file(*restart_writer);

      params->write_restart_file(*restart_writer);

      if (hydro_mask != nullptr) {
        HydroMaskFactory::write_restart_file(*restart_writer, *hydro_mask);
      }

      restart_writer->write(hydro_lastsnap);
      restart_writer->write(hydro_lastrad);
      restart_writer->write(random_seed);

      grid_creator->write_restart_file(*restart_writer);

      if (turbulence_forcing != nullptr) {
        turbulence_forcing->write_restart_file(*restart_writer);
      }

      PhotonSourceDistributionFactory::write_restart_file(*restart_writer,
                                                          *sourcedistribution);

      live_output_manager.write_restart_info(*restart_writer);

      timeline->write_restart_file(*restart_writer);
      restart_writer->write(num_step);
      restart_writer->write(requested_timestep);
      restart_writer->write(has_next_step);
      restart_writer->write(actual_timestep);
      restart_writer->write(current_time);

      delete restart_writer;
      time_logger.end("restart file");
    }

    time_logger.end(num_step_line.str());
    time_logger.output("time_log.txt", true);
  }

  if (stop_simulation) {
    if (log) {
      log->write_status("Prematurely stopping simulation on request.");
    }
    restart_manager.resubmit();
  } else {
    // the final snapshot is always written
    if (write_output) {
      time_logger.start("snapshot");
      writer->write(*grid_creator, hydro_lastsnap, *params, current_time);
      time_logger.end("snapshot");
    }
  }

  cpucycle_tick(program_end);

  time_logger.output("time_log.txt", true);

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

  memory_logger.add_entry("end");

  if (number_of_steps > 0) {
    {
      std::ofstream pfile("program_time.txt");
      pfile << "# rank\tstart\tstop\ttime\n";
      pfile << "0\t" << program_start << "\t" << program_end << "\t"
            << total_timer.value() << "\n";
    }

    {
      std::ofstream mfile("memory_timeline.txt");
      memory_logger.print(mfile, true);
    }
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
  if (external_potential != nullptr) {
    delete external_potential;
  }
  if (hydro_mask != nullptr) {
    delete hydro_mask;
  }
  if (turbulence_forcing != nullptr) {
    delete turbulence_forcing;
  }
  if (radiative_cooling != nullptr) {
    delete radiative_cooling;
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
