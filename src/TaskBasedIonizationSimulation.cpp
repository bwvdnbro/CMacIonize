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
#include "AbundanceModelFactory.hpp"
#include "ContinuousPhotonSourceFactory.hpp"
#include "CrossSectionsFactory.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DensitySubGrid.hpp"
#include "DensitySubGridCreator.hpp"
#include "DiffuseReemissionHandlerFactory.hpp"
#include "DistributedPhotonSource.hpp"
#include "MemorySpace.hpp"
#include "OpenMP.hpp"
#include "ParameterFile.hpp"
#include "PhotonPacketStatistics.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "RecombinationRatesFactory.hpp"
#include "Signals.hpp"
#include "SimulationBox.hpp"
#include "TaskQueue.hpp"
#include "TemperatureCalculator.hpp"
#include "ThreadStats.hpp"
#include "TrackerManager.hpp"

#include <fstream>
#include <sstream>

/*! @brief Stop the serial time timer and start the parallel time timer. */
#define start_parallel_timing_block()                                          \
  _serial_timer.stop();                                                        \
  _parallel_timer.start();

/*! @brief Stop the parallel time timer and start the serial time timer. */
#define stop_parallel_timing_block()                                           \
  _parallel_timer.stop();                                                      \
  _serial_timer.start();

/*! @brief Uncomment to enable stop condition output. */
//#define OUTPUT_STOP_CONDITION

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
 * @brief Write file with queue size information for an iteration.
 *
 * @param iloop Iteration number (added to file names).
 * @param queues Per thread queues.
 * @param general_queue General queue.
 */
inline void output_queues(const unsigned int iloop,
                          std::vector< TaskQueue * > &queues,
                          TaskQueue &general_queue) {

  // first compose the file name
  std::stringstream filename;
  filename << "queues_";
  filename.fill('0');
  filename.width(2);
  filename << iloop;
  filename << ".txt";

  // now output
  // open the file
  std::ofstream ofile(filename.str(), std::ofstream::trunc);

  ofile << "# rank\tqueue\tsize\n";

  // start with the general queue (-1)
  ofile << "0\t-1\t" << general_queue.get_max_queue_size() << "\n";
  general_queue.reset_max_queue_size();

  // now do the other queues
  for (size_t i = 0; i < queues.size(); ++i) {
    TaskQueue &queue = *queues[i];
    ofile << "0\t" << i << "\t" << queue.get_max_queue_size() << "\n";
    queue.reset_max_queue_size();
  }
}

/**
 * @brief Get the next task for the given thread.
 *
 * We try to get a task in different stages:
 *  - first, we try to get a task from the thread's own queue
 *  - next, we try to steal a task from another thread's queue
 *  - if all else fails, we try to get a task from the shared queue
 * If this doesn't yield a task, we give up and assume there are no more tasks
 * available to this thread.
 *
 * @param thread_id Thread id of this thread.
 * @return Index of a task or NO_TASK if no tasks are available.
 */
uint_fast32_t
TaskBasedIonizationSimulation::get_task(const int_fast8_t thread_id) {

  uint_fast32_t task_index = _queues[thread_id]->get_task(*_tasks);
  if (task_index == NO_TASK) {

    // try to steal a task from another thread's queue

    // sort the queues by size
    std::vector< size_t > queue_sizes(_queues.size(), 0);
    for (size_t i = 0; i < _queues.size(); ++i) {
      queue_sizes[i] = _queues[i]->size();
    }
    std::vector< uint_fast32_t > sorti = Utilities::argsort(queue_sizes);

    // now try to steal from the largest queue first
    uint_fast32_t i = 0;
    while (task_index == NO_TASK && i < queue_sizes.size() &&
           queue_sizes[sorti[queue_sizes.size() - i - 1]] > 0) {
      task_index =
          _queues[sorti[queue_sizes.size() - i - 1]]->try_get_task(*_tasks);
      ++i;
    }
    if (task_index != NO_TASK) {
      // stealing means transferring ownership...
      if ((*_tasks)[task_index].get_type() == TASKTYPE_PHOTON_TRAVERSAL) {
        (*_grid_creator->get_subgrid((*_tasks)[task_index].get_subgrid()))
            .set_owning_thread(thread_id);
      }
    } else {
      // get a task from the shared queue
      task_index = _shared_queue->get_task(*_tasks);
    }
  }

  return task_index;
}

/**
 * @brief Output diagnostic information about the task execution (if activated).
 *
 * @param verbose Is verbose output activated?
 * @param log Log to write logging info to.
 * @param thread_id Thread ID (only thread 0 writes diagnostic info).
 * @param verbose_timer Timer used to regulate diagnostic output (output is
 * written every minute).
 * @param num_empty Atomic counter for the number of inactive subgrid buffers.
 * @param num_empty_target Target number of inactive subgrid buffers.
 * @param num_active_buffers Number of active buffers not related to any
 * subgrid.
 * @param num_photon_done Total number of photon packets that has been
 * processed.
 * @param num_photon_target Target number of photon packets to process.
 * @param verbose_last_num_empty Number of empty photon buffers last time
 * around (to detect if we are stuck).
 * @param verbose_last_num_active_buffers Number of active buffers last time
 * around (to detect if we are stuck).
 * @param verbose_last_num_photon_done Number of processed photon packets last
 * time around (to detect if we are stuck).
 * @param shared_queue Shared queue (to output its statistics).
 * @param queues Thread queues (to output their statistics).
 * @param tasks Tasks (to output the unfinished ones).
 * @param grid_creator Subgrids (to access their photon buffers).
 */
inline void task_status(const bool verbose, Log *log,
                        const int_fast32_t thread_id, Timer &verbose_timer,
                        AtomicValue< uint_fast32_t > &num_empty,
                        const uint_fast32_t num_empty_target,
                        AtomicValue< uint_fast32_t > &num_active_buffers,
                        AtomicValue< uint_fast32_t > &num_photon_done,
                        const uint_fast32_t num_photon_target,
                        uint_fast32_t &verbose_last_num_empty,
                        uint_fast32_t &verbose_last_num_active_buffers,
                        uint_fast32_t &verbose_last_num_photon_done,
                        TaskQueue &shared_queue,
                        std::vector< TaskQueue * > &queues,
                        ThreadSafeVector< Task > &tasks,
                        DensitySubGridCreator< DensitySubGrid > &grid_creator) {

  if (verbose && log != nullptr && thread_id == 0) {
    if (verbose_timer.interval() > 60.) {
      const uint_fast32_t current_num_empty = num_empty.value();
      const uint_fast32_t current_num_active_buffers =
          num_active_buffers.value();
      const uint_fast32_t current_num_photon_done = num_photon_done.value();
      log->write_info("num_empty: ", current_num_empty, " (", num_empty_target,
                      "), num_active_buffers: ", current_num_active_buffers,
                      ", num_photon_done: ", current_num_photon_done, " (",
                      num_photon_target, ")");
      if (current_num_empty == verbose_last_num_empty &&
          current_num_active_buffers == verbose_last_num_active_buffers &&
          current_num_photon_done == verbose_last_num_photon_done) {
        // This is curious. We might be deadlocked. Output additional
        // information.
        log->write_info("Shared queue size: ", shared_queue.size());
        log->write_info("Thread queue sizes:");
        for (uint_fast32_t ithread = 0; ithread < queues.size(); ++ithread) {
          log->write_info("queue[", ithread, "]: ", queues[ithread]->size());
        }
        const size_t current_num_tasks = tasks.get_number_of_active_elements();
        log->write_info("Number of unfinished tasks: ", current_num_tasks);
        Task **current_tasks = new Task *[current_num_tasks];
        const size_t number_of_tasks_retrieved =
            tasks.get_active_elements(current_num_tasks, current_tasks);
        for (size_t itask = 0; itask < number_of_tasks_retrieved; ++itask) {
          log->write_info("task[", itask,
                          "]: ", current_tasks[itask]->get_type());
        }
        delete[] current_tasks;
        log->write_info("Subgrid buffers:");
        for (auto gridit = grid_creator.begin();
             gridit != grid_creator.all_end(); ++gridit) {
          DensitySubGrid &this_subgrid = *gridit;
          log->write_info("subgrid[", gridit.get_index(),
                          "]: ", this_subgrid.get_largest_buffer_size());
        }
      }
      verbose_last_num_empty = current_num_empty;
      verbose_last_num_active_buffers = current_num_active_buffers;
      verbose_last_num_photon_done = current_num_photon_done;
      // reset the timer
      verbose_timer.start();
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
 *  - number of iterations: Number of iterations of the photoionization
 *    algorithm to perform (default: 10)
 *  - number of photons: Number of photon packets to use for each iteration of
 *    the photoionization algorithm (default: 1e6)
 *  - output folder: Folder where all output files will be placed (default: .)
 *  - diffuse field: Should the diffuse field be tracked? (default: false)
 *  - source copy level: Copy level for subgrids that contain a source (default:
 *    4)
 *  - enable trackers: Track photon packets travelling through specific
 *    positions? (default: no)
 *
 * @param num_thread Number of shared memory parallel threads to use.
 * @param parameterfile_name Name of the parameter file to use.
 * @param task_plot Output task plot information?
 * @param verbose Output detailed diagnostic output to the standard output?
 * @param output_initial_snapshot Output a snapshot before the initial
 * iteration?
 * @param log Log to write logging info to.
 */
TaskBasedIonizationSimulation::TaskBasedIonizationSimulation(
    const int_fast32_t num_thread, const std::string parameterfile_name,
    const bool task_plot, const bool verbose,
    const bool output_initial_snapshot, Log *log)
    : _parameter_file(parameterfile_name),
      _number_of_iterations(_parameter_file.get_value< uint_fast32_t >(
          "TaskBasedIonizationSimulation:number of iterations", 10)),
      _number_of_photons(_parameter_file.get_value< uint_fast64_t >(
          "TaskBasedIonizationSimulation:number of photons", 1e6)),
      _source_copy_level(_parameter_file.get_value< uint_fast32_t >(
          "TaskBasedIonizationSimulation:source copy level", 4)),
      _simulation_box(_parameter_file),
      _abundance_model(AbundanceModelFactory::generate(_parameter_file, log)),
      _abundances(_abundance_model->get_abundances()), _log(log),
      _task_plot(task_plot), _verbose(verbose),
      _output_initial_snapshot(output_initial_snapshot) {

  set_number_of_threads(num_thread);

  // install signal handlers
  OperatingSystem::install_signal_handlers(true);

  cpucycle_tick(_program_start);
  _total_timer.start();
  _serial_timer.start();

  _memory_log.add_entry("start");
  _time_log.start("simulation start");

  _time_log.start("memory space");
  const size_t number_of_buffers = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:number of buffers", 50000);
  _memory_log.add_entry("memory space");
  _buffers = new MemorySpace(number_of_buffers);
  _memory_log.finalize_entry();
  _time_log.end("memory space");

  _time_log.start("thread queues");
  const size_t queue_size_per_thread = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:queue size per thread", 10000);
  _memory_log.add_entry("thread queues");
  _queues.resize(num_thread);
  for (int_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    std::stringstream queue_name;
    queue_name << "Queue for Thread " << static_cast< int_fast32_t >(ithread);
    _queues[ithread] = new TaskQueue(queue_size_per_thread, queue_name.str());
  }
  _memory_log.finalize_entry();
  _time_log.end("thread queues");

  _time_log.start("shared queue");
  const size_t shared_queue_size = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:shared queue size", 100000);
  _memory_log.add_entry("shared queue");
  _shared_queue = new TaskQueue(shared_queue_size, "Shared queue");
  _memory_log.finalize_entry();
  _time_log.end("shared queue");

  _time_log.start("tasks");
  const size_t number_of_tasks = _parameter_file.get_value< size_t >(
      "TaskBasedIonizationSimulation:number of tasks", 500000);
  _memory_log.add_entry("tasks");
  _tasks = new ThreadSafeVector< Task >(number_of_tasks, "Tasks");
  _memory_log.finalize_entry();
  _time_log.end("tasks");

  _random_generators.resize(num_thread);
  const int_fast32_t random_seed = _parameter_file.get_value< int_fast32_t >(
      "TaskBasedIonizationSimulation:random seed", 42);
  for (uint_fast8_t ithread = 0; ithread < num_thread; ++ithread) {
    _random_generators[ithread].set_seed(random_seed + ithread);
  }

  _time_log.start("grid creator");
  _grid_creator = new DensitySubGridCreator< DensitySubGrid >(
      _simulation_box.get_box(), _parameter_file);
  _time_log.end("grid creator");

  _time_log.start("density function");
  _density_function = DensityFunctionFactory::generate(_parameter_file, _log);
  _time_log.end("density function");

  // set up output
  std::string output_folder =
      Utilities::get_absolute_path(_parameter_file.get_value< std::string >(
          "TaskBasedIonizationSimulation:output folder", "."));
  _density_grid_writer = DensityGridWriterFactory::generate(
      output_folder, _parameter_file, false, _log);

  _photon_source_distribution =
      PhotonSourceDistributionFactory::generate(_parameter_file, _log);
  _photon_source_spectrum = PhotonSourceSpectrumFactory::generate(
      "PhotonSourceSpectrum", _parameter_file, _log);

  // create the continuous UV sources
  _continuous_photon_source = ContinuousPhotonSourceFactory::generate(
      _simulation_box.get_box(), _parameter_file, _log);
  _continuous_photon_source_spectrum = PhotonSourceSpectrumFactory::generate(
      "ContinuousPhotonSourceSpectrum", _parameter_file, _log);

  _cross_sections = CrossSectionsFactory::generate(_parameter_file, _log);
  _recombination_rates =
      RecombinationRatesFactory::generate(_parameter_file, _log);

  _total_luminosity = 0.;
  if (_photon_source_distribution) {
    _total_luminosity += _photon_source_distribution->get_total_luminosity();
  }
  if (_continuous_photon_source) {
    if (_continuous_photon_source->has_total_luminosity()) {
      _total_luminosity += _continuous_photon_source->get_total_luminosity();
    } else {
      _total_luminosity +=
          _continuous_photon_source_spectrum->get_total_flux() *
          _continuous_photon_source->get_total_surface_area();
    }
  }
  // used to calculate both the ionization state and the temperature
  _temperature_calculator = new TemperatureCalculator(
      _total_luminosity, _abundances, _line_cooling_data, *_recombination_rates,
      _charge_transfer_rates, _parameter_file, _log);

  // the second condition is necessary to deal with old parameter files
  if (_parameter_file.get_value< bool >(
          "TaskBasedIonizationSimulation:diffuse field", false) ||
      _parameter_file.has_value("PhotonSource:diffuse field")) {
    _reemission_handler = DiffuseReemissionHandlerFactory::generate(
        *_cross_sections, _parameter_file, _log);
  } else {
    _reemission_handler = nullptr;
  }

  if (_parameter_file.get_value< bool >(
          "TaskBasedIonizationSimulation:enable trackers", false)) {
    _trackers = new TrackerManager(_parameter_file);
  } else {
    _trackers = nullptr;
  }

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  const std::string usedvaluename = parameterfile_name + ".used-values";
  std::ofstream pfile(usedvaluename);
  _parameter_file.print_contents(pfile);
  pfile.close();
  if (_log) {
    _log->write_status("Wrote used parameters to ", usedvaluename, ".");
  }

  _memory_log.add_entry("parameters done");
  _time_log.end("simulation start");
}

/**
 * @brief Destructor.
 */
TaskBasedIonizationSimulation::~TaskBasedIonizationSimulation() {

  uint_fast64_t program_end;
  cpucycle_tick(program_end);

  _memory_log.add_entry("end");
  _serial_timer.stop();
  _total_timer.stop();

  if (_log) {
    _log->write_status("Total serial time: ",
                       Utilities::human_readable_time(_serial_timer.value()),
                       ".");
    _log->write_status("Total parallel time: ",
                       Utilities::human_readable_time(_parallel_timer.value()),
                       ".");
    _log->write_status("Total overall time: ",
                       Utilities::human_readable_time(_total_timer.value()),
                       ".");
    _log->write_status(
        "Total photon shooting time: ",
        Utilities::human_readable_time(_photon_propagation_timer.value()), ".");
    _log->write_status(
        "Total cell update time: ",
        Utilities::human_readable_time(_cell_update_timer.value()), ".");
  }

  if (_task_plot) {
    std::ofstream pfile("program_time.txt");
    pfile << "# rank\tstart\tstop\ttime\n";
    pfile << "0\t" << _program_start << "\t" << program_end << "\t"
          << _total_timer.value() << "\n";
  }

  {
    std::ofstream mfile("memory.txt");
    _memory_log.print(mfile, false);
  }

  {
    std::ofstream mfile("memory_timeline.txt");
    _memory_log.print(mfile, true);
  }

  _time_log.output("time_log.txt", false);

  delete _buffers;
  for (uint_fast8_t ithread = 0; ithread < _queues.size(); ++ithread) {
    delete _queues[ithread];
  }
  delete _shared_queue;
  delete _tasks;
  delete _grid_creator;
  delete _density_function;
  delete _density_grid_writer;
  delete _photon_source_distribution;
  delete _photon_source_spectrum;
  delete _continuous_photon_source_spectrum;
  delete _continuous_photon_source;
  delete _temperature_calculator;
  delete _cross_sections;
  delete _recombination_rates;
  delete _reemission_handler;
  delete _trackers;
  delete _abundance_model;
}

/**
 * @brief Initialize a photoionization simulation.
 *
 * @param density_function DensityFunction to use to initialize the density
 * field.
 */
void TaskBasedIonizationSimulation::initialize(
    DensityFunction *density_function) {

  _time_log.start("initialisation");

  if (density_function == nullptr) {
    density_function = _density_function;
  }

  _time_log.start("density function");
  density_function->initialize();
  _memory_log.add_entry("density function");
  _time_log.end("density function");

  _time_log.start("grid");
  _memory_log.add_entry("grid");
  start_parallel_timing_block();
  _grid_creator->initialize(*density_function);
  stop_parallel_timing_block();

  if (_log) {
    _log->write_status("Task-based structure sizes:");
    _log->write_status("DensitySubGrid: ",
                       Utilities::human_readable_bytes(
                           (*_grid_creator->begin()).get_memory_size()));
    _log->write_status("Single cell: ", Utilities::human_readable_bytes(
                                            sizeof(IonizationVariables)));
    _log->write_status("PhotonBuffer: ",
                       Utilities::human_readable_bytes(sizeof(PhotonBuffer)));
  }

#ifdef VARIABLE_ABUNDANCES
  for (auto gridit = _grid_creator->begin();
       gridit != _grid_creator->original_end(); ++gridit) {
    for (auto cellit = (*gridit).begin(); cellit != (*gridit).end(); ++cellit) {
      cellit.get_ionization_variables().get_abundances().set_abundances(
          _abundances);
    }
  }
#endif

  _memory_log.finalize_entry();
  _time_log.end("grid");

  _memory_log.add_entry("post initialisation");

  density_function->free();
  _memory_log.add_entry("pre run");

  _time_log.end("initialisation");
}

/**
 * @brief Run a photoionization simulation.
 *
 * @param density_grid_writer DensityGridWriter to use to output the results.
 */
void TaskBasedIonizationSimulation::run(
    DensityGridWriter *density_grid_writer) {

  _time_log.start("main run");

  // write the initial state of the grid to an output file (only do this if
  // we are not in library mode)
  if (_density_grid_writer && _output_initial_snapshot) {
    _time_log.start("snapshot");
    _density_grid_writer->write(*_grid_creator, 0, _parameter_file);
    _time_log.end("snapshot");
  }

  if (density_grid_writer == nullptr) {
    density_grid_writer = _density_grid_writer;
  }

  _time_log.start("subgrid copies");
  std::vector< uint_fast8_t > levels(
      _grid_creator->number_of_original_subgrids(), 0);

  // set the copy level of all subgrids containing a source to the given
  // parameter value (for now)
  if (_photon_source_distribution) {
    const photonsourcenumber_t number_of_sources =
        _photon_source_distribution->get_number_of_sources();
    for (photonsourcenumber_t isource = 0; isource < number_of_sources;
         ++isource) {
      const CoordinateVector<> position =
          _photon_source_distribution->get_position(isource);
      DensitySubGridCreator< DensitySubGrid >::iterator gridit =
          _grid_creator->get_subgrid(position);
      levels[gridit.get_index()] = _source_copy_level;
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
          const uint_fast8_t numngbs = _grid_creator->get_neighbours(i, ngbs);
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
  _memory_log.add_entry("subgrid copies");
  _grid_creator->create_copies(levels);
  _memory_log.finalize_entry();
  _time_log.end("subgrid copies");
  if (_log) {
    _log->write_status("Created ",
                       _grid_creator->number_of_actual_subgrids() -
                           _grid_creator->number_of_original_subgrids(),
                       " subgrid copies.");
  }

  DistributedPhotonSource< DensitySubGrid > *photon_source = nullptr;
  uint_fast32_t number_of_discrete_photons = 0;
  if (_photon_source_distribution != nullptr) {
    number_of_discrete_photons = _number_of_photons;
  }

  // per thread execution statistics
  std::vector< ThreadStats > thread_stats(_queues.size());

  const uint_fast32_t number_of_continuous_blocks = _queues.size();
  std::vector< ThreadLock > continuous_source_lock(number_of_continuous_blocks);
  AtomicValue< uint_fast32_t > number_of_continuous_photons(0);
  AtomicValue< uint_fast32_t > continuous_photons_flushed(0);
  std::vector< std::vector< PhotonBuffer > > continuous_buffers(
      number_of_continuous_blocks);
  if (_continuous_photon_source != nullptr) {
    for (uint_fast32_t i = 0; i < number_of_continuous_blocks; ++i) {
      continuous_buffers[i].resize(
          _grid_creator->number_of_original_subgrids());
    }
    number_of_continuous_photons.set(_number_of_photons);
  }

  double discrete_photon_weight = 1.;
  double continuous_photon_weight = 1.;
  if (number_of_discrete_photons > 0 &&
      number_of_continuous_photons.value() > 0) {
    number_of_discrete_photons >>= 1;
    number_of_continuous_photons.set(_number_of_photons -
                                     number_of_discrete_photons);
    const double luminosity_ratio =
        _total_luminosity /
            _photon_source_distribution->get_total_luminosity() -
        1.;
    discrete_photon_weight = 2. / (luminosity_ratio + 1.);
    continuous_photon_weight = 2. * luminosity_ratio / (luminosity_ratio + 1.);
  }

  const uint_fast32_t fixed_number_of_continuous_photons =
      number_of_continuous_photons.value();

  if (_photon_source_distribution != nullptr) {
    photon_source = new DistributedPhotonSource< DensitySubGrid >(
        number_of_discrete_photons, *_photon_source_distribution,
        *_grid_creator);
  }

  _time_log.start("subgrid initialisation");
  {
    AtomicValue< size_t > igrid(0);
    start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    while (igrid.value() < _grid_creator->number_of_actual_subgrids()) {
      const size_t this_igrid = igrid.post_increment();
      if (this_igrid < _grid_creator->number_of_actual_subgrids()) {
        DensitySubGrid &subgrid = *_grid_creator->get_subgrid(this_igrid);
        for (int ingb = 0; ingb < TRAVELDIRECTION_NUMBER; ++ingb) {
          subgrid.set_active_buffer(ingb, NEIGHBOUR_OUTSIDE);
          subgrid.set_owning_thread(get_thread_index());
        }
      }
    }
    stop_parallel_timing_block();
  }
  _time_log.end("subgrid initialisation");

  _time_log.start("photoionization loop");
  for (uint_fast32_t iloop = 0; iloop < _number_of_iterations; ++iloop) {

    std::stringstream iloopstr;
    iloopstr << "loop " << iloop;
    _time_log.start(iloopstr.str());
    _memory_log.add_entry(iloopstr.str());

    if (_log) {
      _log->write_status("Starting loop ", iloop, ".");
    }

    uint_fast64_t iteration_start, iteration_end;
    cpucycle_tick(iteration_start);
    _photon_propagation_timer.start();

    // reset the photon source information
    if (photon_source) {
      photon_source->reset();
    }

    // define photon scattering stats
    PhotonPacketStatistics statistics(5);

    if (_trackers != nullptr && iloop == _number_of_iterations - 1) {
      if (_log) {
        _log->write_status("Adding trackers...");
      }
      _trackers->add_trackers(*_grid_creator);
      _number_of_photons =
          std::max(_number_of_photons, _trackers->get_number_of_photons());
      if (_log) {
        _log->write_status("Done adding trackers.");
      }
    }

    // reset mean intensity counters
    {
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < _grid_creator->number_of_actual_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < _grid_creator->number_of_actual_subgrids()) {
          auto gridit = _grid_creator->get_subgrid(this_igrid);
          (*gridit).reset_intensities();
        }
      }
      stop_parallel_timing_block();
    }

    _time_log.start("diffuse field variables");
    // reset the diffuse field variables
    if (_reemission_handler != nullptr) {
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < _grid_creator->number_of_actual_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < _grid_creator->number_of_actual_subgrids()) {
          auto gridit = _grid_creator->get_subgrid(this_igrid);
          for (auto cellit = (*gridit).begin(); cellit != (*gridit).end();
               ++cellit) {
            IonizationVariables &vars = cellit.get_ionization_variables();
            _reemission_handler->set_reemission_probabilities(vars);
          }
        }
      }
      stop_parallel_timing_block();
    }
    _time_log.end("diffuse field variables");

    _time_log.start("photon source tasks");
    size_t number_of_photons_done = 0;
    if (photon_source) {
      while (number_of_photons_done < number_of_discrete_photons) {
        for (size_t isrc = 0; isrc < photon_source->get_number_of_sources();
             ++isrc) {

          const size_t number_of_photons_this_batch =
              photon_source->get_photon_batch(isrc, PHOTONBUFFER_SIZE);
          if (number_of_photons_this_batch > 0) {
            const size_t new_task = _tasks->get_free_element();
            (*_tasks)[new_task].set_type(TASKTYPE_SOURCE_DISCRETE_PHOTON);
            (*_tasks)[new_task].set_subgrid(isrc);
            (*_tasks)[new_task].set_buffer(number_of_photons_this_batch);
            _shared_queue->add_task(new_task);
            number_of_photons_done += number_of_photons_this_batch;
          }
        }
      }
    }
    if (_continuous_photon_source) {
      number_of_continuous_photons.set(fixed_number_of_continuous_photons);
      const uint_fast32_t batch_size = PHOTONBUFFER_SIZE;
      uint_fast32_t block_index = 0;
      const uint_fast32_t num_batches =
          number_of_continuous_photons.value() / batch_size;
      for (uint_fast32_t ibatch = 0; ibatch < num_batches; ++ibatch) {
        const size_t new_task = _tasks->get_free_element();
        (*_tasks)[new_task].set_type(TASKTYPE_SOURCE_CONTINUOUS_PHOTON);
        (*_tasks)[new_task].set_buffer(batch_size);
        (*_tasks)[new_task].set_subgrid(block_index %
                                        number_of_continuous_blocks);
        (*_tasks)[new_task].set_dependency(
            &continuous_source_lock[block_index % number_of_continuous_blocks]);
        ++block_index;
        _shared_queue->add_task(new_task);
        number_of_photons_done += batch_size;
      }
      const uint_fast32_t num_last_batch =
          number_of_continuous_photons.value() % batch_size;
      if (num_last_batch > 0) {
        const size_t new_task = _tasks->get_free_element();
        (*_tasks)[new_task].set_type(TASKTYPE_SOURCE_CONTINUOUS_PHOTON);
        (*_tasks)[new_task].set_buffer(num_last_batch);
        (*_tasks)[new_task].set_subgrid(block_index %
                                        number_of_continuous_blocks);
        (*_tasks)[new_task].set_dependency(
            &continuous_source_lock[block_index % number_of_continuous_blocks]);
        ++block_index;
        _shared_queue->add_task(new_task);
        number_of_photons_done += num_last_batch;
      }
      number_of_continuous_photons.set(fixed_number_of_continuous_photons);
      continuous_photons_flushed.set(0);
    }
    cmac_assert_message(number_of_photons_done == _number_of_photons,
                        "%zu =/= %zu", number_of_photons_done,
                        _number_of_photons);
    _time_log.end("photon source tasks");

    _time_log.start("photon propagation");
    bool global_run_flag = true;
    const uint_fast32_t num_empty_target =
        TRAVELDIRECTION_NUMBER * _grid_creator->number_of_actual_subgrids();
    AtomicValue< uint_fast32_t > num_empty(num_empty_target);
    AtomicValue< uint_fast32_t > num_active_buffers(0);
    AtomicValue< uint_fast32_t > num_photon_done(0);
    Timer verbose_timer;
    verbose_timer.start();
    uint_fast32_t verbose_last_num_empty = 0;
    uint_fast32_t verbose_last_num_active_buffers = 0;
    uint_fast32_t verbose_last_num_photon_done = 0;
    start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
    {
      // thread initialisation
      const int_fast8_t thread_id = get_thread_index();
      PhotonBuffer local_buffers[TRAVELDIRECTION_NUMBER];
      bool local_buffer_flags[TRAVELDIRECTION_NUMBER];
      for (int_fast8_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
        local_buffers[i].set_direction(
            TravelDirections::output_to_input_direction(i));
        local_buffers[i].reset();
        local_buffer_flags[i] = true;
      }

      // actual run flag
      uint_fast32_t current_index = _shared_queue->get_task(*_tasks);
      while (global_run_flag) {

        if (current_index == NO_TASK) {
          uint_fast32_t threshold_size = PHOTONBUFFER_SIZE;
          while (threshold_size > 0) {
            threshold_size >>= 1;
            for (auto gridit = _grid_creator->begin();
                 gridit != _grid_creator->all_end(); ++gridit) {
              DensitySubGrid &this_subgrid = *gridit;
              if (this_subgrid.get_largest_buffer_size() > threshold_size &&
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

                  const size_t task_index = _tasks->get_free_element();
                  Task &new_task = (*_tasks)[task_index];
                  new_task.set_subgrid(
                      (*_buffers)[non_full_index].get_subgrid_index());
                  new_task.set_buffer(non_full_index);
                  if (largest_index > 0) {
                    DensitySubGrid &subgrid = *_grid_creator->get_subgrid(
                        (*_buffers)[non_full_index].get_subgrid_index());
                    new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);

                    // add dependency
                    new_task.set_dependency(subgrid.get_dependency());

                    const uint_fast32_t queue_index =
                        subgrid.get_owning_thread();
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
                    if (this_subgrid.get_active_buffer(ibuffer) !=
                            NEIGHBOUR_OUTSIDE &&
                        (*_buffers)[this_subgrid.get_active_buffer(ibuffer)]
                                .size() > new_largest_size) {
                      new_largest_index = ibuffer;
                      new_largest_size =
                          (*_buffers)[this_subgrid.get_active_buffer(ibuffer)]
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
                  // no semi-full buffers for this subgrid: release the lock
                  // again
                  this_subgrid.get_dependency()->unlock();
                }
              }
            }
          }
          current_index = get_task(thread_id);
        }

        while (current_index != NO_TASK) {

          // execute task
          uint_fast32_t num_tasks_to_add = 0;
          uint_fast32_t tasks_to_add[TRAVELDIRECTION_NUMBER];
          int_fast32_t queues_to_add[TRAVELDIRECTION_NUMBER];

          Task &task = (*_tasks)[current_index];
          thread_stats[thread_id].start(task.get_type());

          if (task.get_type() == TASKTYPE_SOURCE_DISCRETE_PHOTON) {

            task.start(thread_id);
            num_active_buffers.pre_increment();
            const size_t source_index = task.get_subgrid();

            const size_t num_photon_this_loop = task.get_buffer();
            const size_t subgrid_index =
                photon_source->get_subgrid(source_index);

            // get a free photon buffer in the central queue
            uint_fast32_t buffer_index = (*_buffers).get_free_buffer();
            PhotonBuffer &input_buffer = (*_buffers)[buffer_index];

            // set general buffer information
            input_buffer.grow(num_photon_this_loop);
            input_buffer.set_subgrid_index(subgrid_index);
            input_buffer.set_direction(TRAVELDIRECTION_INSIDE);

            const CoordinateVector<> source_position =
                photon_source->get_position(source_index);

            // draw random photons and store them in the buffer
            for (uint_fast32_t i = 0; i < num_photon_this_loop; ++i) {

              PhotonPacket &photon = input_buffer[i];

              photon.set_type(PHOTONTYPE_PRIMARY);
              photon.set_scatter_counter(0);

              // initial position: we currently assume a single source at the
              // origin
              photon.set_position(source_position);

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
              const CoordinateVector<> direction(sint * cosp, sint * sinp,
                                                 cost);

              photon.set_direction(direction);

              // we currently assume equal weight for all photons
              photon.set_weight(discrete_photon_weight);

              // target optical depth (exponential distribution)
              photon.set_target_optical_depth(-std::log(
                  _random_generators[thread_id].get_uniform_random_double()));

              const double frequency =
                  _photon_source_spectrum->get_random_frequency(
                      _random_generators[thread_id]);
              photon.set_energy(frequency);
              for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                double sigma =
                    _cross_sections->get_cross_section(ion, frequency);
#ifndef VARIABLE_ABUNDANCES
                if (ion != ION_H_n) {
                  sigma *= _abundances.get_abundance(get_element(ion));
                }
#endif
                photon.set_photoionization_cross_section(ion, sigma);
              }
            }

            // add to the queue of the corresponding thread
            DensitySubGrid &subgrid =
                *_grid_creator->get_subgrid(subgrid_index);
            const size_t task_index = _tasks->get_free_element();
            Task &new_task = (*_tasks)[task_index];
            new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
            new_task.set_subgrid(subgrid_index);
            new_task.set_buffer(buffer_index);

            // add dependency for task:
            //  - subgrid
            // (the output buffers belong to the subgrid and do not count as a
            // dependency)
            new_task.set_dependency(subgrid.get_dependency());

            queues_to_add[num_tasks_to_add] = subgrid.get_owning_thread();
            tasks_to_add[num_tasks_to_add] = task_index;
            ++num_tasks_to_add;

            // log the end time of the task
            task.stop();

          } else if (task.get_type() == TASKTYPE_SOURCE_CONTINUOUS_PHOTON) {

            task.start(thread_id);

            const uint_fast32_t source_copy = task.get_subgrid();
            const size_t num_photon_this_loop = task.get_buffer();

            // draw random photons and store them in the continuous buffers
            for (uint_fast32_t i = 0; i < num_photon_this_loop; ++i) {

              auto posdir =
                  _continuous_photon_source->get_random_incoming_direction(
                      _random_generators[thread_id]);

              const size_t subgrid_index =
                  _grid_creator->get_subgrid(posdir.first).get_index();

              PhotonBuffer &active_buffer =
                  continuous_buffers[source_copy][subgrid_index];
              const uint_fast32_t active_index =
                  active_buffer.get_next_free_photon();

              PhotonPacket &photon = active_buffer[active_index];

              photon.set_type(PHOTONTYPE_PRIMARY);
              photon.set_scatter_counter(0);

              // initial position: we currently assume a single source at the
              // origin
              photon.set_position(posdir.first);
              photon.set_direction(posdir.second);

              // we currently assume equal weight for all photons
              photon.set_weight(continuous_photon_weight);

              // target optical depth (exponential distribution)
              photon.set_target_optical_depth(-std::log(
                  _random_generators[thread_id].get_uniform_random_double()));

              const double frequency =
                  _continuous_photon_source_spectrum->get_random_frequency(
                      _random_generators[thread_id]);
              photon.set_energy(frequency);
              for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                double sigma =
                    _cross_sections->get_cross_section(ion, frequency);
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
                uint_fast32_t buffer_index = (*_buffers).get_free_buffer();
                PhotonBuffer &input_buffer = (*_buffers)[buffer_index];

                num_active_buffers.pre_increment();

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
                DensitySubGrid &subgrid =
                    *_grid_creator->get_subgrid(subgrid_index);
                const size_t task_index = _tasks->get_free_element();
                Task &new_task = (*_tasks)[task_index];
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

            cmac_assert(number_of_continuous_photons.value() >=
                        num_photon_this_loop);
            number_of_continuous_photons.pre_subtract(num_photon_this_loop);

            // are we done with all photons?
            if (number_of_continuous_photons.value() == 0) {
              const uint_fast32_t flush =
                  continuous_photons_flushed.pre_increment();
              if (flush == 1) {
                // yes: add tasks to flush all buffers
                for (uint_fast32_t other_copy = 0;
                     other_copy < number_of_continuous_blocks; ++other_copy) {
                  const size_t task_index = _tasks->get_free_element();
                  Task &new_task = (*_tasks)[task_index];
                  new_task.set_type(TASKTYPE_FLUSH_CONTINUOUS_PHOTON_BUFFERS);
                  new_task.set_subgrid(other_copy);
                  new_task.set_buffer(0);

                  // add dependency for task:
                  //  - subgrid
                  // (the output buffers belong to the subgrid and do not
                  // count as a dependency)
                  new_task.set_dependency(&continuous_source_lock[other_copy]);

                  _shared_queue->add_task(task_index);
                }
              }
            }

            // log the end time of the task
            task.stop();

          } else if (task.get_type() ==
                     TASKTYPE_FLUSH_CONTINUOUS_PHOTON_BUFFERS) {

            task.start(thread_id);

            const uint_fast32_t source_copy = task.get_subgrid();
            for (uint_fast32_t i = 0;
                 i < continuous_buffers[source_copy].size(); ++i) {
              if (continuous_buffers[source_copy][i].size() > 0) {
                // yes: send the buffer off!
                uint_fast32_t buffer_index = (*_buffers).get_free_buffer();
                PhotonBuffer &input_buffer = (*_buffers)[buffer_index];

                num_active_buffers.pre_increment();

                // set general buffer information
                input_buffer.grow(continuous_buffers[source_copy][i].size());
                input_buffer.set_subgrid_index(i);
                input_buffer.set_direction(TRAVELDIRECTION_INSIDE);

                // copy over the photons
                for (uint_fast32_t iphoton = 0;
                     iphoton < continuous_buffers[source_copy][i].size();
                     ++iphoton) {
                  input_buffer[iphoton] =
                      continuous_buffers[source_copy][i][iphoton];
                }

                // reset the active buffer
                continuous_buffers[source_copy][i].reset();

                // add to the queue of the corresponding thread
                DensitySubGrid &subgrid = *_grid_creator->get_subgrid(i);
                const size_t task_index = _tasks->get_free_element();
                Task &new_task = (*_tasks)[task_index];
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

            // log the end time of the task
            task.stop();

          } else if (task.get_type() == TASKTYPE_PHOTON_REEMIT) {

            task.start(thread_id);

            const size_t current_buffer_index = task.get_buffer();
            PhotonBuffer &buffer = (*_buffers)[current_buffer_index];

            uint_fast32_t num_photon_done_now = buffer.size();
            DensitySubGrid &subgrid =
                *_grid_creator->get_subgrid(task.get_subgrid());

            // reemission
            uint_fast32_t index = 0;
            for (uint_fast32_t iphoton = 0; iphoton < buffer.size();
                 ++iphoton) {
              PhotonPacket &old_photon = buffer[iphoton];
              const IonizationVariables &ionization_variables =
                  subgrid.get_cell(old_photon.get_position())
                      .get_ionization_variables();
#ifdef HAS_HELIUM
#ifdef VARIABLE_ABUNDANCES
              const double AHe =
                  ionization_variables.get_abundances().get_abundance(
                      ELEMENT_He);
#else
              // the helium abundance is already part of the cross section
              const double AHe = 1.;
#endif
#else
              const double AHe = 0.;
#endif
              PhotonType new_type;
              const double new_frequency = _reemission_handler->reemit(
                  old_photon, AHe, ionization_variables,
                  _random_generators[thread_id], new_type);
              if (new_frequency > 0.) {
                PhotonPacket &new_photon = buffer[index];
                new_photon.set_type(new_type);
                new_photon.set_scatter_counter(
                    old_photon.get_scatter_counter() + 1);
                new_photon.set_position(old_photon.get_position());
                new_photon.set_weight(old_photon.get_weight());

                new_photon.set_energy(new_frequency);
                for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                  double sigma =
                      _cross_sections->get_cross_section(ion, new_frequency);
#ifndef VARIABLE_ABUNDANCES
                  if (ion != ION_H_n) {
                    sigma *= _abundances.get_abundance(get_element(ion));
                  }
#endif
                  new_photon.set_photoionization_cross_section(ion, sigma);
                }

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
                const CoordinateVector<> direction(sint * cosp, sint * sinp,
                                                   cost);

                new_photon.set_direction(direction);

                // target optical depth (exponential distribution)
                new_photon.set_target_optical_depth(-std::log(
                    _random_generators[thread_id].get_uniform_random_double()));

                ++index;
              } else {
                statistics.absorb_photon(old_photon);
              }
            }
            // update the size of the buffer to account for photons that were
            // not reemitted
            buffer.grow(index);

            num_photon_done_now -= buffer.size();
            num_photon_done.pre_add(num_photon_done_now);

            if (index > 0) {
              // there are still photon packets left: generate a traversal
              // task
              const size_t task_index = _tasks->get_free_element();
              Task &new_task = (*_tasks)[task_index];
              new_task.set_type(TASKTYPE_PHOTON_TRAVERSAL);
              new_task.set_subgrid(task.get_subgrid());
              new_task.set_buffer(current_buffer_index);
              new_task.set_dependency(subgrid.get_dependency());

              queues_to_add[num_tasks_to_add] = subgrid.get_owning_thread();
              tasks_to_add[num_tasks_to_add] = task_index;
              ++num_tasks_to_add;
            } else {
              // delete the original buffer, as we are done with it
              _buffers->free_buffer(current_buffer_index);
              cmac_assert_message(num_active_buffers.value() > 0,
                                  "Number of active buffers < 0!");
              num_active_buffers.pre_decrement();
            }

            task.stop();

          } else if (task.get_type() == TASKTYPE_PHOTON_TRAVERSAL) {

            task.start(thread_id);
            uint_fast64_t task_start, task_stop;
            cpucycle_tick(task_start);

            const uint_fast32_t current_buffer_index = task.get_buffer();
            PhotonBuffer &photon_buffer = (*_buffers)[current_buffer_index];
            const uint_fast32_t igrid = photon_buffer.get_subgrid_index();
            DensitySubGrid &this_grid = *_grid_creator->get_subgrid(igrid);

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
            if (_reemission_handler == nullptr) {
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
              } else {
                if (result == 0) {
                  statistics.absorb_photon(photon);
                } else {
                  statistics.escape_photon(photon);
                }
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
                  if ((*_buffers)[add_index].size() == 0) {
                    _buffers->free_buffer(add_index);
                    this_grid.set_active_buffer(i, NEIGHBOUR_OUTSIDE);
                    // we have created a new empty buffer
                    num_empty.pre_increment();
                  } else {
                    this_grid.set_active_buffer(i, add_index);

                    cmac_assert_message(
                        (*_buffers)[add_index].get_subgrid_index() == ngb,
                        "Wrong subgrid");
                    cmac_assert_message(
                        (*_buffers)[add_index].get_direction() ==
                            TravelDirections::output_to_input_direction(i),
                        "Wrong direction");
                  }

                  // YES: create a task for the buffer and add it to the queue
                  // the task type depends on the buffer: photon packets in the
                  // internal buffer were absorbed and could be reemitted,
                  // photon packets in the other buffers left the subgrid and
                  // need to be traversed in the neighbouring subgrid
                  if (i > 0) {
                    DensitySubGrid &subgrid = *_grid_creator->get_subgrid(
                        (*_buffers)[new_index].get_subgrid_index());
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
                        (*_grid_creator->get_subgrid(ngb)).get_owning_thread();
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
            cpucycle_tick(task_stop);
            this_grid.add_computational_cost(task_stop - task_start);
          }
          task.unlock_dependency();
          thread_stats[thread_id].stop(task.get_type());

          // we are done with the task, clean up (if we don't output it)
          if (!_task_plot) {
            _tasks->free_element(current_index);
          }

          for (uint_fast32_t itask = 0; itask < num_tasks_to_add; ++itask) {
            if (queues_to_add[itask] < 0) {
              // general queue
              _shared_queue->add_task(tasks_to_add[itask]);
            } else {
              _queues[queues_to_add[itask]]->add_task(tasks_to_add[itask]);
            }
          }

          // we need to call task_status twice: once inside the task loop
          // (for when thread 0 is happily working) and once outside (for
          // when thread 0 is idling)
          task_status(_verbose, _log, thread_id, verbose_timer, num_empty,
                      num_empty_target, num_active_buffers, num_photon_done,
                      _number_of_photons, verbose_last_num_empty,
                      verbose_last_num_active_buffers,
                      verbose_last_num_photon_done, *_shared_queue, _queues,
                      *_tasks, *_grid_creator);

          current_index = get_task(thread_id);
        }

        task_status(_verbose, _log, thread_id, verbose_timer, num_empty,
                    num_empty_target, num_active_buffers, num_photon_done,
                    _number_of_photons, verbose_last_num_empty,
                    verbose_last_num_active_buffers,
                    verbose_last_num_photon_done, *_shared_queue, _queues,
                    *_tasks, *_grid_creator);

#ifdef OUTPUT_STOP_CONDITION
        cmac_warning("num_empty: %" PRIuFAST32 " (%" PRIuFAST32
                     "), num_active_buffers: %" PRIuFAST32
                     ", num_photon_done: %" PRIuFAST32,
                     num_empty.value(), num_empty_target,
                     num_active_buffers.value(), num_photon_done.value());
#endif
        if (num_empty.value() == num_empty_target &&
            num_active_buffers.value() == 0 &&
            num_photon_done.value() == _number_of_photons) {
          global_run_flag = false;
        } else {
          current_index = get_task(thread_id);
        }
      } // while(global_run_flag)
    }   // parallel region
    stop_parallel_timing_block();
    _time_log.end("photon propagation");

    _time_log.start("update copies");
    start_parallel_timing_block();
    _grid_creator->update_original_counters();
    stop_parallel_timing_block();
    _time_log.end("update copies");
    if (iloop == _number_of_iterations - 1) {
      statistics.print_stats();
    }
    _photon_propagation_timer.stop();

    if (_log != nullptr) {
      _log->write_info("Done shooting photons.");
      _log->write_info("Starting temperature calculation...");
    }
    _cell_update_timer.start();
    _time_log.start("temperature calculation");
    {
      AtomicValue< size_t > igrid(0);
      start_parallel_timing_block();
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#endif
      while (igrid.value() < _grid_creator->number_of_original_subgrids()) {
        const size_t this_igrid = igrid.post_increment();
        if (this_igrid < _grid_creator->number_of_original_subgrids()) {
          auto gridit = _grid_creator->get_subgrid(this_igrid);

          const size_t itask = _tasks->get_free_element();
          Task &task = (*_tasks)[itask];
          task.set_type(TASKTYPE_TEMPERATURE_STATE);
          task.start(get_thread_index());
          thread_stats[get_thread_index()].start(TASKTYPE_TEMPERATURE_STATE);

#ifndef VARIABLE_ABUNDANCES
          // correct the intensity counters for abundance factors
          for (auto cellit = (*gridit).begin(); cellit != (*gridit).end();
               ++cellit) {
            IonizationVariables &vars = cellit.get_ionization_variables();
            for (int_fast32_t ion = 1; ion < NUMBER_OF_IONNAMES; ++ion) {
              const double abundance =
                  _abundances.get_abundance(get_element(ion));
              if (abundance > 0.) {
                vars.set_mean_intensity(ion, vars.get_mean_intensity(ion) /
                                                 abundance);
              }
            }
#ifdef HAS_HELIUM
            vars.set_heating(HEATINGTERM_He,
                             vars.get_heating(HEATINGTERM_He) /
                                 _abundances.get_abundance(ELEMENT_He));
#endif
          }
#endif
          _temperature_calculator->calculate_temperature(
              iloop, _number_of_photons, *gridit);
          task.stop();
          thread_stats[get_thread_index()].stop(TASKTYPE_TEMPERATURE_STATE);

          // clean up (if we don't need the task any more)
          if (!_task_plot) {
            _tasks->free_element(itask);
          }
        }
      }
      stop_parallel_timing_block();
    }
    _time_log.end("temperature calculation");

    _cell_update_timer.stop();

    // output diagnostic information
    {
      uint_fast64_t early_iteration_end;
      cpucycle_tick(early_iteration_end);
      // compose the file name
      std::stringstream filename;
      filename << "diagnostics_";
      filename.fill('0');
      filename.width(2);
      filename << iloop;
      filename << ".txt";

      // now open the file
      std::ofstream ofile(filename.str(), std::ofstream::trunc);
      ofile << "iteration:\n";
      ofile << "  start: " << iteration_start << "\n";
      ofile << "  end: " << early_iteration_end << "\n";
      ofile << "buffers:\n";
      ofile << "  total: " << _buffers->get_total_number_elements() << "\n";
      _buffers->reset_total_number_elements();
      ofile << "  max: " << _buffers->get_max_number_elements() << "\n";
      _buffers->reset_max_number_elements();
      ofile << "tasks:\n";
      ofile << "  total: " << _tasks->get_total_number_taken() << "\n";
      _tasks->reset_total_number_taken();
      ofile << "  max: " << _tasks->get_max_number_taken() << "\n";
      _tasks->reset_max_number_taken();
      ofile << "queues:\n";
      ofile << "  shared:\n";
      ofile << "    total: " << _shared_queue->get_total_queue_size() << "\n";
      _shared_queue->reset_total_queue_size();
      ofile << "    max: " << _shared_queue->get_max_queue_size() << "\n";
      _shared_queue->reset_max_queue_size();
      ofile << "    average: " << _shared_queue->get_average_queue_size()
            << "\n";
      _shared_queue->reset_average_queue_size();
      for (uint_fast32_t i = 0; i < _queues.size(); ++i) {
        ofile << "  thread " << i << ":\n";
        ofile << "    total: " << _queues[i]->get_total_queue_size() << "\n";
        _queues[i]->reset_total_queue_size();
        ofile << "    max: " << _queues[i]->get_max_queue_size() << "\n";
        _queues[i]->reset_max_queue_size();
        ofile << "    average: " << _queues[i]->get_average_queue_size()
              << "\n";
        _queues[i]->reset_average_queue_size();
      }
      ofile << "threads:\n";
      for (uint_fast32_t i = 0; i < thread_stats.size(); ++i) {
        ofile << "  thread " << i << ":\n";
        for (int_fast32_t j = 0; j < TASKTYPE_NUMBER; ++j) {
          ofile << "    task " << j << ":\n";
          ofile << "      number: "
                << thread_stats[i].get_number_of_tasks_executed(j) << "\n";
          ofile << "      time: " << thread_stats[i].get_total_time(j) << "\n";
          ofile << "      squared time: "
                << thread_stats[i].get_total_time_squared(j) << "\n";
        }
        thread_stats[i].reset();
      }
      ofile << "subgrids:\n";
      for (auto it = _grid_creator->begin(); it != _grid_creator->all_end();
           ++it) {
        ofile << "  subgrid " << it.get_index() << ": "
              << (*it).get_computational_cost() << "\n";
        (*it).reset_computational_cost();
      }
    }

    _time_log.start("buffer reset");
    _buffers->reset();
    _time_log.end("buffer reset");

    cmac_assert_message(_buffers->is_empty(), "Number of active buffers: %zu",
                        _buffers->get_number_of_active_buffers());

    // update copies
    _time_log.start("copy update");
    start_parallel_timing_block();
    _grid_creator->update_copy_properties();
    stop_parallel_timing_block();
    _time_log.end("copy update");

    if (_task_plot) {
      _time_log.start("task output");
      cpucycle_tick(iteration_end);
      output_tasks(iloop, *_tasks, iteration_start, iteration_end);
      output_queues(iloop, _queues, *_shared_queue);
      _time_log.end("task output");
    }

    _time_log.start("task reset");
    _tasks->clear();
    _time_log.end("task reset");

    _time_log.end(iloopstr.str());

  } // photoionization loop
  _time_log.end("photoionization loop");

  if (photon_source) {
    delete photon_source;
  }

  if (_trackers != nullptr) {
    if (_log) {
      _log->write_status("Outputting trackers...");
    }
    _trackers->output_trackers();
    if (_log) {
      _log->write_status("Done outputting trackers.");
    }
  }

  _time_log.start("snapshot");
  _density_grid_writer->write(*_grid_creator, _number_of_iterations,
                              _parameter_file);
  _time_log.end("snapshot");

  _time_log.end("main run");
}
