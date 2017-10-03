/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file IonizationSimulation.cpp
 *
 * @brief IonizationSimulation implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "IonizationSimulation.hpp"
#include "ChargeTransferRates.hpp"
#include "ContinuousPhotonSourceFactory.hpp"
#include "CrossSectionsFactory.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DensityMaskFactory.hpp"
#include "DiffuseReemissionHandler.hpp"
#include "IonizationVariablesPropertyAccessors.hpp"
#include "LineCoolingData.hpp"
#include "MPICommunicator.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "RecombinationRatesFactory.hpp"
#include "SimulationBox.hpp"
#include "TemperatureCalculator.hpp"
#include "WorkEnvironment.hpp"
#include <fstream>

/*! @brief Start the serial and total program time timers at the start of a
 *  member function call. */
#define function_start_timers()                                                \
  _serial_timer.start();                                                       \
  _total_timer.start();

/*! @brief Stop the serial and total program time timers at the end of a
 *  member function call. */
#define function_stop_timers()                                                 \
  _serial_timer.stop();                                                        \
  _total_timer.stop();

/*! @brief Stop the serial time timer and start the parallel time timer. */
#define start_parallel_timing_block()                                          \
  _serial_timer.stop();                                                        \
  _parallel_timer.start();

/*! @brief Stop the parallel time timer and start the serial time timer. */
#define stop_parallel_timing_block()                                           \
  _parallel_timer.stop();                                                      \
  _serial_timer.start();

/**
 * @brief Constructor.
 *
 * This method will read the following parameters from the parameter file:
 *  - number of iterations: Number of iterations of the photoionization
 *    algorithm to perform (default: 10)
 *  - number of photons: Number of photon packets to use for each iteration of
 *    the photoionization algorithm (default: 1e5)
 *  - number of photons first loop: Number of photon packets to use for the
 *    first iteration of the photoionization algorithm (default: (number of
 *    photons))
 *  - output folder: Folder where all output files will be placed (default: .)
 *  - random seed: Seed used to initialize the random number generator (default:
 *    42)
 *
 * @param write_output Should this process write output?
 * @param every_iteration_output Write an output file after every iteration of
 * the algorithm?
 * @param num_thread Number of shared memory parallel threads to use.
 * @param parameterfile Name of the parameter file to use.
 * @param mpi_communicator MPICommunicator to use for distributed memory
 * communications.
 * @param log Log to write logging info to.
 */
IonizationSimulation::IonizationSimulation(const bool write_output,
                                           const bool every_iteration_output,
                                           const int num_thread,
                                           const std::string parameterfile,
                                           MPICommunicator *mpi_communicator,
                                           Log *log)
    : _num_thread(WorkEnvironment::set_max_num_threads(num_thread)),
      _every_iteration_output(every_iteration_output),
      _mpi_communicator(mpi_communicator), _log(log),
      _work_distributor(_num_thread), _parameter_file(parameterfile),
      _number_of_iterations(_parameter_file.get_value< uint_fast32_t >(
          "IonizationSimulation:number of iterations", 10)),
      _number_of_photons(_parameter_file.get_value< uint_fast32_t >(
          "IonizationSimulation:number of photons", 1e5)),
      _number_of_photons_init(_parameter_file.get_value< uint_fast32_t >(
          "IonizationSimulation:number of photons first loop",
          _number_of_photons)),
      _abundances(_parameter_file, _log) {

  function_start_timers();

  if (_log) {
    if (_num_thread == 1) {
      _log->write_status("IonizationSimulation will use 1 thread.");
    } else {
      _log->write_status("IonizationSimulation will use ", _num_thread,
                         " threads.");
    }
  }

  // create cross section and recombination rate objects
  _cross_sections = CrossSectionsFactory::generate(_parameter_file, _log);
  _recombination_rates =
      RecombinationRatesFactory::generate(_parameter_file, _log);

  // create the density grid and related objects
  _density_function = DensityFunctionFactory::generate(_parameter_file, _log);
  _density_mask = DensityMaskFactory::generate(_parameter_file, _log);

  const SimulationBox simulation_box(_parameter_file);
  _density_grid = DensityGridFactory::generate(simulation_box, _parameter_file,
                                               false, _log);

  // create the discrete UV sources
  _photon_source_distribution =
      PhotonSourceDistributionFactory::generate(_parameter_file, _log);
  _photon_source_spectrum = PhotonSourceSpectrumFactory::generate(
      "PhotonSourceSpectrum", _parameter_file, _log);

  // sanity checks on discrete sources
  if (_photon_source_distribution != nullptr &&
      _photon_source_spectrum == nullptr) {
    cmac_error("No spectrum provided for the discrete photon sources!");
  }

  // create the continuous UV sources
  _continuous_photon_source = ContinuousPhotonSourceFactory::generate(
      simulation_box.get_box(), _parameter_file, _log);
  _continuous_photon_source_spectrum = PhotonSourceSpectrumFactory::generate(
      "ContinuousPhotonSourceSpectrum", _parameter_file, _log);

  // sanity checks on continuous sources
  if (_continuous_photon_source != nullptr &&
      _continuous_photon_source_spectrum == nullptr) {
    cmac_error("No spectrum provided for the continuous photon sources!");
  }

  // create the actual photon source objects that emits the UV photons
  _photon_source = new PhotonSource(
      _photon_source_distribution, _photon_source_spectrum,
      _continuous_photon_source, _continuous_photon_source_spectrum,
      _abundances, *_cross_sections, _parameter_file, _log);
  const double total_luminosity = _photon_source->get_total_luminosity();

  // set up output
  std::string output_folder =
      Utilities::get_absolute_path(_parameter_file.get_value< std::string >(
          "IonizationSimulation:output folder", "."));
  _density_grid_writer = nullptr;
  if (write_output) {
    _density_grid_writer = DensityGridWriterFactory::generate(
        output_folder, _parameter_file, _log);
  }

  // used to calculate both the ionization state and the temperature
  _temperature_calculator = new TemperatureCalculator(
      total_luminosity, _abundances, _line_cooling_data, *_recombination_rates,
      _charge_transfer_rates, _parameter_file, _log);

  // create ray tracing objects
  int_fast32_t random_seed = _parameter_file.get_value< int_fast32_t >(
      "IonizationSimulation:random seed", 42);
  // make sure every thread on every process has another random seed
  if (_mpi_communicator) {
    random_seed += _mpi_communicator->get_rank() * _num_thread;
  }
  _ionization_photon_shoot_job_market = new IonizationPhotonShootJobMarket(
      *_photon_source, random_seed, *_density_grid, 0, 100, _num_thread);

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output) {
    std::ofstream pfile(output_folder + "/parameters-usedvalues.param");
    _parameter_file.print_contents(pfile);
    pfile.close();
    if (_log) {
      _log->write_status("Wrote used parameters to ", output_folder,
                         "/parameters-usedvalues.param.");
    }
  }

  function_stop_timers();
}

/**
 * @brief Initialize the simulation.
 *
 * @param density_function DensityFunction to use. If no DensityFunction is
 * given, the internal DensityFunction is used.
 */
void IonizationSimulation::initialize(DensityFunction *density_function) {

  function_start_timers();

  if (density_function == nullptr) {
    density_function = _density_function;
  }

  // initialize the density function
  if (_log) {
    _log->write_status("Initializing DensityFunction...");
  }
  density_function->initialize();
  if (_log) {
    _log->write_status("Done.");
  }

  // initialize the actual grid
  std::pair< unsigned long, unsigned long > block;
  if (_mpi_communicator) {
    block = _mpi_communicator->distribute_block(
        0, _density_grid->get_number_of_cells());
  } else {
    block = std::make_pair(0, _density_grid->get_number_of_cells());
  }

  start_parallel_timing_block();
  _density_grid->initialize(block, *density_function);
  stop_parallel_timing_block();

  // _density_grid->initialize initialized:
  // - densities
  // - temperatures
  // - neutral fractions of hydrogen and helium
  // we have to gather these across all processes

  if (_mpi_communicator) {
    start_parallel_timing_block();
    std::pair< DensityGrid::iterator, DensityGrid::iterator > local_chunk =
        _density_grid->get_chunk(block.first, block.second);
    _mpi_communicator->gather< double, NumberDensityPropertyAccessor >(
        _density_grid->begin(), _density_grid->end(), local_chunk.first,
        local_chunk.second, 0);
    _mpi_communicator->gather< double, TemperaturePropertyAccessor >(
        _density_grid->begin(), _density_grid->end(), local_chunk.first,
        local_chunk.second, 0);
    _mpi_communicator
        ->gather< double, IonicFractionPropertyAccessor< ION_H_n > >(
            _density_grid->begin(), _density_grid->end(), local_chunk.first,
            local_chunk.second, 0);
    _mpi_communicator
        ->gather< double, IonicFractionPropertyAccessor< ION_He_n > >(
            _density_grid->begin(), _density_grid->end(), local_chunk.first,
            local_chunk.second, 0);
    stop_parallel_timing_block();
  }

  // if necessary, initialize and apply the density mask
  if (_density_mask != nullptr) {
    if (_log) {
      _log->write_status("Initializing DensityMask...");
    }
    _density_mask->initialize();
    if (_log) {
      _log->write_status("Done initializing mask. Applying mask...");
    }
    _density_mask->apply(*_density_grid);
    if (_log) {
      _log->write_status("Done applying mask.");
    }
  }

  function_stop_timers();
}

/**
 * @brief Run the actual simulation.
 *
 * @param density_grid_writer DensityGridWriter to use for the final output. If
 * an internal DensityGridWriter exists, this one will write more output.
 */
void IonizationSimulation::run(DensityGridWriter *density_grid_writer) {

  function_start_timers();

  // write the initial state of the grid to an output file
  if (_density_grid_writer) {
    _density_grid_writer->write(*_density_grid, 0, _parameter_file);
  }

  std::pair< unsigned long, unsigned long > block;
  if (_mpi_communicator) {
    block = _mpi_communicator->distribute_block(
        0, _density_grid->get_number_of_cells());
  } else {
    block = std::make_pair(0, _density_grid->get_number_of_cells());
  }
  std::pair< DensityGrid::iterator, DensityGrid::iterator > local_chunk =
      _density_grid->get_chunk(block.first, block.second);

  // finally: the actual program loop whereby the density grid is ray traced
  // using photon packets generated by the stellar sources
  unsigned int loop = 0;
  while (loop < _number_of_iterations) {

    if (_log) {
      _log->write_status("Starting loop ", loop, ".");
    }

    unsigned int lnumphoton = _number_of_photons;

    if (loop == 0) {
      // overwrite the number of photons for the first loop (might be useful
      // if more than 1 boundary is periodic, since the initial neutral
      // fractions are very low)
      lnumphoton = _number_of_photons_init;
    }

    _density_grid->reset_grid(*_density_function);
    DiffuseReemissionHandler::set_reemission_probabilities(*_density_grid);
    if (_log) {
      _log->write_status("Start shooting ", lnumphoton, " photons...");
    }

    double typecount[PHOTONTYPE_NUMBER] = {0};

    double totweight = 0.;

    unsigned int local_numphoton = lnumphoton;

    // make sure this process does only part of the total number of photons
    if (_mpi_communicator) {
      local_numphoton = _mpi_communicator->distribute(local_numphoton);
    }

    _ionization_photon_shoot_job_market->set_numphoton(local_numphoton);
    _work_timer.start();
    start_parallel_timing_block();
    _work_distributor.do_in_parallel(*_ionization_photon_shoot_job_market);
    stop_parallel_timing_block();
    _work_timer.stop();

    _ionization_photon_shoot_job_market->update_counters(totweight, typecount);

    // make sure the total weight and typecount is reduced across all
    // processes
    if (_mpi_communicator) {
      start_parallel_timing_block();
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES >(totweight);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, PHOTONTYPE_NUMBER >(
          typecount);
      stop_parallel_timing_block();
    }

    if (_log) {
      _log->write_status("Done shooting photons.");
      _log->write_status(
          100. * typecount[PHOTONTYPE_ABSORBED] / totweight,
          "% of photons were reemitted as non-ionizing photons.");
      _log->write_status(100. * (typecount[PHOTONTYPE_DIFFUSE_HI] +
                                 typecount[PHOTONTYPE_DIFFUSE_HeI]) /
                             totweight,
                         "% of photons were scattered.");
      double escape_fraction =
          (100. * (totweight - typecount[PHOTONTYPE_ABSORBED])) / totweight;
      // since totweight is updated in chunks, while the counters are updated
      // per photon, round off might cause totweight to be slightly smaller
      // than the counter value. This gives (strange looking) negative escape
      // fractions, which we reset to 0 here.
      escape_fraction = std::max(0., escape_fraction);
      _log->write_status("Escape fraction: ", escape_fraction, "%.");
      double escape_fraction_HI =
          (100. * typecount[PHOTONTYPE_DIFFUSE_HI]) / totweight;
      _log->write_status("Diffuse HI escape fraction: ", escape_fraction_HI,
                         "%.");
      double escape_fraction_HeI =
          (100. * typecount[PHOTONTYPE_DIFFUSE_HeI]) / totweight;
      _log->write_status("Diffuse HeI escape fraction: ", escape_fraction_HeI,
                         "%.");
    }

    if (_log) {
      _log->write_status("Calculating ionization state after shooting ",
                         lnumphoton, " photons...");
    }

    // reduce the mean intensity integrals and heating terms across all
    // processes
    start_parallel_timing_block();

    if (_mpi_communicator) {
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_H_n > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_He_n > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_C_p1 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_C_p2 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_N_n > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_N_p1 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_N_p2 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_O_n > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_O_p1 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_Ne_n > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_Ne_p1 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_S_p1 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_S_p2 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 MeanIntensityPropertyAccessor< ION_S_p3 > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 HeatingPropertyAccessor< HEATINGTERM_H > >(
          _density_grid->begin(), _density_grid->end(), 0);
      _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, double,
                                 HeatingPropertyAccessor< HEATINGTERM_He > >(
          _density_grid->begin(), _density_grid->end(), 0);
    }

    _temperature_calculator->calculate_temperature(loop, totweight,
                                                   *_density_grid, block);

    // the calculation above will have changed the ionic fractions, and might
    // have changed the temperatures
    // we have to gather these across all processes

    if (_mpi_communicator) {
      _mpi_communicator->gather< double, TemperaturePropertyAccessor >(
          _density_grid->begin(), _density_grid->end(), local_chunk.first,
          local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_H_n > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_He_n > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_C_p1 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_C_p2 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_N_n > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_N_p1 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_N_p2 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_O_n > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_O_p1 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_Ne_n > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_Ne_p1 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_S_p1 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_S_p2 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
      _mpi_communicator
          ->gather< double, IonicFractionPropertyAccessor< ION_S_p3 > >(
              _density_grid->begin(), _density_grid->end(), local_chunk.first,
              local_chunk.second, 0);
    }

    stop_parallel_timing_block();

    if (_log) {
      _log->write_status("Done calculating ionization state.");
    }

    // calculate emissivities
    // we disabled this, since we now have the post-processing Python library
    // for this
    //    if (loop > 3 && abundances.get_abundance(ELEMENT_He) > 0.) {
    //      emissivity_calculator.calculate_emissivities(*grid);
    //    }

    ++loop;

    if (_density_grid_writer && _every_iteration_output &&
        loop < _number_of_iterations) {
      _density_grid_writer->write(*_density_grid, loop, _parameter_file);
    }
  }

  if (_log && loop == _number_of_iterations) {
    _log->write_status("Maximum number of iterations (", _number_of_iterations,
                       ") reached, stopping.");
  }

  // write final snapshot
  if (_density_grid_writer) {
    _density_grid_writer->write(*_density_grid, _number_of_iterations,
                                _parameter_file);
  }
  if (density_grid_writer) {
    density_grid_writer->write(*_density_grid, _number_of_iterations,
                               _parameter_file);
  }

  if (_log) {
    _log->write_status("Total photon shooting time: ",
                       Utilities::human_readable_time(_work_timer.value()),
                       ".");
  }

  function_stop_timers();
}

/**
 * @brief Destructor.
 *
 * Deallocate memory used by internal variables.
 */
IonizationSimulation::~IonizationSimulation() {

  if (_log) {
    _log->write_status("Total serial time: ",
                       Utilities::human_readable_time(_serial_timer.value()));
    _log->write_status("Total parallel time: ",
                       Utilities::human_readable_time(_parallel_timer.value()));
    _log->write_status("Total overall time: ",
                       Utilities::human_readable_time(_total_timer.value()));
  }

  // we delete the objects in the opposite order in which they were created
  // note that we do not check for nullptrs, as deleting a nullptr is allowed
  // and won't do anything

  // ray tracing objects
  delete _ionization_photon_shoot_job_market;

  // computation objects
  delete _temperature_calculator;

  // snapshot output
  delete _density_grid_writer;

  // actual photon source object
  delete _photon_source;

  // continuous sources
  delete _continuous_photon_source_spectrum;
  delete _continuous_photon_source;

  // discrete sources
  delete _photon_source_spectrum;
  delete _photon_source_distribution;

  // density grid and related objects
  delete _density_grid;
  delete _density_mask;
  delete _density_function;

  // cross sections and recombination rates
  delete _cross_sections;
  delete _recombination_rates;
}
