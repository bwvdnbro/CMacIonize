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
 * @file RadiationHydrodynamicsSimulation.cpp
 *
 * @brief RadiationHydrodynamicsSimulation implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "RadiationHydrodynamicsSimulation.hpp"
#include "Abundances.hpp"
#include "Box.hpp"
#include "ChargeTransferRates.hpp"
#include "CommandLineOption.hpp"
#include "CommandLineParser.hpp"
#include "CompilerInfo.hpp"
#include "Configuration.hpp"
#include "ConfigurationInfo.hpp"
#include "ContinuousPhotonSourceFactory.hpp"
#include "CoordinateVector.hpp"
#include "CrossSectionsFactory.hpp"
#include "DeRijckeRadiativeCooling.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DensityMaskFactory.hpp"
#include "DiffuseReemissionHandler.hpp"
#include "DustSimulation.hpp"
#include "EmissivityCalculator.hpp"
#include "ExternalPotentialFactory.hpp"
#include "FileLog.hpp"
#include "HydroIntegrator.hpp"
#include "HydroMaskFactory.hpp"
#include "IonizationPhotonShootJobMarket.hpp"
#include "IonizationSimulation.hpp"
#include "LineCoolingData.hpp"
#include "MPICommunicator.hpp"
#include "ParameterFile.hpp"
#include "PhotonSource.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "RecombinationRatesFactory.hpp"
#include "SimulationBox.hpp"
#include "TemperatureCalculator.hpp"
#include "TerminalLog.hpp"
#include "TimeLine.hpp"
#include "Timer.hpp"
#include "TreeSelfGravity.hpp"
#include "WorkDistributor.hpp"
#include "WorkEnvironment.hpp"

#include <iostream>
#include <string>

/*! @brief Stop the serial time timer and start the parallel time timer. */
#define start_parallel_timing_block()                                          \
  serial_timer.stop();                                                         \
  parallel_timer.start();

/*! @brief Stop the parallel time timer and start the serial time timer. */
#define stop_parallel_timing_block()                                           \
  parallel_timer.stop();                                                       \
  serial_timer.start();

/**
 * @brief Add RHD simulation specific command line options to the command line
 * parser.
 *
 * The following parameters are added:
 *  - output-time-unit (no abbreviation, optional, string argument): unit used
 *    for time values that are written to the Log.
 *
 * @param parser CommandLineParser that has not yet parsed the command line
 * options.
 */
void RadiationHydrodynamicsSimulation::add_command_line_parameters(
    CommandLineParser &parser) {
  parser.add_option("output-time-unit", 0,
                    "Unit for time stat output to the terminal.",
                    COMMANDLINEOPTION_STRINGARGUMENT, "s");
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
 *  - use mask: Use a mask to keep the hydrodynamics fixed in part of the box?
 *    (default: false)
 *  - use potential: Use an external point mass potential as extra gravitational
 *    force (default: false)
 *  - maximum neutral fraction: Maximum value of the hydrogen neutral fraction
 *    that is allowed at the start of a radiation step (negative values do not
 *    impose an upper limit, default: -1)
 *  - use self gravity: Add the gravitational force due to the gas itself
 *    (default: false)
 *  - use cooling: Add external cooling (default: false)
 *
 * @param parser CommandLineParser that contains the parsed command line
 * arguments.
 * @param write_output Flag indicating whether this process writes output.
 * @param programtimer Total program timer.
 * @param log Log to write logging info to.
 * @return Exit code: 0 on success.
 */
int RadiationHydrodynamicsSimulation::do_simulation(CommandLineParser &parser,
                                                    bool write_output,
                                                    Timer &programtimer,
                                                    Log *log) {

  Timer total_timer;
  Timer serial_timer;
  Timer parallel_timer;

  total_timer.start();
  serial_timer.start();

  bool every_iteration_output =
      parser.get_value< bool >("every-iteration-output");

  // set the maximum number of openmp threads
  WorkEnvironment::set_max_num_threads(
      parser.get_value< int_fast32_t >("threads"));

  // set the unit for time stats terminal output
  std::string output_time_unit =
      parser.get_value< std::string >("output-time-unit");

  // second: initialize the parameters that are read in from static files
  // these files should be configured by CMake and put in a location that is
  // stored in a CMake configured header
  LineCoolingData line_cooling_data;

  // third: read in the parameters of the run from a parameter file. This file
  // should be read by a ParameterFileParser object that acts as a dictionary
  ParameterFile params(parser.get_value< std::string >("params"));

  // fourth: construct the density grid. This should be stored in a separate
  // DensityGrid object with geometrical and physical properties
  DensityFunction *density_function =
      DensityFunctionFactory::generate(params, log);
  DensityMask *density_mask = DensityMaskFactory::generate(params, log);
  CrossSections *cross_sections = CrossSectionsFactory::generate(params, log);
  RecombinationRates *recombination_rates =
      RecombinationRatesFactory::generate(params, log);

  // initialize the simulation box
  const SimulationBox simulation_box(params);

  HydroIntegrator *hydro_integrator =
      new HydroIntegrator(simulation_box, params);
  const double hydro_total_time = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:total time", "1. s");

  double hydro_minimum_timestep = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:minimum timestep", "-1. s");
  if (hydro_minimum_timestep < 0.) {
    hydro_minimum_timestep = 1.e-10 * hydro_total_time;
  }

  double hydro_maximum_timestep = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:maximum timestep", "-1. s");
  if (hydro_maximum_timestep < 0.) {
    hydro_maximum_timestep = 0.1 * hydro_total_time;
  }

  double hydro_snaptime = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:snapshot time", "-1. s");
  if (hydro_snaptime < 0.) {
    hydro_snaptime = 0.1 * hydro_total_time;
  }
  uint_fast32_t hydro_lastsnap = 1;

  const double hydro_radtime = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:radiation time", "-1. s");
  uint_fast32_t hydro_lastrad = 0;

  const double maximum_neutral_fraction = params.get_value< double >(
      "RadiationHydrodynamicsSimulation:maximum neutral fraction", -1.);

  DensityGrid *grid =
      DensityGridFactory::generate(simulation_box, params, true, log);

  // fifth: construct the stellar sources. These should be stored in a
  // separate StellarSources object with geometrical and physical properties.
  PhotonSourceDistribution *sourcedistribution =
      PhotonSourceDistributionFactory::generate(params, log);
  int_fast32_t random_seed = params.get_value< int_fast32_t >(
      "RadiationHydrodynamicsSimulation:random seed", 42);
  PhotonSourceSpectrum *spectrum = PhotonSourceSpectrumFactory::generate(
      "PhotonSourceSpectrum", params, log);

  if (sourcedistribution != nullptr && spectrum == nullptr) {
    cmac_error("No spectrum provided for the discrete photon sources!");
  }
  if (sourcedistribution == nullptr && spectrum != nullptr) {
    cmac_warning("Discrete photon source spectrum provided, but no discrete "
                 "photon source distributions. The given spectrum will be "
                 "ignored.");
  }

  ContinuousPhotonSource *continuoussource =
      ContinuousPhotonSourceFactory::generate(simulation_box.get_box(), params,
                                              log);
  PhotonSourceSpectrum *continuousspectrum =
      PhotonSourceSpectrumFactory::generate("ContinuousPhotonSourceSpectrum",
                                            params, log);

  if (continuoussource != nullptr && continuousspectrum == nullptr) {
    cmac_error("No spectrum provided for the continuous photon sources!");
  }
  if (continuoussource == nullptr && continuousspectrum != nullptr) {
    cmac_warning("Continuous photon source spectrum provided, but no "
                 "continuous photon source. The given spectrum will be "
                 "ignored.");
  }

  Abundances abundances(params, log);

  PhotonSource source(sourcedistribution, spectrum, continuoussource,
                      continuousspectrum, abundances, *cross_sections, params,
                      log);

  // set up output
  std::string output_folder =
      Utilities::get_absolute_path(params.get_value< std::string >(
          "RadiationHydrodynamicsSimulation:output folder", "."));
  DensityGridWriter *writer =
      DensityGridWriterFactory::generate(output_folder, params, true, log);

  uint_fast32_t nloop = params.get_value< uint_fast32_t >(
      "RadiationHydrodynamicsSimulation:number of iterations", 10);

  uint_fast64_t numphoton = params.get_value< uint_fast64_t >(
      "RadiationHydrodynamicsSimulation:number of photons", 1e5);
  uint_fast64_t numphoton1 = params.get_value< uint_fast64_t >(
      "RadiationHydrodynamicsSimulation:number of photons first loop",
      numphoton);
  double Q = source.get_total_luminosity();

  ChargeTransferRates charge_transfer_rates;

  // used to calculate both the ionization state and the temperature
  TemperatureCalculator *temperature_calculator = new TemperatureCalculator(
      Q, abundances, line_cooling_data, *recombination_rates,
      charge_transfer_rates, params, log);

  // optional mask to fix the hydrodynamics in some parts of the box
  HydroMask *mask = nullptr;
  const bool use_mask = params.get_value< bool >(
      "RadiationHydrodynamicsSimulation:use mask", false);
  if (use_mask) {
    mask = HydroMaskFactory::generate(params, log);
  }

  // optional external point mass potential
  ExternalPotential *potential =
      ExternalPotentialFactory::generate(params, log);
  if (params.get_value< bool >("RadiationHydrodynamicsSimulation:use potential",
                               false) &&
      potential == nullptr) {
    cmac_error("No external potential provided!");
  }

  const bool do_self_gravity = params.get_value< bool >(
      "RadiationHydrodynamicsSimulation:use self gravity", false);

  DeRijckeRadiativeCooling *cooling_function = nullptr;
  if (params.get_value< bool >("RadiationHydrodynamicsSimulation:use cooling",
                               false)) {
    cooling_function = new DeRijckeRadiativeCooling();
  }

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output) {
    std::ofstream pfile(output_folder + "/parameters-usedvalues.param");
    params.print_contents(pfile);
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

  // done writing file, now initialize grid
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid->get_number_of_cells());
  grid->initialize(block, *density_function);

  // object used to distribute jobs in a shared memory parallel context
  WorkDistributor< IonizationPhotonShootJobMarket, IonizationPhotonShootJob >
      workdistributor(parser.get_value< int_fast32_t >("threads"));
  const int_fast32_t worksize = workdistributor.get_worksize();
  Timer worktimer;

  if (density_mask != nullptr) {
    log->write_status("Initializing DensityMask...");
    density_mask->initialize(worksize);
    log->write_status("Done initializing mask. Applying mask...");
    density_mask->apply(*grid);
    log->write_status("Done applying mask.");
  }

  TreeSelfGravity *self_gravity = nullptr;
  if (do_self_gravity) {
    self_gravity = new TreeSelfGravity(*grid, 0.25);
  }

  if (log) {
    log->write_status("Program will use ",
                      workdistributor.get_worksize_string(),
                      " for photon shooting.");
  }
  IonizationPhotonShootJobMarket photonshootjobs(source, random_seed, *grid, 0,
                                                 100, worksize);

  // initialize the hydro variables (before we write the initial snapshot)
  hydro_integrator->initialize_hydro_variables(*grid);

  // apply the mask if applicable
  if (use_mask) {
    mask->initialize_mask(*grid);
    mask->apply_mask(*grid, 0, 0);
  }

  if (write_output) {
    writer->write(*grid, 0, params, 0., hydro_integrator->get_internal_units());
  }

  double maximum_timestep = hydro_maximum_timestep;
  if (hydro_radtime > 0.) {
    // make sure the system is evolved hydrodynamically in between successive
    // ionization steps
    maximum_timestep = std::min(maximum_timestep, hydro_radtime);
  }
  TimeLine timeline(0., hydro_total_time, hydro_minimum_timestep,
                    maximum_timestep, log);
  bool has_next_step = true;
  uint_fast32_t num_step = 0;
  while (has_next_step) {
    double requested_timestep = hydro_integrator->get_maximal_timestep(*grid);
    double actual_timestep, current_time;
    has_next_step =
        timeline.advance(requested_timestep, actual_timestep, current_time);
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

    // finally: the actual program loop whereby the density grid is ray traced
    // using photon packets generated by the stellar sources
    uint_fast32_t loop = 0;
    uint_fast32_t nloop_step = nloop;

    // decide whether or not to do the radiation step
    if (hydro_radtime < 0. ||
        (current_time - actual_timestep) >= hydro_lastrad * hydro_radtime) {

      ++hydro_lastrad;
      // update the PhotonSource
      if (sourcedistribution->update(current_time)) {
        source.update(sourcedistribution);
      }

      // reset the neutral fractions if necessary
      if (maximum_neutral_fraction > 0.) {
        for (auto it = grid->begin(); it != grid->end(); ++it) {
          if (it.get_ionization_variables().get_ionic_fraction(ION_H_n) >
              maximum_neutral_fraction) {
            it.get_ionization_variables().set_ionic_fraction(
                ION_H_n, maximum_neutral_fraction);
          }
        }
      }
    } else {
      nloop_step = 0;
    }

    if (log && nloop_step > 0) {
      log->write_status("Starting radiation step...");
    }

    while (loop < nloop_step) {

      if (log) {
        log->write_status("Loop ", loop, " of ", nloop_step, ".");
      }

      uint_fast64_t lnumphoton = numphoton;

      if (loop == 0) {
        // overwrite the number of photons for the first loop (might be useful
        // if more than 1 boundary is periodic, since the initial neutral
        // fractions are very low)
        lnumphoton = numphoton1;
      }

      grid->reset_grid(*density_function);
      DiffuseReemissionHandler::set_reemission_probabilities(*grid);

      double typecount[PHOTONTYPE_NUMBER] = {0};

      double totweight = 0.;

      uint_fast64_t local_numphoton = lnumphoton;

      photonshootjobs.set_numphoton(local_numphoton);
      worktimer.start();
      start_parallel_timing_block();
      workdistributor.do_in_parallel(photonshootjobs);
      stop_parallel_timing_block();
      worktimer.stop();

      photonshootjobs.update_counters(totweight, typecount);

      start_parallel_timing_block();
      temperature_calculator->calculate_temperature(loop, totweight, *grid,
                                                    block);
      stop_parallel_timing_block();

      ++loop;

      if (write_output && every_iteration_output && loop < nloop_step) {
        writer->write(*grid, loop, params, current_time,
                      hydro_integrator->get_internal_units());
      }
    }

    if (log) {
      log->write_status("Done with radiation step.");
    }

    // update the gravitational accelerations if applicable
    if (potential != nullptr) {
      for (auto it = grid->begin(); it != grid->end(); ++it) {
        const CoordinateVector<> a =
            potential->get_acceleration(it.get_cell_midpoint());
        it.get_hydro_variables().set_gravitational_acceleration(a);
      }
    }

    if (self_gravity != nullptr) {
      self_gravity->compute_accelerations(*grid);
    }

    // update the cooling terms if applicable
    if (cooling_function != nullptr) {
      for (auto it = grid->begin(); it != grid->end(); ++it) {
        const double nH = it.get_ionization_variables().get_number_density();
        const double nH2 = nH * nH;
        const double cooling_rate =
            cooling_function->get_cooling_rate(
                it.get_ionization_variables().get_temperature()) *
            nH2 * it.get_volume();
        it.get_hydro_variables().set_energy_term(-cooling_rate);
      }
    }

    hydro_integrator->do_hydro_step(*grid, actual_timestep, serial_timer,
                                    parallel_timer);

    // apply the mask if applicable
    if (use_mask) {
      mask->apply_mask(*grid, actual_timestep, current_time);
    }

    // write snapshot
    // we don't write if this is the last snapshot, because then it is written
    // outside the integration loop
    if (write_output && hydro_lastsnap * hydro_snaptime <= current_time &&
        has_next_step) {
      writer->write(*grid, hydro_lastsnap, params, current_time,
                    hydro_integrator->get_internal_units());
      ++hydro_lastsnap;
    }
  }

  // write snapshot
  if (write_output) {
    writer->write(*grid, hydro_lastsnap, params, hydro_total_time,
                  hydro_integrator->get_internal_units());
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
    log->write_status("Total photon shooting time: ",
                      Utilities::human_readable_time(worktimer.value()), ".");
  }

  if (cooling_function != nullptr) {
    delete cooling_function;
  }
  if (self_gravity != nullptr) {
    delete self_gravity;
  }
  if (potential != nullptr) {
    delete potential;
  }
  if (use_mask) {
    delete mask;
  }
  if (sourcedistribution != nullptr) {
    delete sourcedistribution;
  }
  if (continuoussource != nullptr) {
    delete continuoussource;
  }
  delete density_function;
  if (density_mask != nullptr) {
    delete density_mask;
  }
  delete writer;
  delete grid;
  delete hydro_integrator;
  delete temperature_calculator;
  delete continuousspectrum;
  delete spectrum;

  delete cross_sections;
  delete recombination_rates;

  // we cannot delete the log, since it is still used in the destructor of
  // objects that are destructed at the return of the main program
  // this is not really a problem, as the memory is freed up by the OS anyway
  // delete log;

  return 0;
}
