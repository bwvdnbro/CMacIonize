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
#include "DensityFunctionFactory.hpp"
#include "DensityGridFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DensityMaskFactory.hpp"
#include "DiffuseReemissionHandler.hpp"
#include "DustSimulation.hpp"
#include "EmissivityCalculator.hpp"
#include "FileLog.hpp"
#include "HydroIntegrator.hpp"
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
#include "Timer.hpp"
#include "WorkDistributor.hpp"
#include "WorkEnvironment.hpp"

#include <iostream>
#include <string>

/**
 * @brief Perform an RHD simulation.
 *
 * This method reads the following parameters from the parameter file:
 *  - timestep: Hydrodynamical time step (default: 0.01 s)
 *  - total time: Total simulation time (default: 1. s)
 *  - snapshot time: Time interval between consecutive snapshot dumps (default:
 *    0.1*(total time))
 *  - random seed: Seed for the random number generator (default: 42)
 *  - output folder: Folder where all output files will be placed (default: .)
 *  - number of iterations: Number of iterations of the photoionization
 *    algorithm (default: 10)
 *  - number of photons: Number of photons to use during each iteration of the
 *    photoionization algorithm (default: 1e5)
 *  - number of photons first loop: Number of photons to use during the first
 *    iteration of the photoionization algorithm (default: (number of photons))
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

  bool every_iteration_output =
      parser.get_value< bool >("every-iteration-output");

  // set the maximum number of openmp threads
  WorkEnvironment::set_max_num_threads(parser.get_value< int >("threads"));

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
  const double hydro_timestep = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:timestep", "0.01 s");
  const double hydro_total_time = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:total time", "1. s");
  const unsigned int numstep = hydro_total_time / hydro_timestep;
  double hydro_snaptime = params.get_physical_value< QUANTITY_TIME >(
      "RadiationHydrodynamicsSimulation:snapshot time", "-1. s");
  if (hydro_snaptime < 0.) {
    hydro_snaptime = 0.1 * hydro_total_time;
  }
  unsigned int hydro_lastsnap = 1;

  DensityGrid *grid =
      DensityGridFactory::generate(simulation_box, params, true, log);

  // fifth: construct the stellar sources. These should be stored in a
  // separate StellarSources object with geometrical and physical properties.
  PhotonSourceDistribution *sourcedistribution =
      PhotonSourceDistributionFactory::generate(params, log);
  int random_seed = params.get_value< int >(
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
      DensityGridWriterFactory::generate(output_folder, params, log);

  unsigned int nloop = params.get_value< unsigned int >(
      "RadiationHydrodynamicsSimulation:number of iterations", 10);

  unsigned int numphoton = params.get_value< unsigned int >(
      "RadiationHydrodynamicsSimulation:number of photons", 1e5);
  unsigned int numphoton1 = params.get_value< unsigned int >(
      "RadiationHydrodynamicsSimulation:number of photons first loop",
      numphoton);
  double Q = source.get_total_luminosity();

  ChargeTransferRates charge_transfer_rates;

  // used to calculate both the ionization state and the temperature
  TemperatureCalculator *temperature_calculator = new TemperatureCalculator(
      Q, abundances, line_cooling_data, *recombination_rates,
      charge_transfer_rates, params, log);

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
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid->get_number_of_cells());
  grid->initialize(block, *density_function);

  // grid->initialize initialized:
  // - densities
  // - temperatures
  // - ionic fractions
  // we have to gather these across all processes

  // this is currently BROKEN...
  //  comm.gather(grid->get_number_density_handle());
  //  comm.gather(grid->get_temperature_handle());
  //  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
  //    IonName ion = static_cast< IonName >(i);
  //    comm.gather(grid->get_ionic_fraction_handle(ion));
  //  }

  // object used to distribute jobs in a shared memory parallel context
  WorkDistributor< IonizationPhotonShootJobMarket, IonizationPhotonShootJob >
  workdistributor(parser.get_value< int >("threads"));
  const int worksize = workdistributor.get_worksize();
  Timer worktimer;

  if (density_mask != nullptr) {
    log->write_status("Initializing DensityMask...");
    density_mask->initialize(worksize);
    log->write_status("Done initializing mask. Applying mask...");
    density_mask->apply(*grid);
    log->write_status("Done applying mask.");
  }

  if (log) {
    log->write_status("Program will use ",
                      workdistributor.get_worksize_string(),
                      " for photon shooting.");
  }
  IonizationPhotonShootJobMarket photonshootjobs(source, random_seed, *grid, 0,
                                                 100, worksize);

  if (hydro_integrator != nullptr) {
    // initialize the hydro variables (before we write the initial snapshot)
    hydro_integrator->initialize_hydro_variables(*grid);
  }

  if (write_output) {
    writer->write(*grid, 0, params);
  }

  for (unsigned int istep = 0; istep < numstep; ++istep) {
    if (log) {
      log->write_status("Starting hydro step ", istep, ".");
    }

    // finally: the actual program loop whereby the density grid is ray traced
    // using photon packets generated by the stellar sources
    unsigned int loop = 0;
    while (loop < nloop) {

      if (log) {
        log->write_status("Starting loop ", loop, ".");
      }

      //    if (loop == 3 || loop == 9) {
      //      numphoton *= 10;
      //    }

      unsigned int lnumphoton = numphoton;

      if (loop == 0) {
        // overwrite the number of photons for the first loop (might be useful
        // if more than 1 boundary is periodic, since the initial neutral
        // fractions are very low)
        lnumphoton = numphoton1;
      }

      grid->reset_grid(*density_function);
      DiffuseReemissionHandler::set_reemission_probabilities(*grid);
      if (log) {
        log->write_status("Start shooting ", lnumphoton, " photons...");
      }

      double typecount[PHOTONTYPE_NUMBER] = {0};

      double totweight = 0.;

      unsigned int local_numphoton = lnumphoton;

      photonshootjobs.set_numphoton(local_numphoton);
      worktimer.start();
      workdistributor.do_in_parallel(photonshootjobs);
      worktimer.stop();

      photonshootjobs.update_counters(totweight, typecount);

      if (log) {
        log->write_status("Done shooting photons.");
        log->write_status(
            100. * typecount[PHOTONTYPE_ABSORBED] / totweight,
            "% of photons were reemitted as non-ionizing photons.");
        log->write_status(100. * (typecount[PHOTONTYPE_DIFFUSE_HI] +
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
        log->write_status("Escape fraction: ", escape_fraction, "%.");
        double escape_fraction_HI =
            (100. * typecount[PHOTONTYPE_DIFFUSE_HI]) / totweight;
        log->write_status("Diffuse HI escape fraction: ", escape_fraction_HI,
                          "%.");
        double escape_fraction_HeI =
            (100. * typecount[PHOTONTYPE_DIFFUSE_HeI]) / totweight;
        log->write_status("Diffuse HeI escape fraction: ", escape_fraction_HeI,
                          "%.");
      }

      if (log) {
        log->write_status("Calculating ionization state after shooting ",
                          lnumphoton, " photons...");
      }

      // reduce the mean intensity integrals and heating terms across all
      // processes

      // this code is currently BROKEN...
      //      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      //        IonName ion = static_cast< IonName >(i);
      //        comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(
      //            grid->get_mean_intensity_handle(ion));
      //      }
      //      if (calculate_temperature && loop > 3) {
      //        comm.reduce< MPI_SUM_OF_ALL_PROCESSES
      //        >(grid->get_heating_H_handle());
      //        comm.reduce< MPI_SUM_OF_ALL_PROCESSES
      //        >(grid->get_heating_He_handle());
      //      }

      temperature_calculator->calculate_temperature(loop, totweight, *grid,
                                                    block);

      // the calculation above will have changed the ionic fractions, and might
      // have changed the temperatures
      // we have to gather these across all processes

      // this is currently BROKEN...
      //      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      //        IonName ion = static_cast< IonName >(i);
      //        comm.gather(grid->get_ionic_fraction_handle(ion));
      //      }
      //      if (calculate_temperature && loop > 3) {
      //        comm.gather(grid->get_temperature_handle());
      //      }

      if (log) {
        log->write_status("Done calculating ionization state.");
      }

      // calculate emissivities
      // we disabled this, since we now have the post-processing Python library
      // for this
      //    if (loop > 3 && abundances.get_abundance(ELEMENT_He) > 0.) {
      //      emissivity_calculator.calculate_emissivities(*grid);
      //    }

      ++loop;

      if (write_output && every_iteration_output && loop < nloop) {
        writer->write(*grid, loop, params);
      }
    }

    if (log && loop == nloop) {
      log->write_status("Maximum number of iterations (", nloop,
                        ") reached, stopping.");
    }

    if (hydro_integrator != nullptr) {
      hydro_integrator->do_hydro_step(*grid, hydro_timestep);

      // write snapshot
      if (write_output &&
          hydro_lastsnap * hydro_snaptime < (istep + 1) * hydro_timestep) {
        writer->write(*grid, hydro_lastsnap, params,
                      hydro_lastsnap * hydro_snaptime);
        ++hydro_lastsnap;
      }
    }
  }

  // write snapshot
  if (write_output) {
    if (hydro_integrator == nullptr) {
      writer->write(*grid, nloop, params);
    } else {
      writer->write(*grid, hydro_lastsnap, params,
                    hydro_lastsnap * hydro_snaptime);
    }
  }

  programtimer.stop();
  if (log) {
    log->write_status("Total program time: ",
                      Utilities::human_readable_time(programtimer.value()),
                      ".");
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
  if (density_mask != nullptr) {
    delete density_mask;
  }
  delete writer;
  delete grid;
  if (hydro_integrator != nullptr) {
    delete hydro_integrator;
  }
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
