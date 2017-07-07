/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file CMacIonize.cpp
 *
 * @brief Entrance point of the CMacIonize program
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
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
#include "DensityFunctionFactory.hpp"
#include "DensityGridFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DensityMaskFactory.hpp"
#include "EmissivityCalculator.hpp"
#include "FileLog.hpp"
#include "HydroIntegrator.hpp"
#include "IonizationStateCalculator.hpp"
#include "LineCoolingData.hpp"
#include "MPICommunicator.hpp"
#include "ParameterFile.hpp"
#include "PhotonShootJobMarket.hpp"
#include "PhotonSource.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "TemperatureCalculator.hpp"
#include "TerminalLog.hpp"
#include "Timer.hpp"
#include "VernerCrossSections.hpp"
#include "VernerRecombinationRates.hpp"
#include "WorkDistributor.hpp"
#include "WorkEnvironment.hpp"

#include <iostream>
#include <string>

using namespace std;

/**
 * @brief Entrance point of the program
 *
 * @param argc Number of command line arguments passed on to the program.
 * @param argv Array containing the command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // initialize the MPI communicator and make sure only process 0 writes to the
  // log and output files
  MPICommunicator comm(argc, argv);
  bool write_log = (comm.get_rank() == 0);
  bool write_output = (comm.get_rank() == 0);

  Timer programtimer;

  // first thing we should do: parse the command line arguments
  // we need to define a CommandLineParser object that does this and acts as a
  // dictionary that can be queried
  CommandLineParser parser("CMacIonize");
  parser.add_required_option< string >(
      "params", 'p',
      "Name of the parameter file containing the simulation parameters.");
  parser.add_option("verbose", 'v', "Set the logging level to the lowest "
                                    "possible value to allow more output to be "
                                    "written to the log.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("logfile", 'l', "Output program logs to a file with the "
                                    "given name, instead of to the standard "
                                    "output.",
                    COMMANDLINEOPTION_STRINGARGUMENT, "CMacIonize_run.log");
  parser.add_option("dirty", 'd',
                    "Allow running a dirty code version. This is disabled by "
                    "default, since a dirty code version does not correspond "
                    "to a unique revision number in the code repository, and "
                    "it is therefore impossible to rerun a dirty version with "
                    "exactly the same code afterwards.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("threads", 't', "Number of parallel threads to use.",
                    COMMANDLINEOPTION_INTARGUMENT, "1");
  parser.add_option("dry-run", 'n',
                    "Perform a dry run of the program: this reads the "
                    "parameter file and sets up all the components, but aborts "
                    "before initializing the density grid. This option is "
                    "ideal for checking if a parameter file will work, and to "
                    "check if all input files can be read.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("every-iteration-output", 'e',
                    "Output a snapshot for every iteration of the photon "
                    "traversal algorithm.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("use-version", 'u',
                    "Make sure the code version string matches the given "
                    "string. This is useful when running simulations remotely "
                    "using a workflow system, to ensure that every remote node "
                    "is running the same code version.",
                    COMMANDLINEOPTION_STRINGARGUMENT);
  parser.parse_arguments(argc, argv);

  LogLevel loglevel = LOGLEVEL_STATUS;
  if (parser.get_value< bool >("verbose")) {
    loglevel = LOGLEVEL_INFO;
  }
  Log *log = nullptr;
  if (write_log) {
    if (parser.was_found("logfile")) {
      log = new FileLog(parser.get_value< std::string >("logfile"), loglevel);
    } else {
      // ASCII art generated using
      // http://patorjk.com/software/taag/#p=display&h=2&f=Big&t=CMacIonize
      // all '\' have been manualy escaped, so the actual result looks a bit
      // nicer
      std::string header =
          "  _____ __  __            _____            _\n"
          " / ____|  \\/  |          |_   _|          (_)\n"
          "| |    | \\  / | __ _  ___  | |  ___  _ __  _ _______\n"
          "| |    | |\\/| |/ _` |/ __| | | / _ \\| '_ \\| |_  / _ \\\n"
          "| |____| |  | | (_| | (__ _| || (_) | | | | |/ /  __/\n"
          " \\_____|_|  |_|\\__,_|\\___|_____\\___/|_| |_|_/___\\___|\n";
      log = new TerminalLog(loglevel, header);
    }
  }

  if (parser.get_value< std::string >("use-version") != "") {
    if (CompilerInfo::get_git_version() !=
        parser.get_value< std::string >("use-version")) {
      if (log) {
        log->write_error("Wrong code version (requested ",
                         parser.get_value< std::string >("use-version"),
                         ", got ", CompilerInfo::get_git_version(), ")!");
      }
      cmac_error("Wrong code version (%s != %s)!",
                 parser.get_value< std::string >("use-version").c_str(),
                 CompilerInfo::get_git_version().c_str());
    }
  }

  if (log) {
    log->write_status("This is CMacIonize, version ",
                      CompilerInfo::get_git_version(), ".");
    log->write_status("Code was compiled on ", CompilerInfo::get_full_date(),
                      " using ", CompilerInfo::get_full_compiler_name(), ".");
    log->write_status("Code was compiled for ", CompilerInfo::get_os_name(),
                      ", ", CompilerInfo::get_kernel_name(), " on ",
                      CompilerInfo::get_hardware_name(), " (",
                      CompilerInfo::get_host_name(), ").");
    log->write_status("Configuration options:");
    for (auto it = ConfigurationInfo::begin(); it != ConfigurationInfo::end();
         ++it) {
      log->write_status(it.get_key(), ": ", it.get_value());
    }
  }

  if (log) {
    if (comm.get_size() > 1) {
      log->write_status("Code is running on ", comm.get_size(), " processes.");
    } else {
      log->write_status("Code is running on a single process.");
    }
  }

  if (CompilerInfo::is_dirty()) {
    if (log) {
      log->write_warning(
          "This is a dirty code version (meaning some of the "
          "source files have changed since the code was obtained "
          "from the repository).");
    }
    if (!parser.get_value< bool >("dirty")) {
      if (log) {
        log->write_error("Running a dirty code version is disabled by default. "
                         "If you still want to run this version, add the "
                         "\"--dirty\" flag to the run command.");
      }
      cmac_error("Running a dirty code version is disabled by default.");
    } else {
      if (log) {
        log->write_warning("However, dirty running is enabled.");
      }
    }
  }

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
  ParameterFile params(parser.get_value< string >("params"));

  // fourth: construct the density grid. This should be stored in a separate
  // DensityGrid object with geometrical and physical properties
  DensityFunction *density_function =
      DensityFunctionFactory::generate(params, log);
  DensityMask *density_mask = DensityMaskFactory::generate(params, log);
  VernerCrossSections cross_sections;
  VernerRecombinationRates recombination_rates;

  HydroIntegrator *hydro_integrator = nullptr;
  double hydro_timestep = 0.;
  unsigned int numstep = 1;
  double hydro_snaptime = 0.;
  unsigned int hydro_lastsnap = 1;
  if (params.get_value< bool >("hydro:active", false)) {
    hydro_integrator = new HydroIntegrator(params);
    hydro_timestep =
        params.get_physical_value< QUANTITY_TIME >("hydro:timestep", "0.01 s");
    double hydro_total_time =
        params.get_physical_value< QUANTITY_TIME >("hydro:total_time", "1. s");
    numstep = hydro_total_time / hydro_timestep;
    hydro_snaptime =
        params.get_physical_value< QUANTITY_TIME >("hydro:snaptime", "-1. s");
    if (hydro_snaptime < 0.) {
      hydro_snaptime = 0.1 * hydro_total_time;
    }
  }

  DensityGrid *grid =
      DensityGridFactory::generate(params, *density_function, log);

  // fifth: construct the stellar sources. These should be stored in a
  // separate StellarSources object with geometrical and physical properties.
  PhotonSourceDistribution *sourcedistribution =
      PhotonSourceDistributionFactory::generate(params, log);
  int random_seed = params.get_value< int >("random_seed", 42);
  PhotonSourceSpectrum *spectrum = PhotonSourceSpectrumFactory::generate(
      "photonsourcespectrum", params, log);

  if (sourcedistribution != nullptr && spectrum == nullptr) {
    cmac_error("No spectrum provided for the discrete photon sources!");
  }
  if (sourcedistribution == nullptr && spectrum != nullptr) {
    cmac_warning("Discrete photon source spectrum provided, but no discrete "
                 "photon source distributions. The given spectrum will be "
                 "ignored.");
  }

  ContinuousPhotonSource *continuoussource =
      ContinuousPhotonSourceFactory::generate(params, log);
  PhotonSourceSpectrum *continuousspectrum =
      PhotonSourceSpectrumFactory::generate("continuousphotonsourcespectrum",
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
                      continuousspectrum, abundances, cross_sections, log);

  // set up output
  DensityGridWriter *writer =
      DensityGridWriterFactory::generate(params, *grid, log);

  unsigned int nloop =
      params.get_value< unsigned int >("max_number_iterations", 10);

  unsigned int numphoton =
      params.get_value< unsigned int >("number of photons", 100);
  unsigned int numphoton1 =
      params.get_value< unsigned int >("number of photons init", numphoton);
  double Q = source.get_total_luminosity();

  ChargeTransferRates charge_transfer_rates;

  // used to calculate the ionization state at fixed temperature
  IonizationStateCalculator ionization_state_calculator(
      Q, abundances, recombination_rates, charge_transfer_rates);

  bool calculate_temperature =
      params.get_value< bool >("calculate_temperature", true);

  TemperatureCalculator *temperature_calculator = nullptr;
  if (calculate_temperature) {
    // used to calculate both the ionization state and the temperature
    temperature_calculator = new TemperatureCalculator(
        Q, abundances, params.get_value< double >("pahfac", 1.),
        params.get_value< double >("crfac", 0.),
        params.get_value< double >("crlim", 0.75),
        params.get_physical_value< QUANTITY_LENGTH >("crscale", "1.33333 kpc"),
        line_cooling_data, recombination_rates, charge_transfer_rates, log);
  }

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output) {
    std::string folder = Utilities::get_absolute_path(
        params.get_value< std::string >("densitygridwriter:folder", "."));
    ofstream pfile(folder + "/parameters-usedvalues.param");
    params.print_contents(pfile);
    pfile.close();
    if (log) {
      log->write_status("Wrote used parameters to ", folder,
                        "/parameters-usedvalues.param.");
    }
  }

  if (parser.get_value< bool >("dry-run")) {
    if (log) {
      log->write_warning("Dry run requested. Program will now halt.");
    }
    return 0.;
  }

  // done writing file, now initialize grid
  std::pair< unsigned long, unsigned long > block =
      comm.distribute_block(0, grid->get_number_of_cells());
  grid->initialize(block);

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
  WorkDistributor< PhotonShootJobMarket, PhotonShootJob > workdistributor(
      parser.get_value< int >("threads"));
  const int worksize = workdistributor.get_worksize();
  Timer worktimer;

  if (density_mask != nullptr) {
    log->write_status("Initializing DensityMask...");
    density_mask->initialize(worksize);
    log->write_status("Done initializing mask. Applying mask...");
    density_mask->apply(*grid);
    log->write_status("Done applying mask.");
  }

  // make sure every thread on every process has another random seed
  random_seed += comm.get_rank() * worksize;

  if (log) {
    log->write_status("Program will use ",
                      workdistributor.get_worksize_string(),
                      " for photon shooting.");
  }
  PhotonShootJobMarket photonshootjobs(source, random_seed, *grid, 0, 100,
                                       worksize);

  if (hydro_integrator != nullptr) {
    // initialize the hydro variables (before we write the initial snapshot)
    hydro_integrator->initialize_hydro_variables(*grid);
  }

  if (write_output) {
    writer->write(0, params);
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

      grid->reset_grid();
      if (log) {
        log->write_status("Start shooting ", lnumphoton, " photons...");
      }

      double typecount[PHOTONTYPE_NUMBER] = {0};

      double totweight = 0.;

      unsigned int local_numphoton = lnumphoton;

      // make sure this process does only part of the total number of photons
      local_numphoton = comm.distribute(local_numphoton);

      photonshootjobs.set_numphoton(local_numphoton);
      worktimer.start();
      workdistributor.do_in_parallel(photonshootjobs);
      worktimer.stop();

      photonshootjobs.update_counters(totweight, typecount);

      // make sure the total weight and typecount is reduced across all
      // processes
      comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(totweight);
      comm.reduce< MPI_SUM_OF_ALL_PROCESSES, PHOTONTYPE_NUMBER >(typecount);

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

      if (calculate_temperature && loop > 3) {
        temperature_calculator->calculate_temperature(totweight, *grid, block);
      } else {
        ionization_state_calculator.calculate_ionization_state(totweight, *grid,
                                                               block);
      }

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
        writer->write(loop, params);
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
        writer->write(hydro_lastsnap, params, hydro_lastsnap * hydro_snaptime);
        ++hydro_lastsnap;
      }
    }
  }

  // write snapshot
  if (write_output) {
    if (hydro_integrator == nullptr) {
      writer->write(nloop, params);
    } else {
      writer->write(hydro_lastsnap, params, hydro_lastsnap * hydro_snaptime);
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

  // we cannot delete the log, since it is still used in the destructor of
  // objects that are destructed at the return of the main program
  // this is not really a problem, as the memory is freed up by the OS anyway
  // delete log;

  return 0;
}

/**
 * @mainpage
 *
 * @author Bert Vandenbroucke \n
 *         School of Physics and Astronomy \n
 *         University of St Andrews \n
 *         North Haugh \n
 *         St Andrews \n
 *         Fife \n
 *         KY16 9SS \n
 *         Scotland \n
 *         United Kingdom \n
 *         bv7@st-andrews.ac.uk
 *
 * @section purpose Purpose of the program.
 *
 * CMacIonize is the C++ version of Kenny Wood's photoionization code
 * (Wood, Mathis & Ercolano, 2004, MNRAS, 348, 1337). The name @c CMacIonize is
 * based on the name of Kenny Wood's code (@c mcionize, which stands for Monte
 * Carlo ionization code), but also reflects the fact that it is written in
 * C(++) (not Fortran, like the old code), and the fact that development started
 * during the first week of the main author's post doc in Scotland.
 *
 * The code can be used to perform 3D radiative transfer calculations on a
 * number of different grid structures, using various possible sources of
 * ionizing radiation with various possible spectra. One of the main goals
 * during development was to provide a very general and user-friendly interface,
 * so that the code can be used in a wide variety of possible scenarios, taking
 * input data from a wide variety of possible data formats, with a minimal
 * technical involvement of the end user. This automatically led to a highly
 * modular code design, whereby most functionality is encoded into classes with
 * a limited number of responsibilities that are covered by extensive unit
 * tests. Adding a new functionality should in principle be as simple as adding
 * a new class that implements an already defined general class interface.
 *
 * The code also offers an extensive framework of utility functions and classes
 * that simplify various aspects of new code development, like for example
 * high level wrappers around HDF5 read and write functions, and a very
 * intuitive unit system.
 *
 * The code is currently only parallelized for use on shared memory system using
 * OpenMP, but there are plans to also port it to larger distributed memory
 * systems using MPI. The current OpenMP implementation shows reasonable speed
 * ups on small systems, although no formal scaling tests have been performed
 * yet.
 *
 * @section structure Structure of the program.
 *
 * Due to time constraints, there is no extensive code manual (yet). However,
 * all files, functions, classes... in the source code are fully documented
 * (the Doxygen configuration enforces this), so most of it should be easy to
 * understand. The main program entry point can be found in the file
 * CMacIonize.cpp.
 */
