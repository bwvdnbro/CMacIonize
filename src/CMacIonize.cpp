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
#include "CommandLineOption.hpp"
#include "CommandLineParser.hpp"
#include "CompilerInfo.hpp"
#include "ConfigurationInfo.hpp"
#include "DustSimulation.hpp"
#include "FileLog.hpp"
#include "IonizationSimulation.hpp"
#include "MPICommunicator.hpp"
#include "RadiationHydrodynamicsSimulation.hpp"
#include "TaskBasedIonizationSimulation.hpp"
#include "TaskBasedRadiationHydrodynamicsSimulation.hpp"
#include "TerminalLog.hpp"
#include "Timer.hpp"

#include <string>

/**
 * @brief Entrance point of the program.
 *
 * We first initialize the MPI environment (if active) and start the total
 * program timer. We then parse the command line arguments and perform some
 * basic checks and initialization steps. Once this is done, we branch into one
 * of three possible modes:
 *  - dusty-radiative-transfer mode: enabled by setting the command line option
 *    '--dusty-radiative-transfer'. In this mode, we set up a simple spiral
 *    galaxy model and produce an image with dust exctinction. Used for a first
 *    year lab project at the University of St Andrews.
 *  - rhd mode: enabled by setting the command line option '--rhd'. In this mode
 *    we perform a full radiation hydrodynamics (RHD) simulation.
 *  - default mode: used if no other mode is chosen. In this mode, the
 *    photoionization code is run to post-process an existing density field.
 *
 * The program accepts the following general command line arguments (more are
 * added in RadiationHydrodynamicsSimulation::add_command_line_parameters):
 *  - params ('p', required, string argument): name of the parameter file that
 *    contains all parameters for the run.
 *  - verbose ('v', optional, no argument): set the log level to the lowest
 *    possible value, which means a maximal ammount of information written to
 *    the Log.
 *  - logfile ('l', optional, string argument): if present, write logging info
 *    to a file with the given name, rather than to the terminal window
 *    (default file name: CMacIonize_run.log).
 *  - dirty ('d', optional, no argument): allow the code to run even when
 *    uncommitted changes are detected by git describe (meaning that the run is
 *    potentially irreproducible).
 *  - thread ('t', optional, integer argument): number of shared memory parallel
 *    threads to use while running the code (default: 1). If the given number is
 *    higher than the available number of threads, the maximum available number
 *    is used instead.
 *  - dry-run ('n', optional, no argument): exit the simulation before the main
 *    simulation loop starts, but after all parameters have been parsed and all
 *    objects have been created. Useful to test the parameter file and memory
 *    requirements without actually running the simulation.
 *  - use-version ('u', optional, string argument): enforce the given git
 *    version tag (no default value). If the version tag returned by git
 *    describe does not match this value, the simulation will abort.
 *  - dusty-radiative-transfer (no abbreviation, optional, no argument): run the
 *    code in dusty radiative transfer mode rather than photoionization mode.
 *  - rhd (no abbreviation, optional, no argument): run the code in RHD mode
 *    rather than photoionization mode.
 *  - output-statistics ('s', optional, no argument): output statistics about
 *    the photon packets at the end of each iteration.
 *  - task-based (no abbreviation, optional, no argument): run the code using
 *    a task-based parallel algorithm.
 *  - task-based-rhd (no abbreviation, optional, no argument): run the RHD code
 *    using a task-based parallel algorithm.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // initialize the MPI communicator and make sure only process 0 writes to the
  // log and output files
  MPICommunicator comm(argc, argv);
  bool write_log = (comm.get_rank() == 0);
  bool write_output = (comm.get_rank() == 0);

  // install signal handlers
  OperatingSystem::install_signal_handlers();

  Timer programtimer;

  // first thing we should do: parse the command line arguments
  // we need to define a CommandLineParser object that does this and acts as a
  // dictionary that can be queried
  CommandLineParser parser("CMacIonize");
  parser.add_required_option< std::string >(
      "params", 'p',
      "Name of the parameter file containing the simulation parameters.");
  parser.add_option("verbose", 'v',
                    "Set the logging level to the lowest "
                    "possible value to allow more output to be "
                    "written to the log.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("logfile", 'l',
                    "Output program logs to a file with the "
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
  parser.add_option("dusty-radiative-transfer", 0,
                    "Run a dusty radiative transfer simulation instead of an "
                    "ionization simulation.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("rhd", 0,
                    "Run a radiation hydrodynamics simulation "
                    "instead of an ionization simulation.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("output-statistics", 's',
                    "Output statistical information about the photons.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("task-based", 0,
                    "Run a task-based photoionization simulation.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("task-based-rhd", 0, "Run a task-based RHD simulation.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");

  // add simulation type specific parameters
  RadiationHydrodynamicsSimulation::add_command_line_parameters(parser);

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

  if (parser.get_value< bool >("dusty-radiative-transfer")) {

    if (comm.get_size() > 1) {
      cmac_error("MPI parallel dusty radiative transfer is not supported!");
    }
    return DustSimulation::do_simulation(parser, write_output, programtimer,
                                         log);
  } else if (parser.get_value< bool >("rhd")) {

    if (comm.get_size() > 1) {
      cmac_error("MPI RHD is not (yet) supported!");
    }
    return RadiationHydrodynamicsSimulation::do_simulation(parser, write_output,
                                                           programtimer, log);
  } else if (parser.get_value< bool >("task-based")) {

    if (comm.get_size() > 1) {
      cmac_error("MPI task based algorithm does not exist yet.");
    }

    TaskBasedIonizationSimulation simulation(
        parser.get_value< int_fast32_t >("threads"),
        parser.get_value< std::string >("params"), log);

    if (parser.get_value< bool >("dry-run")) {
      if (log) {
        log->write_warning("Dry run requested. Program will now halt.");
      }
      return 0;
    }

    simulation.initialize();
    simulation.run();

    programtimer.stop();

    size_t memory_usage = OperatingSystem::get_peak_memory_usage();
    comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(memory_usage);
    if (log) {
      log->write_status("Total program time: ",
                        Utilities::human_readable_time(programtimer.value()),
                        ".");
      log->write_status("Peak memory usage: ",
                        Utilities::human_readable_bytes(memory_usage), ".");
    }
    return 0;

  } else if (parser.get_value< bool >("task-based-rhd")) {

    if (comm.get_size() > 1) {
      cmac_error("MPI RHD is not (yet) supported!");
    }
    return TaskBasedRadiationHydrodynamicsSimulation::do_simulation(
        parser, write_output, programtimer, log);
  } else {

    IonizationSimulation simulation(
        write_output, parser.get_value< bool >("every-iteration-output"),
        parser.get_value< bool >("output-statistics"),
        parser.get_value< int_fast32_t >("threads"),
        parser.get_value< std::string >("params"), &comm, log);

    if (parser.get_value< bool >("dry-run")) {
      if (log) {
        log->write_warning("Dry run requested. Program will now halt.");
      }
      return 0;
    }

    simulation.initialize();
    simulation.run();

    programtimer.stop();

    size_t memory_usage = OperatingSystem::get_peak_memory_usage();
    comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(memory_usage);
    if (log) {
      log->write_status("Total program time: ",
                        Utilities::human_readable_time(programtimer.value()),
                        ".");
      log->write_status("Peak memory usage: ",
                        Utilities::human_readable_bytes(memory_usage), ".");
    }
    return 0;
  }
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
