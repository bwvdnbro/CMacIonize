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
#include "ContinuousPhotonSourceFactory.hpp"
#include "CoordinateVector.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "EmissivityCalculator.hpp"
#include "FaucherGiguerePhotonSourceSpectrum.hpp"
#include "FileLog.hpp"
#include "IonizationStateCalculator.hpp"
#include "IterationConvergenceCheckerFactory.hpp"
#include "LineCoolingData.hpp"
#include "ParameterFile.hpp"
#include "PhotonNumberConvergenceCheckerFactory.hpp"
#include "PhotonShootJobMarket.hpp"
#include "PhotonSource.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PlanckPhotonSourceSpectrum.hpp"
#include "TemperatureCalculator.hpp"
#include "TerminalLog.hpp"
#include "Timer.hpp"
#include "VernerCrossSections.hpp"
#include "VernerRecombinationRates.hpp"
#include "WorkDistributor.hpp"
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
  parser.add_option("dirty", 'd', "Allow running a dirty code version.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.add_option("threads", 't', "Number of parallel threads to use.",
                    COMMANDLINEOPTION_INTARGUMENT, "1");
  parser.parse_arguments(argc, argv);

  LogLevel loglevel = LOGLEVEL_STATUS;
  if (parser.get_value< bool >("verbose")) {
    loglevel = LOGLEVEL_INFO;
  }
  Log *log;
  if (parser.was_found("logfile")) {
    log = new FileLog(parser.get_value< std::string >("logfile"), loglevel);
  } else {
    log = new TerminalLog(loglevel);
  }

  log->write_status("This is CMacIonize, version ",
                    CompilerInfo::get_git_version(), ".");
  log->write_status("Code was compiled on ", CompilerInfo::get_full_date(),
                    " using ", CompilerInfo::get_full_compiler_name(), ".");
  log->write_status("Code was compiled for ", CompilerInfo::get_os_name(), ", ",
                    CompilerInfo::get_kernel_name(), " on ",
                    CompilerInfo::get_hardware_name(), " (",
                    CompilerInfo::get_host_name(), ").");

  if (CompilerInfo::is_dirty()) {
    log->write_warning("This is a dirty code version.");
    if (!parser.get_value< bool >("dirty")) {
      log->write_error("Running a dirty code version is disabled by default. "
                       "If you still want to run this version, add the "
                       "\"--dirty\" flag to the run command.");
      cmac_error("Running a dirty code version is disabled by default.");
    } else {
      log->write_warning("However, dirty running is enabled.");
    }
  }

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
  VernerCrossSections cross_sections;
  VernerRecombinationRates recombination_rates;
  DensityGrid *grid =
      DensityGridFactory::generate(params, *density_function, log);

  // fifth: construct the stellar sources. These should be stored in a
  // separate StellarSources object with geometrical and physical properties.
  PhotonSourceDistribution *sourcedistribution =
      PhotonSourceDistributionFactory::generate(params, log);
  int random_seed = params.get_value< int >("random_seed", 42);
  PlanckPhotonSourceSpectrum spectrum(params, log);

  IsotropicContinuousPhotonSource *continuoussource =
      ContinuousPhotonSourceFactory::generate(params, log);
  FaucherGiguerePhotonSourceSpectrum continuousspectrum(params, log);

  Abundances abundances(params, log);

  PhotonSource source(sourcedistribution, &spectrum, continuoussource,
                      &continuousspectrum, abundances, cross_sections, log);

  // set up output
  DensityGridWriter *writer =
      DensityGridWriterFactory::generate(params, *grid, log);

  // set up convergence checking
  PhotonNumberConvergenceChecker *convergence_checker =
      PhotonNumberConvergenceCheckerFactory::generate(*grid, params, log);

  unsigned int nloop =
      params.get_value< unsigned int >("max_number_iterations", 10);

  unsigned int numphoton =
      params.get_value< unsigned int >("number of photons", 100);
  double Q = source.get_total_luminosity();

  ChargeTransferRates charge_transfer_rates;

  // used to calculate the ionization state at fixed temperature
  IonizationStateCalculator ionization_state_calculator(
      Q, abundances, recombination_rates, charge_transfer_rates);
  // used to calculate both the ionization state and the temperature
  TemperatureCalculator temperature_calculator(
      Q, abundances, params.get_value< double >("pahfac", 1.),
      line_cooling_data, recombination_rates, charge_transfer_rates);

  bool calculate_temperature =
      params.get_value< bool >("calculate_temperature", true);

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file
  std::string folder = Utilities::get_absolute_path(
      params.get_value< std::string >("densitygridwriter.folder", "."));
  ofstream pfile(folder + "/parameters-usedvalues.param");
  params.print_contents(pfile);
  pfile.close();

  // finally: the actual program loop whereby the density grid is ray traced
  // using photon packets generated by the stellar sources

  IterationConvergenceChecker *itconvergence_checker =
      IterationConvergenceCheckerFactory::generate(*grid, params, log);

  // object used to distribute jobs in a shared memory parallel context
  WorkDistributor workdistributor(parser.get_value< int >("threads"));
  const int worksize = workdistributor.get_worksize();
  Timer worktimer;

  log->write_status("Program will use ", worksize,
                    " parallel threads for photon shooting.");

  PhotonShootJobMarket photonshootjobs(source, random_seed, *grid, 0, 100,
                                       worksize);

  writer->write(0, params);
  unsigned int loop = 0;
  while (loop < nloop && !itconvergence_checker->is_converged()) {

    log->write_status("Starting loop ", loop, ".");

    // run the number of photons by the IterationConvergenceChecker to allow for
    // corrections
    numphoton = itconvergence_checker->get_next_number_of_photons(numphoton);

    //    if (loop == 3 || loop == 9) {
    //      numphoton *= 10;
    //    }

    unsigned int lnumphoton = numphoton;
    grid->reset_grid();
    log->write_status("Start shooting photons...");
    log->write_status("Initial sub step number: ", lnumphoton, ".");

    double typecount[PHOTONTYPE_NUMBER] = {0};

    unsigned int numsubstep = 0;
    unsigned int totnumphoton = 0;
    double totweight = 0.;
    while (!convergence_checker->is_converged(totnumphoton)) {
      log->write_info("Substep ", numsubstep);

      photonshootjobs.set_numphoton(lnumphoton);
      worktimer.start();
      workdistributor.do_in_parallel(photonshootjobs);
      worktimer.stop();

      totnumphoton += lnumphoton;
      lnumphoton = convergence_checker->get_number_of_photons_next_substep(
          lnumphoton, totnumphoton);

      photonshootjobs.update_counters(totweight, typecount);

      ++numsubstep;
    }
    lnumphoton = totnumphoton;
    log->write_status("Done shooting photons.");
    log->write_status(100. * typecount[PHOTONTYPE_ABSORBED] / totweight,
                      "% of photons were reemitted as non-ionizing photons.");
    log->write_status(100. * (typecount[PHOTONTYPE_DIFFUSE_HI] +
                              typecount[PHOTONTYPE_DIFFUSE_HeI]) /
                          totweight,
                      "% of photons were scattered.");
    double escape_fraction =
        (100. * (totweight - typecount[PHOTONTYPE_ABSORBED])) / totweight;
    // since totweight is updated in chunks, while the counters are updated
    // per photon, round off might cause totweight to be slightly smaller than
    // the counter value. This gives (strange looking) negative escape
    // fractions, which we reset to 0 here.
    escape_fraction = std::max(0., escape_fraction);
    log->write_status("Escape fraction: ", escape_fraction, "%.");
    double escape_fraction_HI =
        (100. * typecount[PHOTONTYPE_DIFFUSE_HI]) / totweight;
    log->write_status("Diffuse HI escape fraction: ", escape_fraction_HI, "%.");
    double escape_fraction_HeI =
        (100. * typecount[PHOTONTYPE_DIFFUSE_HeI]) / totweight;
    log->write_status("Diffuse HeI escape fraction: ", escape_fraction_HeI,
                      "%.");

    log->write_status("Calculating ionization state after shooting ",
                      lnumphoton, " photons...");
    if (calculate_temperature && loop > 3) {
      temperature_calculator.calculate_temperature(totweight, *grid);
    } else {
      ionization_state_calculator.calculate_ionization_state(totweight, *grid);
    }
    log->write_status("Done calculating ionization state.");

    // calculate emissivities
    // we disabled this, since we now have the post-processing Python library
    // for this
    //    if (loop > 3 && abundances.get_abundance(ELEMENT_He) > 0.) {
    //      emissivity_calculator.calculate_emissivities(*grid);
    //    }

    // use the current number of photons as a guess for the new number
    numphoton = convergence_checker->get_new_number_of_photons(lnumphoton);

// print out a curve that shows the evolution of chi2
#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_CHI2_CURVE
#pragma message "Outputting chi2 curve!"
    ((ChiSquaredPhotonNumberConvergenceChecker *)convergence_checker)
        ->output_chi2_curve(loop);
#endif

    ++loop;
  }

  if (loop == nloop) {
    log->write_status("Maximum number of iterations (", nloop,
                      ") reached, stopping.");
  }

  // write snapshot
  writer->write(loop - 1, params);

  programtimer.stop();
  log->write_status("Total program time: ",
                    Utilities::human_readable_time(programtimer.value()), ".");
  log->write_status("Total photon shooting time: ",
                    Utilities::human_readable_time(worktimer.value()), ".");

  if (sourcedistribution != nullptr) {
    delete sourcedistribution;
  }
  if (continuoussource != nullptr) {
    delete continuoussource;
  }
  delete density_function;
  delete writer;
  delete grid;
  delete itconvergence_checker;
  delete convergence_checker;

  // we cannot delete the log, since it is still used in the destructor of
  // objects that are destructed at the return of the main program
  // this is not really a problem, as the memory is freed up by the OS anyway
  // delete log;

  return 0;
}
