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
#include "Box.hpp"
#include "CommandLineOption.hpp"
#include "CommandLineParser.hpp"
#include "CoordinateVector.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGrid.hpp"
#include "DensityGridWriter.hpp"
#include "FileLog.hpp"
#include "LineCoolingData.hpp"
#include "ParameterFile.hpp"
#include "PhotonSource.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PlanckPhotonSourceSpectrum.hpp"
#include "TerminalLog.hpp"
#include "VernerCrossSections.hpp"
#include "VernerRecombinationRates.hpp"
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
  DensityGrid grid(params, *density_function, recombination_rates, log);

  // fifth: construct the stellar sources. These should be stored in a
  // separate StellarSources object with geometrical and physical properties.
  PhotonSourceDistribution *sourcedistribution =
      PhotonSourceDistributionFactory::generate(params, log);
  PlanckPhotonSourceSpectrum spectrum;
  PhotonSource source(*sourcedistribution, spectrum, cross_sections, log);

  // finally: the actual program loop whereby the density grid is ray traced
  // using photon packets generated by the stellar sources
  // it would be nice to have some interactor classes that can do this

  // this should be an iteration
  DensityGridWriter writer("snapshot", grid);
  for (unsigned int loop = 0; loop < 10; ++loop) {
    grid.reset_grid();
    unsigned int numphoton = 1000000;
    source.set_number_of_photons(numphoton);
    log->write_status("Start shooting photons...");
    unsigned int typecount[PHOTONTYPE_NUMBER] = {0};
    for (unsigned int i = 0; i < numphoton; ++i) {
      if (!(i % 100000)) {
        log->write_status("Photon ", i, " of ", numphoton, ".");
      }
      Photon photon = source.get_random_photon();
      double tau = -std::log(Utilities::random_double());
      while (grid.interact(photon, tau)) {
        if (!source.reemit(photon, grid.get_cell_values(grid.get_cell_indices(
                                       photon.get_position())),
                           params.get_value< double >("helium_abundance"))) {
          break;
        }
        tau = -std::log(Utilities::random_double());
      }
      ++typecount[photon.get_type()];
    }
    log->write_status("Done shooting photons.");
    log->write_status(typecount[PHOTONTYPE_ABSORBED],
                      " photons were reemitted as non-ionizing photons.");
    log->write_status(typecount[PHOTONTYPE_DIFFUSE_HI] +
                          typecount[PHOTONTYPE_DIFFUSE_HeI],
                      " photons were scattered.");
    double escape_fraction =
        (100. * (numphoton - typecount[PHOTONTYPE_ABSORBED])) / numphoton;
    log->write_status("Escape fraction: ", escape_fraction, "%.");
    double escape_fraction_HI =
        (100. * typecount[PHOTONTYPE_DIFFUSE_HI]) / numphoton;
    log->write_status("Diffuse HI escape fraction: ", escape_fraction_HI, "%.");
    double escape_fraction_HeI =
        (100. * typecount[PHOTONTYPE_DIFFUSE_HeI]) / numphoton;
    log->write_status("Diffuse HeI escape fraction: ", escape_fraction_HeI,
                      "%.");
    grid.calculate_ionization_state(numphoton);

    // write snapshot
    writer.write(loop);
  }

  ofstream pfile("parameters-usedvalues.param");
  params.print_contents(pfile);

  // idea: the photons themselves should be a class, storing a current
  // position/cell, current direction, and the photon package energy
  // photons are then emitted by the StellarSources object
  // there is another object, a PhotonInteractor, that selects a random
  // optical depth for a photon and then moves it to a new position until
  // it leaves the system

  delete sourcedistribution;
  delete density_function;

  // we cannot delete the log, since it is still used in the destructor of
  // objects that are destructed at the return of the main program
  // this is not really a problem, as the memory is freed up by the OS anyway
  // delete log;

  return 0;
}
