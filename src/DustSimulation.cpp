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
 * @file DustSimulation.cpp
 *
 * @brief DustSimulation implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DustSimulation.hpp"
#include "Abundances.hpp"
#include "Box.hpp"
#include "CartesianDensityGrid.hpp"
#include "CommandLineParser.hpp"
#include "CompilerInfo.hpp"
#include "Configuration.hpp"
#include "ConfigurationInfo.hpp"
#include "CoordinateVector.hpp"
#include "DustPhotonShootJobMarket.hpp"
#include "DustScattering.hpp"
#include "Log.hpp"
#include "MonochromaticPhotonSourceSpectrum.hpp"
#include "ParameterFile.hpp"
#include "PhotonSource.hpp"
#include "SimulationBox.hpp"
#include "SpiralGalaxyContinuousPhotonSource.hpp"
#include "SpiralGalaxyDensityFunction.hpp"
#include "Timer.hpp"
#include "VernerCrossSections.hpp"
#include "WorkDistributor.hpp"
#include "WorkEnvironment.hpp"

#include <iostream>
#include <string>

/**
 * @brief Perform a dusty radiative transfer simulation.
 *
 * This method reads the following parameters from the parameter file:
 *  - random seed: Seed for the random number generator (default: 42)
 *  - output folder: Folder where all output files will be placed (default: .)
 *  - number of photons: Number of photons to use (default: 5e5)
 *
 * @param parser CommandLineParser that contains the parsed command line
 * arguments.
 * @param write_output Flag indicating whether this process writes output.
 * @param programtimer Total program timer.
 * @param log Log to write logging info to.
 * @return Exit code: 0 on success.
 */
int DustSimulation::do_simulation(CommandLineParser &parser, bool write_output,
                                  Timer &programtimer, Log *log = nullptr) {

  // set the maximum number of openmp threads
  WorkEnvironment::set_max_num_threads(
      parser.get_value< int_fast32_t >("threads"));

  // third: read in the parameters of the run from a parameter file. This file
  // should be read by a ParameterFileParser object that acts as a dictionary
  ParameterFile params(parser.get_value< std::string >("params"));

  // fourth: construct the density grid. This should be stored in a separate
  // DensityGrid object with geometrical and physical properties
  SpiralGalaxyDensityFunction density_function(params, log);

  const SimulationBox simulation_box(params);
  CartesianDensityGrid grid(simulation_box, params, false, log);

  int_fast32_t random_seed =
      params.get_value< int_fast32_t >("DustSimulation:random seed", 42);

  SpiralGalaxyContinuousPhotonSource continuoussource(simulation_box.get_box(),
                                                      params, log);
  MonochromaticPhotonSourceSpectrum continuousspectrum(13.6, 1., log);

  Abundances abundances(0., 0., 0., 0., 0., 0., log);
  VernerCrossSections cross_sections;
  PhotonSource source(nullptr, nullptr, &continuoussource, &continuousspectrum,
                      abundances, cross_sections, false, log);

  DustScattering dust_scattering(params, log);

  // set up output
  std::string output_folder = Utilities::get_absolute_path(
      params.get_value< std::string >("DustSimulation:output folder", "."));
  CCDImage dust_image(output_folder, params, log);

  uint_fast64_t numphoton = params.get_value< uint_fast64_t >(
      "DustSimulation:number of photons", 5e5);

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output) {
    std::string pfilename = output_folder + "/dust-parameters-usedvalues.param";
    std::ofstream pfile(pfilename);
    params.print_contents(pfile);
    pfile.close();
    if (log) {
      log->write_status("Wrote used parameters to ", pfilename, ".");
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
  density_function.initialize();
  if (log) {
    log->write_status("Done.");
  }

  // done writing file, now initialize grid
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  // object used to distribute jobs in a shared memory parallel context
  WorkDistributor< DustPhotonShootJobMarket, DustPhotonShootJob >
      dust_workdistributor(parser.get_value< int_fast32_t >("threads"));
  const int_fast32_t worksize = dust_workdistributor.get_worksize();
  Timer worktimer;

  if (log) {
    log->write_status("Program will use ",
                      dust_workdistributor.get_worksize_string(),
                      " for photon shooting.");
  }
  DustPhotonShootJobMarket dustphotonshootjobs(
      source, dust_scattering, random_seed, grid, 0, dust_image, 100, worksize);

  if (log) {
    log->write_status("Start shooting ", numphoton, " photons...");
  }

  dustphotonshootjobs.set_numphoton(numphoton);
  worktimer.start();
  dust_workdistributor.do_in_parallel(dustphotonshootjobs);
  worktimer.stop();
  dustphotonshootjobs.update_image(dust_image);

  if (log) {
    log->write_status("Done shooting photons.");
  }

  if (log) {
    log->write_status("Saving final image...");
  }
  dust_image.save(1. / numphoton);
  if (log) {
    log->write_status("Done saving image.");
  }

  programtimer.stop();
  if (log) {
    log->write_status("Total program time: ",
                      Utilities::human_readable_time(programtimer.value()),
                      ".");
    log->write_status("Total photon shooting time: ",
                      Utilities::human_readable_time(worktimer.value()), ".");
  }

  return 0;
}
