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
 * @file timeSPHArrayInterface.cpp
 *
 * @brief Timing test for the SPHArrayInterface.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "CMILibrary.hpp"
#include "TimingTools.hpp"
#include "YAMLDictionary.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>

/**
 * @brief Add the parameters to the parameter list.
 *
 * @param do_lloyd Whether or not to do Lloyd iterations.
 * @param npart Number of particles.
 * @param generator_file Name of the file containing the Voronoi generator
 * positions.
 * @param parameters Parameter list to add to.
 */
inline void add_parameters(const bool do_lloyd, const uint_fast32_t npart,
                           const std::string generator_file,
                           YAMLDictionary &parameters) {

  // Abundances
  parameters.add_value("Abundances:carbon", "0.");
  parameters.add_value("Abundances:helium", "0.");
  parameters.add_value("Abundances:neon", "0.");
  parameters.add_value("Abundances:nitrogen", "0.");
  parameters.add_value("Abundances:oxygen", "0.");
  parameters.add_value("Abundances:sulphur", "0.");
  // TemperatureCalculator
  parameters.add_value("TemperatureCalculator:do temperature calculation",
                       "false");
  // PhotonSource
  parameters.add_value("PhotonSource:diffuse field", "false");
  // SimulationBox
  parameters.add_value("SimulationBox:anchor", "[0. m, 0. m, 0. m]");
  parameters.add_value("SimulationBox:sides", "[1. m, 1. m, 1. m]");
  // DensityGrid
  parameters.add_value("DensityGrid:type", "Voronoi");
  parameters.add_value("DensityGrid:grid type", "Old");
  if (do_lloyd) {
    parameters.add_value("DensityGrid:number of Lloyd iterations", "5");
  }
  // DensityGrid:VoronoiGeneratorDistribution
  std::stringstream paramstream;
  paramstream << npart;
  parameters.add_value(
      "DensityGrid:VoronoiGeneratorDistribution:number of positions",
      paramstream.str());
  parameters.add_value("DensityGrid:VoronoiGeneratorDistribution:type", "SPH");
  parameters.add_value("DensityGrid:VoronoiGeneratorDistribution:filename",
                       generator_file);
  // DensityFunction
  parameters.add_value("DensityFunction:type", "Homogeneous");
  // PhotonSourceDistribution
  parameters.add_value("PhotonSourceDistribution:type", "SingleStar");
  parameters.add_value("PhotonSourceDistribution:position",
                       "[0.5 m, 0.5 m, 0.5 m]");
  parameters.add_value("PhotonSourceDistribution:luminosity", "1.e19 s^-1");
  // PhotonSourceSpectrum
  parameters.add_value("PhotonSourceSpectrum:type", "Monochromatic");
  parameters.add_value("PhotonSourceSpectrum:frequency", "13.6 eV");
  // IonizationSimulation
  parameters.add_value("IonizationSimulation:number of iterations", "10");
  parameters.add_value("IonizationSimulation:number of photons", "1e5");
  // DensityGridWriter
  parameters.add_value("DensityGridWriter:type", "AsciiFile");
  parameters.add_value("DensityGridWriter:prefix", "output");
  parameters.add_value("DensityGridWriter:padding", "3");
  // CrossSections
  parameters.add_value("CrossSections:type", "FixedValue");
  parameters.add_value("CrossSections:hydrogen_0", "6.3e-18 cm^2");
  parameters.add_value("CrossSections:helium_0", "0. m^2");
  parameters.add_value("CrossSections:carbon_1", "0. m^2");
  parameters.add_value("CrossSections:carbon_2", "0. m^2");
  parameters.add_value("CrossSections:nitrogen_0", "0. m^2");
  parameters.add_value("CrossSections:nitrogen_1", "0. m^2");
  parameters.add_value("CrossSections:nitrogen_2", "0. m^2");
  parameters.add_value("CrossSections:oxygen_0", "0. m^2");
  parameters.add_value("CrossSections:oxygen_1", "0. m^2");
  parameters.add_value("CrossSections:neon_0", "0. m^2");
  parameters.add_value("CrossSections:neon_1", "0. m^2");
  parameters.add_value("CrossSections:sulphur_1", "0. m^2");
  parameters.add_value("CrossSections:sulphur_2", "0. m^2");
  parameters.add_value("CrossSections:sulphur_3", "0. m^2");
  // RecombinationRates
  parameters.add_value("RecombinationRates:type", "FixedValue");
  parameters.add_value("RecombinationRates:hydrogen_1", "2.7e-13 cm^3 s^-1");
  parameters.add_value("RecombinationRates:helium_1", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:carbon_2", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:carbon_3", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:nitrogen_1", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:nitrogen_2", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:nitrogen_3", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:oxygen_1", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:oxygen_2", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:neon_1", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:neon_2", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:sulphur_2", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:sulphur_3", "0. m^3 s^-1");
  parameters.add_value("RecombinationRates:sulphur_4", "0. m^3 s^-1");
}

/**
 * @brief Timing test for the SPHArrayInterface.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  timingtools_init("timeSPHArrayInterface", argc, argv);

  const uint_fast32_t numpart = 1000;
  // create input arrays
  std::vector< double > x(numpart, 0.);
  std::vector< double > y(numpart, 0.);
  std::vector< double > z(numpart, 0.);
  std::vector< double > h(numpart, 0.);
  std::vector< double > m(numpart, 0.);
  std::vector< double > xH(numpart, 0.);
  for (size_t i = 0; i < numpart; ++i) {
    x[i] = Utilities::random_double();
    y[i] = Utilities::random_double();
    z[i] = Utilities::random_double();
    h[i] = 0.2;
    m[i] = 0.001;
  }

  // create text file with cell generator positions
  std::ofstream posfile("generator_positions.txt");
  for (uint_fast32_t i = 0; i < numpart; ++i) {
    posfile << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
  }
  posfile.close();

  timingtools_start_timing_block("M over V mapping") {
    timingtools_start_timing();

    // create parameter file: generator position file, whether or not Lloyd,
    // type of mapping
    YAMLDictionary parameters;
    add_parameters(false, numpart, "generator_positions.txt", parameters);
    std::ofstream parfile("parameters.param");
    parameters.print_contents(parfile);
    parfile.close();

    // call init
    cmi_init("parameters.param", timingtools_num_threads, 1., 1., "M_over_V",
             false);

    // call library function
    cmi_compute_neutral_fraction_dp(&x[0], &y[0], &z[0], &h[0], &m[0], &xH[0],
                                    numpart);

    // clean up library
    cmi_destroy();

    std::rename("ionization-simulation-time-log.txt", "M-over-V-time-log.txt");

    timingtools_stop_timing();
  }
  timingtools_end_timing_block("M over V mapping");

  timingtools_start_timing_block("centroid mapping, no Lloyd iterations") {
    timingtools_start_timing();

    // create parameter file: generator position file, whether or not Lloyd,
    // type of mapping
    YAMLDictionary parameters;
    add_parameters(false, numpart, "generator_positions.txt", parameters);
    std::ofstream parfile("parameters.param");
    parameters.print_contents(parfile);
    parfile.close();

    // call init
    cmi_init("parameters.param", timingtools_num_threads, 1., 1., "centroid",
             false);

    // call library function
    cmi_compute_neutral_fraction_dp(&x[0], &y[0], &z[0], &h[0], &m[0], &xH[0],
                                    numpart);

    // clean up library
    cmi_destroy();

    std::rename("ionization-simulation-time-log.txt",
                "centroid-no-Lloyd-time-log.txt");

    timingtools_stop_timing();
  }
  timingtools_end_timing_block("centroid mapping, no Lloyd iterations");

  timingtools_start_timing_block("centroid mapping, Lloyd iterations") {
    timingtools_start_timing();

    // create parameter file: generator position file, whether or not Lloyd,
    // type of mapping
    YAMLDictionary parameters;
    add_parameters(true, numpart, "generator_positions.txt", parameters);
    std::ofstream parfile("parameters.param");
    parameters.print_contents(parfile);
    parfile.close();

    // call init
    cmi_init("parameters.param", timingtools_num_threads, 1., 1., "centroid",
             false);

    // call library function
    cmi_compute_neutral_fraction_dp(&x[0], &y[0], &z[0], &h[0], &m[0], &xH[0],
                                    numpart);

    // clean up library
    cmi_destroy();

    std::rename("ionization-simulation-time-log.txt",
                "centroid-Lloyd-time-log.txt");

    timingtools_stop_timing();
  }
  timingtools_end_timing_block("centroid mapping, Lloyd iterations");

  timingtools_start_timing_block("Petkova mapping, no Lloyd iterations") {
    timingtools_start_timing();

    // create parameter file: generator position file, whether or not Lloyd,
    // type of mapping
    YAMLDictionary parameters;
    add_parameters(false, numpart, "generator_positions.txt", parameters);
    std::ofstream parfile("parameters.param");
    parameters.print_contents(parfile);
    parfile.close();

    // call init
    cmi_init("parameters.param", timingtools_num_threads, 1., 1., "Petkova",
             false);

    // call library function
    cmi_compute_neutral_fraction_dp(&x[0], &y[0], &z[0], &h[0], &m[0], &xH[0],
                                    numpart);

    // clean up library
    cmi_destroy();

    std::rename("ionization-simulation-time-log.txt",
                "Petkova-no-Lloyd-time-log.txt");

    timingtools_stop_timing();
  }
  timingtools_end_timing_block("Petkova mapping, no Lloyd iterations");

  timingtools_start_timing_block("Petkova mapping, Lloyd iterations") {
    timingtools_start_timing();

    // create parameter file: generator position file, whether or not Lloyd,
    // type of mapping
    YAMLDictionary parameters;
    add_parameters(true, numpart, "generator_positions.txt", parameters);
    std::ofstream parfile("parameters.param");
    parameters.print_contents(parfile);
    parfile.close();

    // call init
    cmi_init("parameters.param", timingtools_num_threads, 1., 1., "Petkova",
             false);

    // call library function
    cmi_compute_neutral_fraction_dp(&x[0], &y[0], &z[0], &h[0], &m[0], &xH[0],
                                    numpart);

    // clean up library
    cmi_destroy();

    std::rename("ionization-simulation-time-log.txt",
                "Petkova-Lloyd-time-log.txt");

    timingtools_stop_timing();
  }
  timingtools_end_timing_block("Petkova mapping, Lloyd iterations");

  return 0;
}
