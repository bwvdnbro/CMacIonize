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
 * @file EmissivityCalculationSimulation.cpp
 *
 * @brief Program to compute the emissivities based on a CMacIonize snapshot:
 * implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "EmissivityCalculationSimulation.hpp"
#include "AbundanceModelFactory.hpp"
#include "Abundances.hpp"
#include "CommandLineParser.hpp"
#include "EmissivityCalculator.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Timer.hpp"
#include "WorkEnvironment.hpp"

/**
 * @brief Add program specific command line parameters.
 *
 * @param parser CommandLineParser.
 */
void EmissivityCalculationSimulation::add_command_line_parameters(
    CommandLineParser &parser) {
  parser.add_option< std::string >("file", 'f',
                                   "Name of the input/output file.", "");
}

/**
 * @brief Do the simulation.
 *
 * @param parser CommandLineParser.
 * @param write_output Should output be written?
 * @param programtimer Total program time timer.
 * @param log Log to write logging info to.
 * @return Exit code: 0 on success.
 */
int EmissivityCalculationSimulation::do_simulation(CommandLineParser &parser,
                                                   bool write_output,
                                                   Timer &programtimer,
                                                   Log *log) {

  // set the maximum number of openmp threads
  WorkEnvironment::set_max_num_threads(
      parser.get_value< int_fast32_t >("threads"));
  const std::string parameterfile_name =
      parser.get_value< std::string >("params");
  ParameterFile params(parameterfile_name);
  bool do_line[NUMBER_OF_EMISSIONLINES];
  for (int_fast32_t i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
    do_line[i] = params.get_value< bool >(
        "EmissivityValues:" + EmissivityValues::get_name(i), false);
  }

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output) {
    const std::string usedvaluename = parameterfile_name + ".used-values";
    std::ofstream pfile(usedvaluename);
    params.print_contents(pfile);
    pfile.close();
    if (log) {
      log->write_status("Wrote used parameters to ", usedvaluename, ".");
    }
  }

  const std::string input_file_name = parser.get_value< std::string >("file");
  if (input_file_name == "") {
    cmac_error("No input file name provided (--file)!");
  }

  if (log) {
    log->write_status("Reading file \"", input_file_name, "\"...");
  }

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(input_file_name, HDF5Tools::HDF5FILEMODE_APPEND);

  if (log) {
    log->write_status("Reading simulation parameters...");
  }

  // read grid parameters
  HDF5Tools::HDF5Group group = HDF5Tools::open_group(file, "/Parameters");
  std::vector< std::string > parameternames =
      HDF5Tools::get_attribute_names(group);
  ParameterFile simulation_parameters;
  for (auto it = parameternames.begin(); it != parameternames.end(); ++it) {
    std::string attname = *it;
    std::string attvalue =
        HDF5Tools::read_attribute< std::string >(group, attname);
    simulation_parameters.add_value(attname, attvalue);
  }
  HDF5Tools::close_group(group);

  if (log) {
    log->write_status("Done reading parameters.");
  }

  // units
  double unit_number_density_in_SI = 1.;
  double unit_temperature_in_SI = 1.;
  if (HDF5Tools::group_exists(file, "/Units")) {

    if (log) {
      log->write_status("Reading simulation units...");
    }

    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Units");
    const double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit length in cgs (U_L)");
    const double unit_temperature_in_cgs = HDF5Tools::read_attribute< double >(
        units, "Unit temperature in cgs (U_T)");
    const double unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    unit_number_density_in_SI =
        1. / unit_length_in_SI / unit_length_in_SI / unit_length_in_SI;
    // K is K
    unit_temperature_in_SI = unit_temperature_in_cgs;
    HDF5Tools::close_group(units);

    if (log) {
      log->write_status("Done reading units.");
    }
  }

  if (log) {
    log->write_status("Setting up emissivity calculator...");
  }

  // make sure we read the abundances in the right way:
  //  - old parameter file: directly
  //  - new paramter file: using the AbundanceModel
  Abundances abundances;
  if (simulation_parameters.has_value("Abundances:helium")) {
    abundances = Abundances(simulation_parameters);
  } else {
    const AbundanceModel *abundance_model =
        AbundanceModelFactory::generate(params, log);
    abundances = abundance_model->get_abundances();
    delete abundance_model;
  }
  EmissivityCalculator calculator(abundances);

  if (log) {
    log->write_status("Done setting up.");
  }

  HDF5Tools::HDF5Group parttype0 = HDF5Tools::open_group(file, "/PartType0");

  if (log) {
    log->write_status("Creating output datasets...");
  }

  const CoordinateVector< uint_fast32_t > number_of_cells =
      simulation_parameters.get_value< CoordinateVector< uint_fast32_t > >(
          "DensityGrid:number of cells", CoordinateVector< uint_fast32_t >(-1));
  const uint_fast32_t total_number_of_cells =
      number_of_cells.x() * number_of_cells.y() * number_of_cells.z();

  for (int_fast32_t i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
    if (do_line[i]) {
      if (HDF5Tools::group_exists(parttype0, EmissivityValues::get_name(i))) {
        if (log) {
          log->write_warning("Dataset \"", EmissivityValues::get_name(i),
                             "\" already exists! Values will be overwritten!");
        }
      } else {
        HDF5Tools::create_dataset< double >(
            parttype0, EmissivityValues::get_name(i), total_number_of_cells);
      }
    }
  }

  if (log) {
    log->write_status("Done creating datasets.");
  }

  if (log) {
    log->write_status("Starting emissivity calculation...");
  }

  {
    std::vector< double > number_density =
        HDF5Tools::read_dataset< double >(parttype0, "NumberDensity");
    std::vector< double > temperature =
        HDF5Tools::read_dataset< double >(parttype0, "Temperature");
    const size_t size = number_density.size();
    cmac_assert(size == temperature.size());
    cmac_assert(size == total_number_of_cells);

    std::vector< std::vector< double > > neutral_fractions(NUMBER_OF_IONNAMES);
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      // skip ionic fractions that do not exist
      if (HDF5Tools::group_exists(parttype0,
                                  "NeutralFraction" + get_ion_name(i))) {
        neutral_fractions[i] = HDF5Tools::read_dataset< double >(
            parttype0, "NeutralFraction" + get_ion_name(i));
      } else {
        cmac_error("Missing ionic fractions for \"%s\"!",
                   get_ion_name(i).c_str());
      }
      cmac_assert(size == neutral_fractions[i].size());
    }

    // prepare output arrays
    std::vector< double > emissivities[NUMBER_OF_EMISSIONLINES];
    for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
      if (do_line[line]) {
        emissivities[line].resize(size);
      }
    }

#ifdef HAVE_OPENMP
#pragma omp parallel for default(shared)
#endif
    for (size_t i = 0; i < size; ++i) {
      IonizationVariables ionization_variables;
      ionization_variables.set_number_density(number_density[i] *
                                              unit_number_density_in_SI);
      ionization_variables.set_temperature(temperature[i] *
                                           unit_temperature_in_SI);
      for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
        ionization_variables.set_ionic_fraction(ion, neutral_fractions[ion][i]);
      }
      double output[NUMBER_OF_EMISSIONLINES];
      calculator.calculate_emissivities(ionization_variables, do_line, output);
      for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
        if (do_line[line]) {
          emissivities[line][i] = output[line];
        }
      }
    }
    for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
      if (do_line[line]) {
        HDF5Tools::append_dataset< double >(
            parttype0, EmissivityValues::get_name(line), 0, emissivities[line]);
      }
    }
  }

  if (log) {
    log->write_status("Finished emissivity calculation.");
  }

  HDF5Tools::close_group(parttype0);

  HDF5Tools::close_file(file);

  if (log) {
    log->write_status("Closed file.");
  }

  programtimer.stop();
  if (log) {
    log->write_status("Total program time: ",
                      Utilities::human_readable_time(programtimer.value()),
                      ".");
  }

  // we cannot delete the log, since it is still used in the destructor of
  // objects that are destructed at the return of the main program
  // this is not really a problem, as the memory is freed up by the OS anyway
  // delete log;

  return 0;
}
