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
#include "ContinuumEmission.hpp"
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

  const bool full_spectrum =
      params.get_value< bool >("Spectrum:compute full spectrum", false);
  const double minimum_wavelength =
      params.get_physical_value< QUANTITY_LENGTH >(
          "Spectrum:minimum wavelength", "1000. angstrom");
  const double maximum_wavelength =
      params.get_physical_value< QUANTITY_LENGTH >(
          "Spectrum:maximum wavelength", "2.e-4 m");
  const uint_fast32_t number_of_spectral_bins =
      params.get_value< uint_fast32_t >("Spectrum:number of bins", 1000u);

  bool compute_line[NUMBER_OF_EMISSIONLINES];
  bool output_line[NUMBER_OF_EMISSIONLINES];
  for (int_fast32_t i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
    output_line[i] = params.get_value< bool >(
        "EmissivityValues:" + EmissivityValues::get_name(i), false);
    // if we want the full spectrum, we need to compute all lines that are not
    // pseudo-lines
    compute_line[i] =
        output_line[i] ||
        (full_spectrum && EmissivityValues::get_central_wavelength(i) > 0.);
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
  //  - new parameter file: using the AbundanceModel
  Abundances abundances;
  if (simulation_parameters.has_value("Abundances:helium")) {
    abundances = Abundances(simulation_parameters);
  } else {
    const AbundanceModel *abundance_model =
        AbundanceModelFactory::generate(simulation_parameters, log);
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
    if (output_line[i]) {
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
      if (compute_line[line]) {
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
      calculator.calculate_emissivities(ionization_variables, compute_line,
                                        output);
      for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
        if (compute_line[line]) {
          emissivities[line][i] = output[line];
        }
      }
    }

    // output the lines
    for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
      if (output_line[line]) {
        HDF5Tools::append_dataset< double >(
            parttype0, EmissivityValues::get_name(line), 0, emissivities[line]);
      }
    }

    // output the full spectrum
    if (full_spectrum) {

      // set up an empty spectrum
      std::vector< double > spectrum(
          temperature.size() * number_of_spectral_bins, 0.);

      const double clight =
          PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED);

      // construct the wavelength bins (we store the bin edges in frequency
      // space, since this makes more sense
      std::vector< double > frequency_edges(number_of_spectral_bins + 1, 0.);
      std::vector< double > wavelength_edges(number_of_spectral_bins + 1, 0.);
      std::vector< double > wavelengths_mid(number_of_spectral_bins, 0.);
      const double log_minimum_wavelength = std::log10(minimum_wavelength);
      const double log_maximum_wavelength = std::log10(maximum_wavelength);
      const double minimum_frequency = clight / maximum_wavelength;
      const double maximum_frequency = clight / minimum_wavelength;
      const double dlog_wavelength =
          (log_maximum_wavelength - log_minimum_wavelength) /
          number_of_spectral_bins;
      wavelength_edges[0] = minimum_wavelength;
      frequency_edges[0] = clight / wavelength_edges[0];
      for (uint_fast32_t ibin = 0; ibin < number_of_spectral_bins; ++ibin) {
        wavelength_edges[ibin + 1] = std::pow(
            10., log_minimum_wavelength + (ibin + 1.) * dlog_wavelength);
        wavelengths_mid[ibin] =
            0.5 * (wavelength_edges[ibin] + wavelength_edges[ibin + 1]);
        frequency_edges[ibin + 1] = clight / wavelength_edges[ibin + 1];
      }
      // set up a sorted copy of the frequencies for later use
      std::vector< double > inverse_frequency_edges(frequency_edges);
      std::sort(inverse_frequency_edges.begin(), inverse_frequency_edges.end());

      // precompute the Doppler broadening factor
      const double doppler_factor = std::sqrt(
          2. *
          PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN) /
          (PhysicalConstants::get_physical_constant(
               PHYSICALCONSTANT_PROTON_MASS) *
           clight * clight));

      // loop over all lines
      for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
        const double lambda_line =
            EmissivityValues::get_central_wavelength(line);
        if (lambda_line > 0.) {
          const double nu_line = clight / lambda_line;
          for (uint_fast32_t i = 0; i < temperature.size(); ++i) {
            // delta_nu is the standard deviation of the Gaussian line profile
            const double delta_nu =
                doppler_factor * std::sqrt(temperature[i]) * nu_line;
            // we only include a 3 sigma interval around the central wavelength
            const double nu_min = nu_line - 3. * delta_nu;
            const double nu_max = nu_line + 3. * delta_nu;
            if (nu_max > minimum_frequency && nu_min < maximum_frequency) {
              // the line profile overlaps with the desired spectral range
              // add the line to the corresponding bins
              // first, figure out the indices of the overlapping bins
              const uint_fast32_t imin =
                  number_of_spectral_bins -
                  Utilities::locate(nu_min, &inverse_frequency_edges[0],
                                    number_of_spectral_bins + 1);
              const uint_fast32_t imax =
                  number_of_spectral_bins -
                  Utilities::locate(nu_max, &inverse_frequency_edges[0],
                                    number_of_spectral_bins + 1);
              cmac_assert(imin > 0);
              cmac_assert(imin < number_of_spectral_bins + 1);
              cmac_assert(imax > 0);
              cmac_assert(imax < number_of_spectral_bins + 1);
              if (imin == imax) {
                // the line completely overlaps with a single bin
                // add the full strength to the spectrum
                spectrum[i * number_of_spectral_bins + imin - 1] +=
                    emissivities[line][i];
              } else {
                // multiple bins overlap, loop over them
                // note that imin > imax, since the frequencies are sorted
                // backwards
                for (uint_fast32_t ibin = imax; ibin <= imin; ++ibin) {
                  // determine the fraction of the line strength to contribute
                  // to this bin
                  const double lower_limit =
                      std::max(frequency_edges[ibin], nu_min);
                  const double upper_limit =
                      std::min(frequency_edges[ibin - 1], nu_max);
                  const double cdf_min =
                      0.5 * (1. + std::erf((lower_limit - nu_line) /
                                           (std::sqrt(2.) * delta_nu)));
                  const double cdf_max =
                      0.5 * (1. + std::erf((upper_limit - nu_line) /
                                           (std::sqrt(2.) * delta_nu)));
                  spectrum[i * number_of_spectral_bins + ibin - 1] +=
                      (cdf_max - cdf_min) * emissivities[line][i];
                }
              }
            }
          }
        }
      }

      // normalise the line spectrum by dividing by the wavelength bin size
      for (uint_fast32_t ibin = 0; ibin < number_of_spectral_bins; ++ibin) {
        const double dlambda =
            wavelength_edges[ibin + 1] - wavelength_edges[ibin];
        for (uint_fast32_t i = 0; i < temperature.size(); ++i) {
          spectrum[i * number_of_spectral_bins + ibin] /= dlambda;
        }
      }

      // add the continuum
      const double pifac = 0.25 / M_PI;
      ContinuumEmission continuum_emission;
      for (uint_fast32_t ibin = 0; ibin < number_of_spectral_bins; ++ibin) {
        const double lambda = wavelengths_mid[ibin];
        const double dlambda =
            wavelength_edges[ibin + 1] - wavelength_edges[ibin];
        const double dnu = frequency_edges[ibin] - frequency_edges[ibin + 1];
        cmac_assert(dnu > 0.);
        for (uint_fast32_t i = 0; i < temperature.size(); ++i) {
          const double T = temperature[i];
          const double nH = number_density[i];
          const double xH = neutral_fractions[ION_H_n][i];
          const double xHe = neutral_fractions[ION_He_n][i];
          const double nHp = (1. - xH) * nH;
          const double nHep =
              (1. - xHe) * abundances.get_abundance(ELEMENT_He) * nH;
          const double ne = nHp + nHep;
          // convert Hz^-1 to m^-1 by multiplying with dnu/dlambda
          spectrum[i * number_of_spectral_bins + ibin] +=
              pifac * ne * nHp *
              (continuum_emission.gamma_HI(lambda, T) +
               ContinuumEmission::gamma_2q(lambda, T, ne, nHp)) *
              dnu / dlambda;
        }
      }

      if (HDF5Tools::group_exists(parttype0, "Spectrum wavelengths")) {
        if (log) {
          log->write_warning("Dataset \"Spectrum wavelengths\" already exists! "
                             "Values will be overwritten!");
        }
      } else {
        HDF5Tools::create_dataset< double >(parttype0, "Spectrum wavelengths",
                                            number_of_spectral_bins, true);
      }
      HDF5Tools::append_dataset(parttype0, "Spectrum wavelengths", 0,
                                wavelengths_mid);
      if (HDF5Tools::group_exists(parttype0, "Spectrum")) {
        if (log) {
          log->write_warning("Dataset \"Spectrum\" already exists! Values will "
                             "be overwritten!");
        }
      } else {
        HDF5Tools::create_datatable< double >(parttype0, "Spectrum",
                                              temperature.size(),
                                              number_of_spectral_bins, true);
      }
      for (uint_fast32_t i = 0; i < temperature.size(); ++i) {
        std::vector< double > row_copy(
            spectrum.begin() + i * number_of_spectral_bins,
            spectrum.begin() + (i + 1) * number_of_spectral_bins);
        HDF5Tools::fill_row< double >(parttype0, "Spectrum", i, row_copy);
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
