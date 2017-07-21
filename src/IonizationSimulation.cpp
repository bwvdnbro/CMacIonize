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
 * @file IonizationSimulation.cpp
 *
 * @brief IonizationSimulation implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "IonizationSimulation.hpp"
#include "ChargeTransferRates.hpp"
#include "ContinuousPhotonSourceFactory.hpp"
#include "DensityFunctionFactory.hpp"
#include "DensityGridFactory.hpp"
#include "DensityGridWriterFactory.hpp"
#include "DensityMaskFactory.hpp"
#include "IonizationPhotonShootJobMarket.hpp"
#include "IonizationStateCalculator.hpp"
#include "LineCoolingData.hpp"
#include "MPICommunicator.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistributionFactory.hpp"
#include "PhotonSourceSpectrumFactory.hpp"
#include "TemperatureCalculator.hpp"
#include "VernerCrossSections.hpp"
#include "VernerRecombinationRates.hpp"
#include "WorkEnvironment.hpp"
#include <fstream>

/**
 * @brief Constructor.
 *
 * @param write_output Should this process write output?
 * @param every_iteration_output Write an output file after every iteration of
 * the algorithm?
 * @param num_threads Number of shared memory parallel threads to use.
 * @param parameterfile Name of the parameter file to use.
 * @param dry_run Do a dry_run?
 * @param comm MPICommunicator to use for distributed memory communications.
 * @param log Log to write logging info to.
 */
IonizationSimulation::IonizationSimulation(const bool write_output,
                                           const bool every_iteration_output,
                                           const int num_threads,
                                           const std::string parameterfile,
                                           const bool dry_run,
                                           MPICommunicator &comm, Log *log)
    : _log(log) {

  // set the maximum number of openmp threads
  WorkEnvironment::set_max_num_threads(num_threads);

  // second: initialize the parameters that are read in from static files
  // these files should be configured by CMake and put in a location that is
  // stored in a CMake configured header
  LineCoolingData line_cooling_data;

  // third: read in the parameters of the run from a parameter file. This file
  // should be read by a ParameterFile object that acts as a dictionary
  ParameterFile params(parameterfile);

  // fourth: construct the density grid. This should be stored in a separate
  // DensityGrid object with geometrical and physical properties
  DensityFunction *density_function =
      DensityFunctionFactory::generate(params, log);
  DensityMask *density_mask = DensityMaskFactory::generate(params, log);
  VernerCrossSections cross_sections;
  VernerRecombinationRates recombination_rates;

  DensityGrid *grid = DensityGridFactory::generate(params, log);

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
    std::ofstream pfile(folder + "/parameters-usedvalues.param");
    params.print_contents(pfile);
    pfile.close();
    if (log) {
      log->write_status("Wrote used parameters to ", folder,
                        "/parameters-usedvalues.param.");
    }
  }

  if (dry_run) {
    if (log) {
      log->write_warning("Dry run requested. Program will now halt.");
    }
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
      comm.distribute_block(0, grid->get_number_of_cells());
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
      workdistributor(num_threads);
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
  IonizationPhotonShootJobMarket photonshootjobs(source, random_seed, *grid, 0,
                                                 100, worksize);

  if (write_output) {
    writer->write(0, params);
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
      log->write_status(100. * typecount[PHOTONTYPE_ABSORBED] / totweight,
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

  // write snapshot
  if (write_output) {
    writer->write(nloop, params);
  }

  if (log) {
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
  delete temperature_calculator;
  delete continuousspectrum;
  delete spectrum;
}
