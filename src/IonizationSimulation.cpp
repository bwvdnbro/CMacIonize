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
 * @param num_thread Number of shared memory parallel threads to use.
 * @param parameterfile Name of the parameter file to use.
 * @param mpi_communicator MPICommunicator to use for distributed memory
 * communications.
 * @param log Log to write logging info to.
 */
IonizationSimulation::IonizationSimulation(const bool write_output,
                                           const bool every_iteration_output,
                                           const int num_thread,
                                           const std::string parameterfile,
                                           MPICommunicator *mpi_communicator,
                                           Log *log)
    : _num_thread(WorkEnvironment::set_max_num_threads(num_thread)),
      _every_iteration_output(every_iteration_output),
      _mpi_communicator(mpi_communicator), _log(log),
      _work_distributor(_num_thread), _parameter_file(parameterfile),
      _number_of_iterations(_parameter_file.get_value< unsigned int >(
          "max_number_iterations", 10)),
      _number_of_photons(
          _parameter_file.get_value< unsigned int >("number of photons", 100)),
      _number_of_photons_init(_parameter_file.get_value< unsigned int >(
          "number of photons init", _number_of_photons)),
      _abundances(_parameter_file, _log) {

  if (_log) {
    if (_num_thread == 1) {
      _log->write_status("IonizationSimulation will use 1 thread.");
    } else {
      _log->write_status("IonizationSimulation will use ", _num_thread,
                         " threads.");
    }
  }

  // create the density grid and related objects
  _density_function = DensityFunctionFactory::generate(_parameter_file, _log);
  _density_mask = DensityMaskFactory::generate(_parameter_file, _log);

  _density_grid = DensityGridFactory::generate(_parameter_file, _log);

  // create the discrete UV sources
  _photon_source_distribution =
      PhotonSourceDistributionFactory::generate(_parameter_file, _log);
  _photon_source_spectrum = PhotonSourceSpectrumFactory::generate(
      "photonsourcespectrum", _parameter_file, _log);

  // sanity checks on discrete sources
  if (_photon_source_distribution != nullptr &&
      _photon_source_spectrum == nullptr) {
    cmac_error("No spectrum provided for the discrete photon sources!");
  }
  if (_photon_source_distribution == nullptr &&
      _photon_source_spectrum != nullptr) {
    cmac_warning("Discrete photon source spectrum provided, but no discrete "
                 "photon source distributions. The given spectrum will be "
                 "ignored.");
  }

  // create the continuous UV sources
  _continuous_photon_source =
      ContinuousPhotonSourceFactory::generate(_parameter_file, _log);
  _continuous_photon_source_spectrum = PhotonSourceSpectrumFactory::generate(
      "continuousphotonsourcespectrum", _parameter_file, _log);

  // sanity checks on continuous sources
  if (_continuous_photon_source != nullptr &&
      _continuous_photon_source_spectrum == nullptr) {
    cmac_error("No spectrum provided for the continuous photon sources!");
  }
  if (_continuous_photon_source == nullptr &&
      _continuous_photon_source_spectrum != nullptr) {
    cmac_warning("Continuous photon source spectrum provided, but no "
                 "continuous photon source. The given spectrum will be "
                 "ignored.");
  }

  // create the actual photon source objects that emits the UV photons
  _photon_source = new PhotonSource(
      _photon_source_distribution, _photon_source_spectrum,
      _continuous_photon_source, _continuous_photon_source_spectrum,
      _abundances, _cross_sections, _log);
  const double total_luminosity = _photon_source->get_total_luminosity();

  // set up output
  _density_grid_writer = nullptr;
  if (write_output) {
    _density_grid_writer = DensityGridWriterFactory::generate(
        _parameter_file, *_density_grid, _log);
  }

  // computation objects

  // used to calculate the ionization state at fixed temperature
  _ionization_state_calculator = new IonizationStateCalculator(
      total_luminosity, _abundances, _recombination_rates,
      _charge_transfer_rates);

  bool calculate_temperature =
      _parameter_file.get_value< bool >("calculate_temperature", true);

  _temperature_calculator = nullptr;
  if (calculate_temperature) {
    // used to calculate both the ionization state and the temperature
    _temperature_calculator = new TemperatureCalculator(
        total_luminosity, _abundances,
        _parameter_file.get_value< double >("pahfac", 1.),
        _parameter_file.get_value< double >("crfac", 0.),
        _parameter_file.get_value< double >("crlim", 0.75),
        _parameter_file.get_physical_value< QUANTITY_LENGTH >("crscale",
                                                              "1.33333 kpc"),
        _line_cooling_data, _recombination_rates, _charge_transfer_rates, _log);
  }

  // create ray tracing objects
  int random_seed = _parameter_file.get_value< int >("random_seed", 42);
  // make sure every thread on every process has another random seed
  random_seed += _mpi_communicator->get_rank() * _num_thread;
  _ionization_photon_shoot_job_market = new IonizationPhotonShootJobMarket(
      *_photon_source, random_seed, *_density_grid, 0, 100, _num_thread);

  // we are done reading the parameter file
  // now output all parameters (also those for which default values were used)
  // to a reference parameter file (only rank 0 does this)
  if (write_output) {
    std::string folder =
        Utilities::get_absolute_path(_parameter_file.get_value< std::string >(
            "densitygridwriter:folder", "."));
    std::ofstream pfile(folder + "/parameters-usedvalues.param");
    _parameter_file.print_contents(pfile);
    pfile.close();
    if (_log) {
      _log->write_status("Wrote used parameters to ", folder,
                         "/parameters-usedvalues.param.");
    }
  }
}

/**
 * @brief Initialize the simulation.
 *
 * @param density_function DensityFunction to use. If no DensityFunction is
 * given, the internal DensityFunction is used.
 */
void IonizationSimulation::initialize(DensityFunction *density_function) {

  if (density_function == nullptr) {
    density_function = _density_function;
  }

  // initialize the density function
  if (_log) {
    _log->write_status("Initializing DensityFunction...");
  }
  density_function->initialize();
  if (_log) {
    _log->write_status("Done.");
  }

  // initialize the actual grid
  std::pair< unsigned long, unsigned long > block =
      _mpi_communicator->distribute_block(0,
                                          _density_grid->get_number_of_cells());
  _density_grid->initialize(block, *density_function);

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

  // if necessary, initialize and apply the density mask
  if (_density_mask != nullptr) {
    if (_log) {
      _log->write_status("Initializing DensityMask...");
    }
    _density_mask->initialize();
    if (_log) {
      _log->write_status("Done initializing mask. Applying mask...");
    }
    _density_mask->apply(*_density_grid);
    if (_log) {
      _log->write_status("Done applying mask.");
    }
  }
}

/**
 * @brief Run the actual simulation.
 */
void IonizationSimulation::run() {
  // write the initial state of the grid to an output file
  if (_density_grid_writer) {
    _density_grid_writer->write(0, _parameter_file);
  }

  std::pair< unsigned long, unsigned long > block =
      _mpi_communicator->distribute_block(0,
                                          _density_grid->get_number_of_cells());
  // finally: the actual program loop whereby the density grid is ray traced
  // using photon packets generated by the stellar sources
  unsigned int loop = 0;
  while (loop < _number_of_iterations) {

    if (_log) {
      _log->write_status("Starting loop ", loop, ".");
    }

    unsigned int lnumphoton = _number_of_photons;

    if (loop == 0) {
      // overwrite the number of photons for the first loop (might be useful
      // if more than 1 boundary is periodic, since the initial neutral
      // fractions are very low)
      lnumphoton = _number_of_photons_init;
    }

    _density_grid->reset_grid(*_density_function);
    if (_log) {
      _log->write_status("Start shooting ", lnumphoton, " photons...");
    }

    double typecount[PHOTONTYPE_NUMBER] = {0};

    double totweight = 0.;

    unsigned int local_numphoton = lnumphoton;

    // make sure this process does only part of the total number of photons
    local_numphoton = _mpi_communicator->distribute(local_numphoton);

    _ionization_photon_shoot_job_market->set_numphoton(local_numphoton);
    _work_timer.start();
    _work_distributor.do_in_parallel(*_ionization_photon_shoot_job_market);
    _work_timer.stop();

    _ionization_photon_shoot_job_market->update_counters(totweight, typecount);

    // make sure the total weight and typecount is reduced across all
    // processes
    _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES >(totweight);
    _mpi_communicator->reduce< MPI_SUM_OF_ALL_PROCESSES, PHOTONTYPE_NUMBER >(
        typecount);

    if (_log) {
      _log->write_status("Done shooting photons.");
      _log->write_status(
          100. * typecount[PHOTONTYPE_ABSORBED] / totweight,
          "% of photons were reemitted as non-ionizing photons.");
      _log->write_status(100. * (typecount[PHOTONTYPE_DIFFUSE_HI] +
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
      _log->write_status("Escape fraction: ", escape_fraction, "%.");
      double escape_fraction_HI =
          (100. * typecount[PHOTONTYPE_DIFFUSE_HI]) / totweight;
      _log->write_status("Diffuse HI escape fraction: ", escape_fraction_HI,
                         "%.");
      double escape_fraction_HeI =
          (100. * typecount[PHOTONTYPE_DIFFUSE_HeI]) / totweight;
      _log->write_status("Diffuse HeI escape fraction: ", escape_fraction_HeI,
                         "%.");
    }

    if (_log) {
      _log->write_status("Calculating ionization state after shooting ",
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

    if (_temperature_calculator && loop > 3) {
      _temperature_calculator->calculate_temperature(totweight, *_density_grid,
                                                     block);
    } else {
      _ionization_state_calculator->calculate_ionization_state(
          totweight, *_density_grid, block);
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

    if (_log) {
      _log->write_status("Done calculating ionization state.");
    }

    // calculate emissivities
    // we disabled this, since we now have the post-processing Python library
    // for this
    //    if (loop > 3 && abundances.get_abundance(ELEMENT_He) > 0.) {
    //      emissivity_calculator.calculate_emissivities(*grid);
    //    }

    ++loop;

    if (_density_grid_writer && _every_iteration_output &&
        loop < _number_of_iterations) {
      _density_grid_writer->write(loop, _parameter_file);
    }
  }

  if (_log && loop == _number_of_iterations) {
    _log->write_status("Maximum number of iterations (", _number_of_iterations,
                       ") reached, stopping.");
  }

  // write final snapshot
  if (_density_grid_writer) {
    _density_grid_writer->write(_number_of_iterations, _parameter_file);
  }
}

/**
 * @brief Destructor.
 *
 * Deallocate memory used by internal variables.
 */
IonizationSimulation::~IonizationSimulation() {
  if (_log) {
    _log->write_status("Total photon shooting time: ",
                       Utilities::human_readable_time(_work_timer.value()),
                       ".");
  }

  // we delete the objects in the opposite order in which they were created
  // note that we do not check for nullptrs, as deleting a nullptr is allowed
  // and won't do anything

  // ray tracing objects
  delete _ionization_photon_shoot_job_market;

  // computation objects
  delete _temperature_calculator;
  delete _ionization_state_calculator;

  // snapshot output
  delete _density_grid_writer;

  // actual photon source object
  delete _photon_source;

  // continuous sources
  delete _continuous_photon_source_spectrum;
  delete _continuous_photon_source;

  // discrete sources
  delete _photon_source_spectrum;
  delete _photon_source_distribution;

  // density grid and related objects
  delete _density_grid;
  delete _density_mask;
  delete _density_function;
}
