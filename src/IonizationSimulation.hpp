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
 * @file IonizationSimulation.hpp
 *
 * @brief Ionization radiative transfer simulation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONSIMULATION_HPP
#define IONIZATIONSIMULATION_HPP

#include "Abundances.hpp"
#include "ChargeTransferRates.hpp"
#include "IonizationPhotonShootJobMarket.hpp"
#include "LineCoolingData.hpp"
#include "ParameterFile.hpp"
#include "Timer.hpp"

#include <string>

class ContinuousPhotonSource;
class CrossSections;
class DensityFunction;
class DensityGrid;
class DensityGridWriter;
class DensityMask;
class Log;
class MPICommunicator;
class PhotonSource;
class PhotonSourceDistribution;
class PhotonSourceSpectrum;
class RecombinationRates;
class TemperatureCalculator;

/**
 * @brief Ionization radiative transfer simulation.
 */
class IonizationSimulation {
private:
  /// simulation wide parameters

  /*! @brief Number of shared memory parallel threads available for use. */
  const int_fast32_t _num_thread;

  /*! @brief Should the simulation write a snapshot file for every iteration of
   *  the algorithm? */
  const bool _every_iteration_output;

  /*! @brief Should the simulation output statistical information about the
   *  photons? */
  const bool _output_statistics;

  /// objects provided to the constructor, they are not owned by the simulation

  /*! @brief MPI communicator object. */
  MPICommunicator *_mpi_communicator;

  /*! @brief Log to write logging info to. */
  Log *_log;

  /// parameter-less objects owned by the simulation

  /*! @brief Object used to distribute jobs in a shared memory parallel
   *  context. */
  WorkDistributor< IonizationPhotonShootJobMarket, IonizationPhotonShootJob >
      _work_distributor;

  /*! @brief Data values for line cooling. */
  const LineCoolingData _line_cooling_data;

  /*! @brief Charge transfer rates. */
  const ChargeTransferRates _charge_transfer_rates;

  /// parameter file

  /*! @brief Parameter file. */
  ParameterFile _parameter_file;

  /// simulation wide parameters that require parameters from the parameter
  /// file

  /*! @brief Number of iterations of the ray tracing loop. */
  const uint_fast32_t _number_of_iterations;

  /*! @brief Number of photons used in every iteration of the ray tracing
   *  loop. */
  const uint_fast64_t _number_of_photons;

  /*! @brief Number of photons used during the first iteration of the ray
   *  tracing loop (this number is usually set to a lower value because the
   *  initial neutral fractions are set to very low values, making the grid very
   *  transparent). */
  const uint_fast64_t _number_of_photons_init;

  /// objects owned by the simulation that require parameters
  /// these have to be declared and initialized after the parameter file has
  /// been read

  /// non pointer objects owned by the simulation

  /*! @brief Abundances. */
  const Abundances _abundances;

  /// pointer objects owned by the simulation. These have to be deleted in the
  /// destructor.

  /*! @brief Cross sections for photoionization. */
  CrossSections *_cross_sections;

  /*! @brief Recombination rates. */
  RecombinationRates *_recombination_rates;

  /*! @brief Density function that sets the density field. */
  DensityFunction *_density_function;

  /*! @brief Mask that is added to the density field after the density function
   *  has been applied. */
  DensityMask *_density_mask;

  /*! @brief Grid used for the photoionization simulation. */
  DensityGrid *_density_grid;

  /*! @brief Distribution of discrete UV light sources. */
  PhotonSourceDistribution *_photon_source_distribution;

  /*! @brief Spectrum for the discrete UV sources. */
  PhotonSourceSpectrum *_photon_source_spectrum;

  /*! @brief Continuous source of UV light. */
  ContinuousPhotonSource *_continuous_photon_source;

  /*! @brief Spectrum for the continuous UV source. */
  PhotonSourceSpectrum *_continuous_photon_source_spectrum;

  /*! @brief Photon source that is used to actually get UV photons. */
  PhotonSource *_photon_source;

  /*! @brief Object used to write snapshot output. */
  DensityGridWriter *_density_grid_writer;

  /*! @brief Object used to compute the combined ionization and temperature
   *  balance at the end of a ray tracing step. */
  TemperatureCalculator *_temperature_calculator;

  /*! @brief Object used to do the ray tracing of photons through the grid in
   *  parallel. */
  IonizationPhotonShootJobMarket *_ionization_photon_shoot_job_market;

  /// internal timers

  /*! @brief Timer to quantify time spent in ray tracing. */
  Timer _work_timer;

  /*! @brief Timer to quantify time spent in serial parts of the code .*/
  Timer _serial_timer;

  /*! @brief Timer to quantify time spent in parallel parts of the code .*/
  Timer _parallel_timer;

  /*! @brief Timer to quantity total time spent in the ionization code .*/
  Timer _total_timer;

public:
  IonizationSimulation(const bool write_output,
                       const bool every_iteration_output,
                       const bool output_statistics,
                       const int_fast32_t num_thread,
                       const std::string parameterfile,
                       MPICommunicator *mpi_communicator = nullptr,
                       Log *log = nullptr);

  void initialize(DensityFunction *density_function = nullptr);
  void run(DensityGridWriter *density_grid_writer = nullptr);

  ~IonizationSimulation();
};

#endif // IONIZATIONSIMULATION_HPP
