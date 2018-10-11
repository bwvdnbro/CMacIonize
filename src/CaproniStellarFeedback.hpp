/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file CaproniStellarFeedback.hpp
 *
 * @brief Dwarf galaxy stellar feedback prescription based on the SN rates in
 * Caproni et al. (2017).
 *
 * For more details, see Vandenbroucke et al., in prep.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CAPRONISTELLARFEEDBACK_HPP
#define CAPRONISTELLARFEEDBACK_HPP

#include "CaproniStellarRoutines.hpp"
#include "DensityGrid.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"

#include <fstream>

/**
 * @brief Dwarf galaxy stellar feedback prescription based on the SN rates in
 * Caproni et al. (2017).
 */
class CaproniStellarFeedback {
private:
  /*! @brief RandomGenerator used to generate random positions. */
  RandomGenerator _random_generator;

  /*! @brief Time the last SN event occured (in s). */
  double _time_since_last;

  /*! @brief Output file to write log output to. */
  std::ofstream *_output_file;

public:
  /**
   * @brief Get the instanteneous supernova rate at the given time.
   *
   * We use a 8th order polynomial fit to the Caproni et al. (2017) data.
   *
   * @param time Current simulation time (in s).
   * @return Instanteneous supernova rate (in s^-1).
   */
  inline static double get_SN_rate(const double time) {

    if (time > 6.4e16) {
      cmac_error("Time value outside fit validity range!");
    }

    // fit coefficients
    const double a[9] = {
        6.89799700195e-143, -1.91715448814e-125, 2.19302246787e-108,
        -1.32201654133e-91, 4.43164838505e-75,   -7.80981517111e-59,
        5.39845077482e-43,  9.70796661139e-28,   -8.44606535214e-14};

    double rate = a[0] * time + a[1];
    rate = rate * time + a[2];
    rate = rate * time + a[3];
    rate = rate * time + a[4];
    rate = rate * time + a[5];
    rate = rate * time + a[6];
    rate = rate * time + a[7];
    rate = rate * time + a[8];

    // make sure the return value is positive
    return std::max(0., rate);
  }

  /**
   * @brief Constructor.
   *
   * @param seed Seed for the RandomGenerator.
   * @param write_output Output SN events to log file?
   * @param log Log to write logging info to.
   */
  inline CaproniStellarFeedback(const uint_fast32_t seed,
                                const bool write_output, Log *log = nullptr)
      : _random_generator(seed), _time_since_last(0.), _output_file(nullptr) {

    if (write_output) {
      _output_file = new std::ofstream("Caproni_SN_events.txt");
      *_output_file << "# time (s)\tx (m)\ty (m)\tz (m)\n";
      _output_file->flush();
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read:
   *  - seed: Seed for the random number generator (default: 42)
   *  - write output: Output SN events to a log file? (default: false)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline CaproniStellarFeedback(ParameterFile &params, Log *log = nullptr)
      : CaproniStellarFeedback(
            params.get_value< uint_fast32_t >("StellarFeedback:seed", 42),
            params.get_value< bool >("StellarFeedback:write output", false),
            log) {}

  inline ~CaproniStellarFeedback() {
    if (_output_file != nullptr) {
      delete _output_file;
    }
  }

  /**
   * @brief Add stellar feedback at the given time.
   *
   * @param grid DensityGrid to operate on.
   * @param time Current simulation time (in s).
   */
  inline void add_stellar_feedback(DensityGrid &grid, const double time) {

    // compute the number of SN events during the last time step
    const double dt = time - _time_since_last;
    const uint_fast32_t N_SN = std::round(get_SN_rate(time) * dt);
    for (uint_fast32_t i = 0; i < N_SN; ++i) {
      const CoordinateVector<> position =
          CaproniStellarRoutines::generate_source_position(time,
                                                           _random_generator);
      DensityGrid::iterator cell = grid.get_cell(position);
      cell.get_hydro_variables().set_energy_term(
          cell.get_hydro_variables().get_energy_term() + 1.e44);

      if (_output_file != nullptr) {
        *_output_file << time << "\t" << position.x() << "\t" << position.y()
                      << "\t" << position.z() << "\n";
      }
    }

    if (N_SN > 0) {
      _time_since_last = time;
    }

    if (_output_file != nullptr) {
      _output_file->flush();
    }
  }
};

#endif // CAPRONISTELLARFEEDBACK_HPP
