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
 * @file CaproniPhotonSourceDistribution.hpp
 *
 * @brief Dwarf galaxy PhotonSourceDistribution based on the SN rates in
 * Caproni et al. (2017).
 *
 * The number of UV luminous O and B stars is determined by assuming that they
 * are a subset of the Caproni et al. (2017) SNII events (assuming a Chabrier/
 * Kroupa @f$M^{-2.3}@f$ high-mass tail for the IMF). The life times of the
 * sources are determined randomly from the high-mass PARSEC stellar evolution
 * tracks (Tang et al., 2014), while their luminosities are based on the
 * Sternberg et al. (2013) stellar models. Positions are based on a rough
 * polynomial fit to the Caproni et al. (2017) SN event radii.
 *
 * To allow more numerical flexibility, we introduced scale parameters that can
 * be used to change the normalisation of the OB star number function and the
 * stellar UV luminosity.
 *
 * For more details, see Vandenbroucke et al., in prep.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CAPRONIPHOTONSOURCEDISTRIBUTION_HPP
#define CAPRONIPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <vector>

/**
 * @brief Dwarf galaxy PhotonSourceDistribution based on the SN rates in
 * Caproni et al. (2017).
 */
class CaproniPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Normalisation factor for the OB number function. */
  const double _number_function_norm;

  /*! @brief Normalisation factor for the UV luminosity function. */
  const double _UV_luminosity_norm;

  /*! @brief Pseudo-random number generator. */
  RandomGenerator _random_generator;

  /*! @brief Positions of the sources (in m). */
  std::vector< CoordinateVector<> > _source_positions;

  /*! @brief Remaining lifetime of the sources (in s). */
  std::vector< double > _source_lifetimes;

  /*! @brief Luminosity of the sources (in s^-1). */
  std::vector< double > _source_luminosities;

  /*! @brief Total luminosity of all sources (in s^-1). */
  double _total_source_luminosity;

  /*! @brief Output file for the sources (if applicable). */
  std::ofstream *_output_file;

  /*! @brief Update time interval (in s). */
  const double _update_interval;

  /*! @brief Number of updates since the start of the simulation. */
  uint_fast32_t _number_of_updates;

  /*! @brief Indices of the sources (if output is enabled). */
  std::vector< uint_fast32_t > _source_indices;

  /*! @brief Index of the next source to add (if output is enabled). */
  uint_fast32_t _next_index;

  /**
   * @brief Generate a new source position.
   *
   * @return New source position (in m).
   */
  inline CoordinateVector<> generate_source_position() {
    return CoordinateVector<>();
  }

public:
  /**
   * @brief Constructor.
   *
   * @param number_function_norm Normalisation factor for the OB number
   * function.
   * @param UV_luminosity_norm Normalisation factor for the UV luminosity
   * function.
   * @param seed Seed for the pseudo-random number generator.
   * @param update_interval Time interval in between successive source
   * distribution updates (in s).
   * @param starting_time Start time of the simulation. The distribution is
   * evolved forward in time to this point before it is used (in s).
   * @param output_sources Should the source positions be written to a file?
   */
  inline CaproniPhotonSourceDistribution(const double number_function_norm,
                                         const double UV_luminosity_norm,
                                         const int_fast32_t seed,
                                         const double update_interval,
                                         const double starting_time,
                                         bool output_sources = false)
      : _number_function_norm(number_function_norm),
        _UV_luminosity_norm(UV_luminosity_norm), _random_generator(seed),
        _output_file(nullptr), _update_interval(update_interval),
        _number_of_updates(1), _next_index(0) {

    // generate sources
    _total_source_luminosity = 0.;

    if (output_sources) {
      _output_file = new std::ofstream("Caproni_source_positions.txt");
      *_output_file << "#time (s)\tx (m)\ty (m)\tz "
                       "(m)\tevent\tindex\tluminosity (s^-1)\tlifetime (s)\n";
      for (uint_fast32_t i = 0; i < _source_positions.size(); ++i) {
        _source_indices.push_back(_next_index);
        ++_next_index;
        const CoordinateVector<> &pos = _source_positions[i];
        *_output_file << 0. << "\t" << pos.x() << "\t" << pos.y() << "\t"
                      << pos.z() << "\t1\t" << _source_indices[i] << "\t"
                      << _source_luminosities[i] << "\t" << _source_lifetimes[i]
                      << "\n";
      }
      _output_file->flush();
    }

    // make sure the distribution is evolved up to the right starting time
    while (_number_of_updates * _update_interval <= starting_time) {

      const double total_time = _number_of_updates * _update_interval;
      // first clear out sources that do no longer exist
      size_t i = 0;
      while (i < _source_lifetimes.size()) {
        _source_lifetimes[i] -= _update_interval;
        if (_source_lifetimes[i] <= 0.) {
          // remove the element
          if (_output_file != nullptr) {
            *_output_file << total_time << "\t0.\t0.\t0.\t2\t"
                          << _source_indices[i] << "\t0\t0\n";
            _source_indices.erase(_source_indices.begin() + i);
          }

          _source_positions.erase(_source_positions.begin() + i);
          _source_lifetimes.erase(_source_lifetimes.begin() + i);
          _source_luminosities.erase(_source_luminosities.begin() + i);
        } else {
          // check the next element
          ++i;
        }
      }

      // now check if new sources need to be generated

      // update total source luminosity

      if (_output_file != nullptr) {
        _output_file->flush();
      }

      ++_number_of_updates;
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - number function norm: Normalisation factor for the OB number function
   *    (default: 1.)
   *  - UV luminosity norm: Normalisation factor for the UV luminosity function
   *    (default: 1.)
   *  - random seed: Random seed used to initialize the random generator that
   *    is used to sample the individual positions (default: 42)
   *  - update interval: Time interval in between successive distribution
   *    updates (default: 0.01 Gyr)
   *  - starting time: Starting time of the simulation. The distribution is
   *    evolved forward in time to this point before it is used
   *    (default: 0. Gyr)
   *  - output sources: Whether or not to write the source positions to a file
   *    (default: false)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  CaproniPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : CaproniPhotonSourceDistribution(
            params.get_value< double >(
                "PhotonSourceDistribution:number function norm", 1.),
            params.get_value< double >(
                "PhotonSourceDistribution:UV luminosity norm", 1.),
            params.get_value< int_fast32_t >(
                "PhotonSourceDistribution:random seed", 42),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:update interval", "0.01 Gyr"),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:starting time", "0. Gyr"),
            params.get_value< bool >("PhotonSourceDistribution:output sources",
                                     false)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~CaproniPhotonSourceDistribution() {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * The PhotonSourceDistribution will return exactly this number of valid
   * and unique positions by successive application of operator().
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_sources() const {
    return _source_positions.size();
  }

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {
    return _source_positions[index];
  }

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the photon source, used to determine how many photons are
   * emitted from this particular source.
   */
  virtual double get_weight(photonsourcenumber_t index) const {
    return _source_luminosities[index] / _total_source_luminosity;
  }

  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */
  virtual double get_total_luminosity() const {
    return _total_source_luminosity;
  }

  /**
   * @brief Update the distribution after the system moved to the given time.
   *
   * @param simulation_time Current simulation time (in s).
   * @return True if the distribution changed, false otherwise.
   */
  virtual bool update(const double simulation_time) {

    bool changed = false;
    while (_number_of_updates * _update_interval <= simulation_time) {

      const double total_time = _number_of_updates * _update_interval;
      // first clear out sources that do no longer exist
      size_t i = 0;
      while (i < _source_lifetimes.size()) {
        _source_lifetimes[i] -= _update_interval;
        if (_source_lifetimes[i] <= 0.) {
          // remove the element
          if (_output_file != nullptr) {
            *_output_file << total_time << "\t0.\t0.\t0.\t2\t"
                          << _source_indices[i] << "\t0\t0\n";
            _source_indices.erase(_source_indices.begin() + i);
          }

          _source_positions.erase(_source_positions.begin() + i);
          _source_lifetimes.erase(_source_lifetimes.begin() + i);
          _source_luminosities.erase(_source_luminosities.begin() + i);
          changed = true;
        } else {
          // check the next element
          ++i;
        }
      }

      // now check if new sources need to be generated
      // don't forget to update "changed variable"!
      // update total source luminosity!

      if (_output_file != nullptr) {
        _output_file->flush();
      }

      ++_number_of_updates;
    }

    return changed;
  }
};

#endif // CAPRONIPHOTONSOURCEDISTRIBUTION_HPP
