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
 * @file DiscPatchPhotonSourceDistribution.hpp
 *
 * @brief Disc patch PhotonSourceDistribution.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISCPATCHPHOTONSOURCEDISTRIBUTION_HPP
#define DISCPATCHPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <vector>

/**
 * @brief Disc patch PhotonSourceDistribution.
 */
class DiscPatchPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Lifetime of a source (in s). */
  const double _source_lifetime;

  /*! @brief Ionising luminosity of a single source (in s^-1). */
  const double _source_luminosity;

  /*! @brief Probability of a new source being created (in s^-1). */
  const double _source_probability;

  /*! @brief Average number of sources at any given time. */
  const uint_fast32_t _average_number_of_sources;

  /*! @brief x component of the anchor of the rectangular disk (in m). */
  const double _anchor_x;

  /*! @brief y component of the anchor of the rectangular disk (in m). */
  const double _anchor_y;

  /*! @brief x side length of the rectangular disk (in m). */
  const double _sides_x;

  /*! @brief y side length of the rectangular disk (in m). */
  const double _sides_y;

  /*! @brief Origin of the Gaussian disk height distribution (in m). */
  const double _origin_z;

  /*! @brief Scale height of the Gaussian disk height distribution (in m). */
  const double _scaleheight_z;

  /*! @brief Pseudo-random number generator. */
  RandomGenerator _random_generator;

  /*! @brief Positions of the sources (in m). */
  std::vector< CoordinateVector<> > _source_positions;

  /*! @brief Remaining lifetime of the sources (in s). */
  std::vector< double > _source_lifetimes;

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
    const double x =
        _anchor_x + _random_generator.get_uniform_random_double() * _sides_x;
    const double y =
        _anchor_y + _random_generator.get_uniform_random_double() * _sides_y;
    // we use the Box-Muller method to sample the Gaussian
    const double z =
        _scaleheight_z *
            std::sqrt(-2. *
                      std::log(_random_generator.get_uniform_random_double())) *
            std::cos(2. * M_PI *
                     _random_generator.get_uniform_random_double()) +
        _origin_z;

    return CoordinateVector<>(x, y, z);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param source_lifetime Lifetime of a source (in s).
   * @param source_luminosity Ionising luminosity of a single source (in s^-1).
   * @param average_number Average number of sources at any given time.
   * @param anchor_x x component of the anchor of the rectangular disk (in m).
   * @param sides_x x side length of the rectangular disk (in m).
   * @param anchor_y  y component of the anchor of the rectangular disk (in m).
   * @param sides_y y side length of the rectangular disk (in m).
   * @param origin_z Origin of the Gaussian disk height distribution (in m).
   * @param scaleheight_z Scale height of the Gaussian disk height distribution
   * (in m).
   * @param seed Seed for the pseudo-random number generator.
   * @param update_interval Time interval in between successive source
   * distribution updates (in s).
   * @param output_sources Should the source positions be written to a file?
   */
  inline DiscPatchPhotonSourceDistribution(
      const double source_lifetime, const double source_luminosity,
      const uint_fast32_t average_number, const double anchor_x,
      const double sides_x, const double anchor_y, const double sides_y,
      const double origin_z, const double scaleheight_z,
      const int_fast32_t seed, const double update_interval,
      bool output_sources = false)
      : _source_lifetime(source_lifetime),
        _source_luminosity(source_luminosity),
        _source_probability(update_interval / source_lifetime),
        _average_number_of_sources(average_number), _anchor_x(anchor_x),
        _anchor_y(anchor_y), _sides_x(sides_x), _sides_y(sides_y),
        _origin_z(origin_z), _scaleheight_z(scaleheight_z),
        _random_generator(seed), _output_file(nullptr),
        _update_interval(update_interval), _number_of_updates(1),
        _next_index(0) {

    // generate sources
    for (uint_fast32_t i = 0; i < _average_number_of_sources; ++i) {
      const double lifetime =
          _random_generator.get_uniform_random_double() * _source_lifetime;
      _source_lifetimes.push_back(lifetime);
      _source_positions.push_back(generate_source_position());
    }

    if (output_sources) {
      _output_file = new std::ofstream("DiscPatch_source_positions.txt");
      *_output_file << "#time (s)\tx (m)\ty (m)\tz (m)\tevent\tindex\n";
      for (uint_fast32_t i = 0; i < _source_positions.size(); ++i) {
        _source_indices.push_back(_next_index);
        ++_next_index;
        const CoordinateVector<> &pos = _source_positions[i];
        *_output_file << 0. << "\t" << pos.x() << "\t" << pos.y() << "\t"
                      << pos.z() << "\t1\t" << _source_indices[i] << "\n";
      }
      _output_file->flush();
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - source lifetime: Lifetime of a source (default: 20. Myr)
   *  - source luminosity: Ionising luminosity of a single source
   *    (default: 1.e48 s^-1)
   *  - average number of sources: Average number of sources (default: 24)
   *  - anchor x: X position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides x: X side length of the 2D disc (default: 2. kpc)
   *  - anchor y: Y position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides y: Y side length of the 2D disc (default: 2. kpc)
   *  - origin z: Origin of the exponential disc profile in the z direction
   *    (default: 0. pc)
   *  - scaleheight z: Vertical scale height of the exponential disc profile
   *    (default: 63. pc)
   *  - random seed: Random seed used to initialize the random generator that
   *    is used to sample the individual positions (default: 42)
   *  - update interval: Time interval in between successive distribution
   *    updates (default: 0.1 Myr)
   *  - output sources: Whether or not to write the source positions to a file
   *    (default: false)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  DiscPatchPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : DiscPatchPhotonSourceDistribution(
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:source lifetime", "20. Myr"),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "PhotonSourceDistribution:source luminosity", "3.125e49 s^-1"),
            params.get_value< photonsourcenumber_t >(
                "PhotonSourceDistribution:average number of sources", 24),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:anchor x", "-1. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:sides x", "2. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:anchor y", "-1. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:sides y", "2. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:origin z", "0. pc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:scaleheight z", "63. pc"),
            params.get_value< int_fast32_t >(
                "PhotonSourceDistribution:random seed", 42),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:update interval", "0.1 Myr"),
            params.get_value< bool >("PhotonSourceDistribution:output sources",
                                     false)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~DiscPatchPhotonSourceDistribution() {}

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
    return 1. / get_number_of_sources();
  }

  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */
  virtual double get_total_luminosity() const {
    return _source_luminosity * get_number_of_sources();
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
                          << _source_indices[i] << "\n";
            _source_indices.erase(_source_indices.begin() + i);
          }

          _source_positions.erase(_source_positions.begin() + i);
          _source_lifetimes.erase(_source_lifetimes.begin() + i);
          changed = true;
        } else {
          // check the next element
          ++i;
        }
      }

      // now check if new sources need to be generated
      for (uint_fast32_t i = 0; i < _average_number_of_sources; ++i) {
        double x = _random_generator.get_uniform_random_double();
        if (x <= _source_probability) {
          // bingo: create a new source
          // the source could have been created at any given time during the
          // past
          // time step
          const double offset =
              _random_generator.get_uniform_random_double() * _update_interval;
          _source_lifetimes.push_back(_source_lifetime - offset);
          _source_positions.push_back(generate_source_position());
          if (_output_file != nullptr) {
            _source_indices.push_back(_next_index);
            ++_next_index;
            const CoordinateVector<> &pos = _source_positions.back();
            *_output_file << total_time << "\t" << pos.x() << "\t" << pos.y()
                          << "\t" << pos.z() << "\t1\t"
                          << _source_indices.back() << "\n";
          }
          changed = true;
        }
      }

      if (_output_file != nullptr) {
        _output_file->flush();
      }

      ++_number_of_updates;
    }

    return changed;
  }
};

#endif // DISCPATCHPHOTONSOURCEDISTRIBUTION_HPP
