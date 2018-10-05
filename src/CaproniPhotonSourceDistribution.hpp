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
#include <cmath>
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

public:
  /// private functions that are public to make them accessible from the unit
  /// test

  /**
   * @brief Get the UV luminosity and life time for a random OB star distributed
   * according to the IMF.
   *
   * @param luminosity Output UV luminosity of the random star (in s^-1).
   * @param lifetime Output life time of the star (in s).
   * @return Mass of the star (in Msol; only for test purposes).
   */
  inline double get_random_star(double &luminosity, double &lifetime) {

    const double alphap1 = -1.3;
    const double alphap1inv = 1. / alphap1;
    const double mlowterm = std::pow(20., alphap1);
    const double mrangefac = std::pow(100., alphap1) - mlowterm;

    // get a random mass distributed according to the IMF
    const double u = _random_generator.get_uniform_random_double();
    const double mstar = std::pow(u * mrangefac + mlowterm, alphap1inv);

    // get the corresponding luminosity
    // we use a 3th order polynomial fit (and extrapolation) for the OB
    // luminosity data from Sternberg et al. (2013) (OB_LCV.dat)
    const double a[4] = {-8.85154170718e+43, 2.21555601476e+46,
                         -4.25455875963e+47, 8.55819263554e+47};
    luminosity = a[0] * mstar + a[1];
    luminosity = luminosity * mstar + a[2];
    luminosity = luminosity * mstar + a[3];

    // get the corresponding life time
    // we use a 10th order polynomial fit to the stellar life time data from
    // Tang et al. (2014), the Z0.017Y0.279 model for masses in the range
    // [20, 100] Msol. The life time was computed by taking the difference
    // between the stellar ages at the start of phase 4 (near the ZAM) and
    // phase 8 (base of the RGB).
    const double la[11] = {
        0.00645837785727,   -3.99457743218,     1089.30492119,
        -172343.582909,     17515331.6211,      -1195587654.33,
        55641063153.0,      -1.75287305075e+12, 3.61957288273e+13,
        -4.54128336105e+14, 2.91526506459e+15};
    lifetime = la[0] * mstar + la[1];
    lifetime = lifetime * mstar + la[2];
    lifetime = lifetime * mstar + la[3];
    lifetime = lifetime * mstar + la[4];
    lifetime = lifetime * mstar + la[5];
    lifetime = lifetime * mstar + la[6];
    lifetime = lifetime * mstar + la[7];
    lifetime = lifetime * mstar + la[8];
    lifetime = lifetime * mstar + la[9];
    lifetime = lifetime * mstar + la[10];

    return mstar;
  }

  /**
   * @brief Get the expected number of stars for the given time.
   *
   * This number is based on the polynomial fit to the Caproni et al. (2017)
   * derived OB number function.
   *
   * @param t Current simulation time (in s).
   * @return Expected number of OB stars at this time.
   */
  inline static uint_fast32_t get_number_of_stars(const double t) {

    // coefficients for the polynomial
    const double a[10] = {-4.02787220841e-146, 1.60421512996e-128,
                          -2.61501324962e-111, 2.28378552108e-94,
                          -1.16252321273e-77,  3.47502087069e-61,
                          -5.70743788013e-45,  4.0416648576e-29,
                          1.81489889811e-15,   28.7796833735};

    // evaluate polynomial in O(n) time using Horner's method
    double result = a[0] * t + a[1];
    result = result * t + a[2];
    result = result * t + a[3];
    result = result * t + a[4];
    result = result * t + a[5];
    result = result * t + a[6];
    result = result * t + a[7];
    result = result * t + a[8];
    result = result * t + a[9];

    return result;
  }

  /**
   * @brief Get a random galactic radius corresponding to the given time.
   *
   * The radius is based on a polynomial fit to the Caproni et al. (2017) SN
   * location data.
   *
   * @param t Current simulation time (in s).
   * @return Random galactic radius (in m).
   */
  inline double get_galactic_radius(const double t) {

    // fit coefficients;
    const double a[10] = {-2.47715326891e-128, 6.89912786829e-111,
                          -7.94169102884e-94,  4.87832127858e-77,
                          -1.72187755687e-60,  3.502399403e-44,
                          -3.90591687341e-28,  2.1013797052e-12,
                          -2757.1605573,       4.61751485207e+18};

    double ravg = a[0] * t + a[1];
    ravg = ravg * t + a[2];
    ravg = ravg * t + a[3];
    ravg = ravg * t + a[4];
    ravg = ravg * t + a[5];
    ravg = ravg * t + a[6];
    ravg = ravg * t + a[7];
    ravg = ravg * t + a[8];
    ravg = ravg * t + a[9];

    const double gauss =
        std::sqrt(-2. *
                  std::log(_random_generator.get_uniform_random_double())) *
        std::cos(2. * M_PI * _random_generator.get_uniform_random_double());

    return ravg + 3.086e18 * gauss;
  }

  /**
   * @brief Generate a new source position.
   *
   * @param t Current simulation time (in s).
   * @return New source position (in m).
   */
  inline CoordinateVector<> generate_source_position(const double t) {

    // get a random radius for this time
    const double r = get_galactic_radius(t);

    // get a random direction
    const double cost = 2. * _random_generator.get_uniform_random_double() - 1.;
    const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
    const double phi =
        2. * M_PI * _random_generator.get_uniform_random_double();
    const double cosp = std::cos(phi);
    const double sinp = std::sin(phi);

    return CoordinateVector<>(r * sint * cosp, r * sint * sinp, r * cost);
  }

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
   * @param log Log to write logging info to.
   */
  inline CaproniPhotonSourceDistribution(const double number_function_norm,
                                         const double UV_luminosity_norm,
                                         const int_fast32_t seed,
                                         const double update_interval,
                                         const double starting_time,
                                         bool output_sources = false,
                                         Log *log = nullptr)
      : _number_function_norm(number_function_norm),
        _UV_luminosity_norm(UV_luminosity_norm), _random_generator(seed),
        _output_file(nullptr),
        _update_interval(std::min(update_interval, 9.9e13)),
        _number_of_updates(1), _next_index(0) {

    if (log != nullptr && update_interval > 9.9e13) {
      log->write_warning("CaproniPhotonSourceDistribution update interval "
                         "larger than the minimum source life time. Reset to "
                         "smaller value!");
    }

    if (starting_time > 6.4e16) {
      cmac_error("Starting time larger than number function validity range!");
    }

    // generate sources
    _total_source_luminosity = 0.;
    const uint_fast32_t nexp = _number_function_norm * get_number_of_stars(0.);
    for (uint_fast32_t i = 0; i < nexp; ++i) {
      const CoordinateVector<> position = generate_source_position(0.);
      double luminosity, lifetime;
      get_random_star(luminosity, lifetime);
      luminosity *= _UV_luminosity_norm;
      // reduce the lifetime with a random amount
      lifetime *= _random_generator.get_uniform_random_double();
      _source_positions.push_back(position);
      _source_lifetimes.push_back(lifetime);
      _source_luminosities.push_back(luminosity);
      _total_source_luminosity += luminosity;
    }

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

          _total_source_luminosity -= _source_luminosities[i];

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
      const uint_fast32_t nexp =
          _number_function_norm * get_number_of_stars(total_time);
      for (uint_fast32_t i = _source_positions.size(); i < nexp; ++i) {
        const CoordinateVector<> position =
            generate_source_position(total_time);
        double luminosity, lifetime;
        get_random_star(luminosity, lifetime);
        luminosity *= _UV_luminosity_norm;
        // reduce the lifetime with a random fraction of the update interval
        lifetime -=
            _random_generator.get_uniform_random_double() * _update_interval;
        _source_positions.push_back(position);
        _source_lifetimes.push_back(lifetime);
        _source_luminosities.push_back(luminosity);
        _total_source_luminosity += luminosity;

        if (_output_file != nullptr) {
          _source_indices.push_back(_next_index);
          ++_next_index;
          const CoordinateVector<> &pos = _source_positions[i];
          *_output_file << 0. << "\t" << pos.x() << "\t" << pos.y() << "\t"
                        << pos.z() << "\t1\t" << _source_indices[i] << "\t"
                        << _source_luminosities[i] << "\t"
                        << _source_lifetimes[i] << "\n";
        }
      }

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
                                     false),
            log) {}

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

    if (simulation_time > 6.4e16) {
      cmac_error("Simulation time larger than number function validity range!");
    }

    bool changed = false;
    while (_number_of_updates * _update_interval <= simulation_time) {

      const double total_time = _number_of_updates * _update_interval;
      // first clear out sources that do no longer exist
      size_t i = 0;
      while (i < _source_lifetimes.size()) {
        _source_lifetimes[i] -= _update_interval;
        if (_source_lifetimes[i] <= 0.) {

          _total_source_luminosity -= _source_luminosities[i];

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
      const uint_fast32_t nexp =
          _number_function_norm * get_number_of_stars(total_time);
      for (uint_fast32_t i = _source_positions.size(); i < nexp; ++i) {
        const CoordinateVector<> position =
            generate_source_position(total_time);
        double luminosity, lifetime;
        get_random_star(luminosity, lifetime);
        luminosity *= _UV_luminosity_norm;
        // reduce the lifetime with a random fraction of the update interval
        lifetime -=
            _random_generator.get_uniform_random_double() * _update_interval;
        _source_positions.push_back(position);
        _source_lifetimes.push_back(lifetime);
        _source_luminosities.push_back(luminosity);
        _total_source_luminosity += luminosity;

        if (_output_file != nullptr) {
          _source_indices.push_back(_next_index);
          ++_next_index;
          const CoordinateVector<> &pos = _source_positions[i];
          *_output_file << total_time << "\t" << pos.x() << "\t" << pos.y()
                        << "\t" << pos.z() << "\t1\t" << _source_indices[i]
                        << "\t" << _source_luminosities[i] << "\t"
                        << _source_lifetimes[i] << "\n";
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

#endif // CAPRONIPHOTONSOURCEDISTRIBUTION_HPP
