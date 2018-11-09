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
 * Sternberg et al. (2003) stellar models. Positions are based on a rough
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
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <fstream>
#include <unistd.h>
#include <vector>

/**
 * @brief Dwarf galaxy PhotonSourceDistribution based on the SN rates in
 * Caproni et al. (2017).
 */
class CaproniPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /// factors that change the energetic output

  /*! @brief Normalisation factor for the stellar number function. */
  const double _number_function_norm;

  /*! @brief Normalisation factor for the UV luminosity function. */
  const double _UV_luminosity_norm;

  /*! @brief Boost factor for the feedback energy. */
  const double _boost_factor;

  /// IMF related functions

  /*! @brief @f$A@f$ factor used to generate stellar masses distributed
   *  according to the IMF. */
  const double _IMF_A;

  /*! @brief @f$A@f$ factor used to generate stellar masses distributed
   *  according to the IMF. */
  const double _IMF_B;

  /*! @brief @f$B@f$ factor used to generate stellar masses distributed
   *  according to the IMF. */
  const double _IMF_C;

  /*! @brief Lower mass limit for stars that are UV luminous (in Msol). */
  const double _OB_mass_limit_in_Msol;

  /// variables that control the random sampling

  /*! @brief Update time interval (in s). */
  const double _update_interval;

  /*! @brief Pseudo-random number generator. */
  RandomGenerator _random_generator;

  /*! @brief Number of updates since the start of the simulation. */
  uint_fast32_t _number_of_updates;

  /// stellar properties

  /*! @brief Positions of all stars (in m). */
  std::vector< CoordinateVector<> > _source_positions;

  /*! @brief Remaining lifetime of all stars (in s). */
  std::vector< double > _source_lifetimes;

  /*! @brief UV luminosity of all stars (in s^-1). */
  std::vector< double > _source_luminosities;

  /// UV source properties

  /*! @brief Indices of the O/B stars in the stellar vectors. */
  std::vector< uint_fast32_t > _OB_indices;

  /*! @brief Total luminosity of all sources (in s^-1). */
  double _total_source_luminosity;

  /*! @brief Flag that tells us whether the last population update changed
   *  anything to the O/B star population. */
  bool _Oflag;

  /// (optional) source output

  /*! @brief Output file for the sources (if applicable). */
  std::ofstream *_output_file;

  /*! @brief Indices of the sources (if output is enabled). */
  std::vector< uint_fast32_t > _source_indices;

  /*! @brief Index of the next source to add (if output is enabled). */
  uint_fast32_t _next_index;

public:
  /// semi-private functions that are not private to make them accessible from
  /// the unit test

  /**
   * @brief Get the expected number of stars for the given time.
   *
   * This number is based on the polynomial fit to the Caproni et al. (2017)
   * derived stellar number function.
   *
   * @param t Current simulation time (in s).
   * @return Expected number of massive stars at this time.
   */
  inline uint_fast32_t get_number_of_stars(const double t) const {

    // coefficients for the polynomial
    const double a[10] = {-4.0728750557e-145,  1.56644058448e-127,
                          -2.49066778113e-110, 2.1349428803e-93,
                          -1.07057674726e-76,  3.15595029622e-60,
                          -5.09080808197e-44,  3.44307404864e-28,
                          2.21557198304e-13,   431.31515864};

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
   * @brief Generate a random stellar mass based on the given uniform random
   * number.
   *
   * If the high-mass end of the IMF has the functional form
   * @f[
   *   \xi{}(M) = M^{\alpha{}},
   * @f]
   * with @f$\alpha{}@f$ a parameter, and if this functional form is valid for
   * the entire mass range @f$[M_L, M_U]@f$, then a random mass @f$M_r@f$
   * distributed according to the IMF is given by
   * @f[
   *   M_r = \left[u_r \left(M_U^{\alpha{} + 1} - M_L^{\alpha{} + 1}\right) +
   *               M_L^{\alpha{} + 1}\right]^{\frac{1}{\alpha{} + 1}},
   * @f]
   * with @f$u_r@f$ a uniform random number in the range @f$[0,1]@f$.
   *
   * We precompute the constants
   * @f[
   *   A = M_L^{\alpha{} + 1},
   * @f]
   * @f[
   *   B = M_U^{\alpha{} + 1} - A
   * @f]
   * and
   * @f[
   *   C = \frac{1}{\alpha{} + 1},
   * @f]
   * so that this equation reduces to
   * @f[
   *   M_r = \left(A + B u\right)^C.
   * @f]
   *
   * @param u Uniform random number in the range @f$[0, 1]@f$.
   * @return Random stellar mass distributed according to the IMF.
   */
  inline double get_random_stellar_mass(const double u) const {
    return std::pow(_IMF_A + _IMF_B * u, _IMF_C);
  }

  /**
   * @brief Get the stellar life time for a star with the given stellar mass
   * in solar masses.
   *
   * We use a double power-law fit to the stellar life time data from Tang et
   * al. (2014), the Z0.017Y0.279 model for masses in the range
   * @f$[20, 100]~{\rm{}M}_\odot{}@f$. The life time was computed by taking the
   * difference between the stellar ages at the start of phase 4 (near the ZAM)
   * and phase 8 (base of the RGB).
   *
   * @param M_in_Msol Mass of the star (in Msol).
   * @return Life time of the star (in s).
   */
  inline double get_stellar_lifetime(const double M_in_Msol) const {

    const double la[5] = {7.55609422e+13, 1.03371798e+16, -1.31168267e+00,
                          1.11162246e+18, -3.81030835e+00};
    return la[0] + la[1] * std::pow(M_in_Msol, la[2]) +
           la[3] * std::pow(M_in_Msol, la[4]);
  }

  /**
   * @brief Get the UV luminosity for a star with the given stellar mass in
   * solar masses.
   *
   * Only stars with a mass larger than the lower limit for UV luminous stars
   * are considered. We use a third order polynomial fit (and extrapolation) to
   * the OB luminosity data from Sternberg et al. (2003) (OB_LCV.dat).
   *
   * @param M_in_Msol Stellar mass (in Msol).
   * @return Corresponding UV luminosity (in s^-1).
   */
  inline double get_stellar_UV_luminosity(const double M_in_Msol) const {

    if (M_in_Msol < _OB_mass_limit_in_Msol) {
      return 0.;
    } else {
      const double a[4] = {-8.85154170718e+43, 2.21555601476e+46,
                           -4.25455875963e+47, 8.55819263554e+47};
      double luminosity = a[0] * M_in_Msol + a[1];
      luminosity = luminosity * M_in_Msol + a[2];
      luminosity = luminosity * M_in_Msol + a[3];
      return luminosity;
    }
  }

  /**
   * @brief Get a random galactic radius corresponding to the given time.
   *
   * The radius is based on a polynomial fit to the Caproni et al. (2017) SN
   * location data.
   *
   * @param t Current simulation time (in s).
   * @param random_generator RandomGenerator to use.
   * @return Random galactic radius (in m).
   */
  inline double get_galactic_radius(const double t,
                                    RandomGenerator &random_generator) const {

    // fit coefficients;
    const double a[10] = {-2.6765175763e-128, 7.3980382167e-111,
                          -8.44806044068e-94, 5.14530890285e-77,
                          -1.79963859431e-60, 3.62542751726e-44,
                          -4.00468400669e-28, 2.14263386338e-12,
                          -2968.65760812,     5.693802974e+18};

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
                  std::log(random_generator.get_uniform_random_double())) *
        std::cos(2. * M_PI * random_generator.get_uniform_random_double());

    return ravg + 3.086e18 * gauss;
  }

  /**
   * @brief Generate a new source position.
   *
   * @param t Current simulation time (in s).
   * @param random_generator RandomGenerator to use.
   * @return New source position (in m).
   */
  inline CoordinateVector<>
  generate_source_position(const double t,
                           RandomGenerator &random_generator) const {

    // get a random radius for this time
    const double r = get_galactic_radius(t, random_generator);

    // get a random direction
    const double cost = 2. * random_generator.get_uniform_random_double() - 1.;
    const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
    const double phi = 2. * M_PI * random_generator.get_uniform_random_double();
    const double cosp = std::cos(phi);
    const double sinp = std::sin(phi);

    return CoordinateVector<>(r * sint * cosp, r * sint * sinp, r * cost);
  }

  /**
   * @brief Evolve the stellar population to the given time.
   *
   * @param time Time to evolve to (in s).
   * @return std::vector<> containing all the positions were stars exploded.
   */
  inline std::vector< CoordinateVector<> >
  evolve_stellar_population(const double time) {

    std::vector< CoordinateVector<> > SN_positions;

    _Oflag = false;
    while (_number_of_updates * _update_interval <= time) {

      const double total_time = _number_of_updates * _update_interval;
      // first clear out sources that do no longer exist
      size_t i = 0;
      while (i < _source_lifetimes.size()) {
        _source_lifetimes[i] -= _update_interval;
        if (_source_lifetimes[i] <= 0.) {

          SN_positions.push_back(_source_positions[i]);
          if (_source_luminosities[i] > 0.) {
            _total_source_luminosity -= _source_luminosities[i];
            _Oflag = true;
          }

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
            generate_source_position(0., _random_generator);
        const double mass_in_Msol = get_random_stellar_mass(
            _random_generator.get_uniform_random_double());
        double luminosity = get_stellar_UV_luminosity(mass_in_Msol);
        double lifetime = get_stellar_lifetime(mass_in_Msol);
        // boost the luminosity
        luminosity *= _UV_luminosity_norm;
        // reduce the lifetime with a random fraction of the update interval
        // as the star could be born at any time during this interval
        lifetime -=
            _random_generator.get_uniform_random_double() * _update_interval;
        _source_positions.push_back(position);
        _source_lifetimes.push_back(lifetime);
        _source_luminosities.push_back(luminosity);
        if (luminosity > 0.) {
          _total_source_luminosity += luminosity;
          _Oflag = true;
        }

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

    return SN_positions;
  }

  /// member functions

  /**
   * @brief Constructor.
   *
   * @param number_function_norm Normalisation factor for the stellar number
   * function.
   * @param UV_luminosity_norm Normalisation factor for the UV luminosity
   * function.
   * @param SN_mass_limit Lower mass limit for stars that undergo SN explosions
   * (in kg).
   * @param OB_mass_limit Lower mass limit for stars that are UV luminous (in
   * kg).
   * @param stellar_mass_limit Upper limit for the mass of a star (in kg).
   * @param IMF_slope Power law slope for the high-mass end of the IMF.
   * @param seed Seed for the pseudo-random number generator.
   * @param update_interval Time interval in between successive source
   * distribution updates (in s).
   * @param starting_time Start time of the simulation. The distribution is
   * evolved forward in time to this point before it is used (in s).
   * @param boost_factor Boost factor for stellar feedback energy.
   * @param output_sources Should the source positions be written to a file?
   * @param log Log to write logging info to.
   */
  inline CaproniPhotonSourceDistribution(
      const double number_function_norm, const double UV_luminosity_norm,
      const double SN_mass_limit, const double OB_mass_limit,
      const double stellar_mass_limit, const double IMF_slope,
      const int_fast32_t seed, const double update_interval,
      const double starting_time, const double boost_factor,
      bool output_sources = false, Log *log = nullptr)
      : _number_function_norm(number_function_norm),
        _UV_luminosity_norm(UV_luminosity_norm), _boost_factor(boost_factor),
        _IMF_A(
            std::pow(SN_mass_limit / PhysicalConstants::get_physical_constant(
                                         PHYSICALCONSTANT_SOLAR_MASS),
                     IMF_slope + 1.)),
        _IMF_B(std::pow(stellar_mass_limit /
                            PhysicalConstants::get_physical_constant(
                                PHYSICALCONSTANT_SOLAR_MASS),
                        IMF_slope + 1.) -
               _IMF_A),
        _IMF_C(1. / (IMF_slope + 1.)),
        _OB_mass_limit_in_Msol(OB_mass_limit /
                               PhysicalConstants::get_physical_constant(
                                   PHYSICALCONSTANT_SOLAR_MASS)),
        _update_interval(std::min(update_interval, 9.9e13)),
        _random_generator(seed), _number_of_updates(1), _output_file(nullptr),
        _next_index(0) {

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
    const uint_fast32_t nexp =
        _number_function_norm *
        CaproniPhotonSourceDistribution::get_number_of_stars(0.);
    for (uint_fast32_t i = 0; i < nexp; ++i) {
      const CoordinateVector<> position =
          generate_source_position(0., _random_generator);
      const double mass_in_Msol = get_random_stellar_mass(
          _random_generator.get_uniform_random_double());
      double luminosity = get_stellar_UV_luminosity(mass_in_Msol);
      double lifetime = get_stellar_lifetime(mass_in_Msol);
      // boost the luminosity
      luminosity *= _UV_luminosity_norm;
      // reduce the lifetime with a random amount to account for the fact that
      // we sample the star at an arbitrary phase in its life
      lifetime *= _random_generator.get_uniform_random_double();
      // add the star
      _source_positions.push_back(position);
      _source_lifetimes.push_back(lifetime);
      _source_luminosities.push_back(luminosity);
      // check if the star is an O/B star
      if (luminosity > 0.) {
        // add it to the UV source list
        _total_source_luminosity += luminosity;
      }
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

    evolve_stellar_population(starting_time);

    // store the O/B star indices
    for (uint_fast32_t i = 0; i < _source_positions.size(); ++i) {
      if (_source_luminosities[i] > 0.) {
        _OB_indices.push_back(i);
      }
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - number function norm: Normalisation factor for the stellar number
   *    function (default: 1.)
   *  - UV luminosity norm: Normalisation factor for the UV luminosity function
   *    (default: 1.)
   *  - SN mass limit: Lower mass limit for stars that undergo SN explosions
   *    (default: 8. Msol)
   *  - OB mass limit: Lower mass limit for stars that are UV luminous (default:
   *    20. Msol)
   *  - stellar mass limit: Upper limit for the mass of a star (default: 100.
   *    Msol)
   *  - IMF slope: Power law slope for the high-mass end of the IMF (default:
   *    -2.3)
   *  - random seed: Random seed used to initialize the random generator that
   *    is used to sample the individual positions (default: 42)
   *  - update interval: Time interval in between successive distribution
   *    updates (default: 0.01 Gyr)
   *  - starting time: Starting time of the simulation. The distribution is
   *    evolved forward in time to this point before it is used
   *    (default: 0. Gyr)
   *  - boost factor: Boost factor for stellar feedback energy (default: 1.)
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
            params.get_physical_value< QUANTITY_MASS >(
                "PhotonSourceDistribution:SN mass limit", "8. Msol"),
            params.get_physical_value< QUANTITY_MASS >(
                "PhotonSourceDistribution:OB mass limit", "20. Msol"),
            params.get_physical_value< QUANTITY_MASS >(
                "PhotonSourceDistribution:stellar mass limit", "100. Msol"),
            params.get_value< double >("PhotonSourceDistribution:IMF slope",
                                       -2.3),
            params.get_value< int_fast32_t >(
                "PhotonSourceDistribution:random seed", 42),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:update interval", "0.01 Gyr"),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:starting time", "0. Gyr"),
            params.get_value< double >("PhotonSourceDistribution:boost factor",
                                       1.),
            params.get_value< bool >("PhotonSourceDistribution:output sources",
                                     false),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~CaproniPhotonSourceDistribution() {
    if (_output_file != nullptr) {
      delete _output_file;
    }
  }

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * The PhotonSourceDistribution will return exactly this number of valid
   * and unique positions by successive application of operator().
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_sources() const {
    return _OB_indices.size();
  }

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {
    return _source_positions[_OB_indices[index]];
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
    return _source_luminosities[_OB_indices[index]] / _total_source_luminosity;
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

    // if stellar feedback is active, the system will already be in the
    // correct state and this function will not do anything
    const bool old_Oflag = _Oflag;
    evolve_stellar_population(simulation_time);
    const bool new_Oflag = _Oflag;

    // clear the Oflag for the next time step
    _Oflag = false;

    if (old_Oflag | new_Oflag) {
      _OB_indices.clear();
      for (uint_fast32_t i = 0; i < _source_positions.size(); ++i) {
        if (_source_luminosities[i] > 0.) {
          _OB_indices.push_back(i);
        }
      }
    }

    // we have to check both the new and the old flag, as the population
    // may be evolved by either this function or the stellar feedback
    // function
    return old_Oflag | new_Oflag;
  }

  /**
   * @brief Add stellar feedback at the given time.
   *
   * @param grid DensityGrid to operate on.
   * @param time Current simulation time (in s).
   */
  virtual void add_stellar_feedback(DensityGrid &grid, const double time) {

    std::vector< CoordinateVector<> > SN_positions =
        evolve_stellar_population(time);

    for (uint_fast32_t i = 0; i < SN_positions.size(); ++i) {
      const CoordinateVector<> position = SN_positions[i];
      DensityGrid::iterator cell = grid.get_cell(position);
      cell.get_hydro_variables().set_energy_term(
          cell.get_hydro_variables().get_energy_term() + _boost_factor * 1.e44);
    }
  }

  /**
   * @brief Write the distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    restart_writer.write(_number_function_norm);
    restart_writer.write(_UV_luminosity_norm);
    restart_writer.write(_boost_factor);
    restart_writer.write(_IMF_A);
    restart_writer.write(_IMF_B);
    restart_writer.write(_IMF_C);
    restart_writer.write(_OB_mass_limit_in_Msol);
    restart_writer.write(_update_interval);
    _random_generator.write_restart_file(restart_writer);
    restart_writer.write(_number_of_updates);
    {
      const auto size = _source_positions.size();
      restart_writer.write(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i].write_restart_file(restart_writer);
      }
    }
    {
      const auto size = _source_lifetimes.size();
      restart_writer.write(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        restart_writer.write(_source_lifetimes[i]);
      }
    }
    {
      const auto size = _source_luminosities.size();
      restart_writer.write(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        restart_writer.write(_source_luminosities[i]);
      }
    }
    {
      const auto size = _OB_indices.size();
      restart_writer.write(size);
      for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
        restart_writer.write(_OB_indices[i]);
      }
    }
    restart_writer.write(_total_source_luminosity);
    restart_writer.write(_Oflag);
    const bool has_output = (_output_file != nullptr);
    restart_writer.write(has_output);
    if (has_output) {
      // store current position in the std::ofstream
      // we want to be able to continue writing from that point
      const auto filepos = _output_file->tellp();
      restart_writer.write(filepos);
      {
        const auto size = _source_indices.size();
        restart_writer.write(size);
        for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
          restart_writer.write(_source_indices[i]);
        }
      }
      restart_writer.write(_next_index);
    }
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline CaproniPhotonSourceDistribution(RestartReader &restart_reader)
      : _number_function_norm(restart_reader.read< double >()),
        _UV_luminosity_norm(restart_reader.read< double >()),
        _boost_factor(restart_reader.read< double >()),
        _IMF_A(restart_reader.read< double >()),
        _IMF_B(restart_reader.read< double >()),
        _IMF_C(restart_reader.read< double >()),
        _OB_mass_limit_in_Msol(restart_reader.read< double >()),
        _update_interval(restart_reader.read< double >()),
        _random_generator(restart_reader),
        _number_of_updates(restart_reader.read< uint_fast32_t >()) {

    {
      const std::vector< CoordinateVector<> >::size_type size =
          restart_reader.read< std::vector< CoordinateVector<> >::size_type >();
      _source_positions.resize(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i] = CoordinateVector<>(restart_reader);
      }
    }
    {
      const std::vector< double >::size_type size =
          restart_reader.read< std::vector< double >::size_type >();
      _source_lifetimes.resize(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        _source_lifetimes[i] = restart_reader.read< double >();
      }
    }
    {
      const std::vector< double >::size_type size =
          restart_reader.read< std::vector< double >::size_type >();
      _source_luminosities.resize(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        _source_luminosities[i] = restart_reader.read< double >();
      }
    }
    {
      const std::vector< uint_fast32_t >::size_type size =
          restart_reader.read< std::vector< uint_fast32_t >::size_type >();
      _OB_indices.resize(size);
      for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
        _OB_indices[i] = restart_reader.read< uint_fast32_t >();
      }
    }
    _total_source_luminosity = restart_reader.read< double >();
    _Oflag = restart_reader.read< bool >();
    const bool has_output = restart_reader.read< bool >();
    if (has_output) {
      const std::streampos filepos = restart_reader.read< std::streampos >();
      // truncate the original file to the size we were at
      if (truncate("Caproni_source_positions.txt", filepos) != 0) {
        cmac_error("Error while truncating output file!");
      }
      // now open the file in append mode
      _output_file =
          new std::ofstream("Caproni_source_positions.txt", std::ios_base::app);
      {
        const std::vector< uint_fast32_t >::size_type size =
            restart_reader.read< std::vector< uint_fast32_t >::size_type >();
        _source_indices.resize(size);
        for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
          _source_indices[i] = restart_reader.read< uint_fast32_t >();
        }
      }
      _next_index = restart_reader.read< uint_fast32_t >();
    }
  }
};

#endif // CAPRONIPHOTONSOURCEDISTRIBUTION_HPP
