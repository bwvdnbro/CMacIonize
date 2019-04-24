/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *                    Nina Sartorio (sartorio.nina@gmail.com)
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
 * @file AlveliusTurbulenceForcing.hpp
 *
 * @brief Turbulence forcing using the method of Alvelius (1999).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 * @author Nina Sartorio (sartorio.nina@gmail.com)
 */
#ifndef ALVELIUSTURBULENCEFORCING_HPP
#define ALVELIUSTURBULENCEFORCING_HPP

#include "HydroDensitySubGrid.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"

/**
 * @brief Turbulence forcing using the method of Alvelius (1999).
 */
class AlveliusTurbulenceForcing {
private:
  /*! @brief Real amplitudes of the forcing (in m s^-2). */
  std::vector< CoordinateVector<> > _amplitudes_real;

  /*! @brief Imaginary amplitudes of the forcing (in m s^-2). */
  std::vector< CoordinateVector<> > _amplitudes_imaginary;

  /*! @brief A table that has the x, y and z wavenumber components of the
   *  different driving modes (in m^-1). */
  std::vector< CoordinateVector<> > _ktable;

  /*! @brief Direction unit vectors describing the direction of the first force
   *  term for every mode. */
  std::vector< CoordinateVector<> > _e1;

  /*! @brief Direction unit vectors describing the direction of the second force
   *  term for every mode. */
  std::vector< CoordinateVector<> > _e2;

  /*! @brief The forcing for each mode (in m s^-2). */
  std::vector< double > _kforce;

  /*! @brief Random Generator for turbulence used to generate random forces. */
  RandomGenerator _random_generator;

  /*! @brief Driving time step (in s). */
  const double _time_step;

  /*! @brief Number of driving steps since the start of the simulation. */
  uint_fast32_t _number_of_driving_steps;

  /**
   * @brief Function gets the real and imaginary parts of the amplitudes
   * Aran and Bran of the unit vector e1 and e2, respectively, as in Eq. 11.
   * Alvelius (1999).
   *
   * @param RandGen Random number generator.
   * @param RealRand Output array to store the real parts of the amplitudes of
   * the forcing.
   * @param ImRand Output array to store the imaginary parts of the amplitudes
   * of the forcing.
   */
  static void get_random_factors(RandomGenerator &RandGen, double *RealRand,
                                 double *ImRand) {

    const double phi = 2. * M_PI * RandGen.get_uniform_random_double();
    const double ga = std::sin(phi);
    const double gb = std::cos(phi);
    const double theta1 = 2. * M_PI * RandGen.get_uniform_random_double();
    const double theta2 = 2. * M_PI * RandGen.get_uniform_random_double();
    RealRand[0] = std::cos(theta1) * ga;
    ImRand[0] = std::sin(theta1) * ga;
    RealRand[1] = std::cos(theta2) * gb;
    ImRand[1] = std::sin(theta2) * gb;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param kmin Minimum wave number (in m^-1).
   * @param kmax Maximum wave number (in m^-1).
   * @param kforcing Wave number of highest forcing (in m^-1).
   * @param concentration_factor Width of the spectral function (in m^-2).
   * @param power_forcing Input power (in m^2 s^-3).
   * @param seed Seed for the random generator.
   * @param dtfor Forcing time step (in s).
   */
  AlveliusTurbulenceForcing(double kmin, double kmax, double kforcing,
                            double concentration_factor, double power_forcing,
                            int_fast32_t seed, double dtfor)
      : _random_generator(seed), _time_step(dtfor),
        _number_of_driving_steps(0) {

    /* The force spectrum here prescribed is  Gaussian in shape:
     * F(k) = amplitude*exp^((k-kforcing)^2/concentration_factor)
     *        amplitude*gaussian_exponential
     */
    double spectra_sum = 0.;
    const double cinv = 1. / concentration_factor;
    // M_1_PI is equal to 1/pi
    const double inv2pi = 0.5 * M_1_PI;

    /*
     * Iterate over all possible wavenumbers for
     * x, y and z to obtain all the modes for the forcing.
     * Obtain all the spectra for te forcing as well as the direction
     * in which the forcing is going to be applied (given by unit vectors e1 and
     * e2)
     */
    for (double k1 = 0; k1 <= kmax; k1 += 1.) {
      for (double k2 = 0; k2 <= kmax; k2 += 1.) {
        for (double k3 = 0; k3 <= kmax; k3 += 1.) {
          const double pwrk1 = k1 * k1;
          const double pwrk2 = k2 * k2;
          const double pwrk3 = k3 * k3;
          const double kk = pwrk1 + pwrk2 + pwrk3;
          const double k = std::sqrt(kk);
          if (k <= kmax && k >= kmin) {
            const double kdiff = (k - kforcing);
            const double sqrtk12 = std::sqrt(pwrk1 + pwrk2);
            const double invkk = 1. / kk;
            const double invk = 1. / k;
            if (sqrtk12 > 0.) {
              const double invsqrtk12 = 1. / sqrtk12;
              _e1.push_back(
                  CoordinateVector<>(k2 * invsqrtk12, -k1 * invsqrtk12, 0.));
              _e2.push_back(CoordinateVector<>(k1 * k3 * invsqrtk12 * invk,
                                               k2 * k3 * invsqrtk12 * invk,
                                               -sqrtk12 * invk));
            } else {
              const double sqrtk13 = std::sqrt(pwrk1 + pwrk3);
              const double invsqrtk13 = 1. / sqrtk13;
              _e1.push_back(
                  CoordinateVector<>(-k3 * invsqrtk13, 0., k1 * invsqrtk13));
              _e2.push_back(CoordinateVector<>(k1 * k2 * invsqrtk13 * invk,
                                               -sqrtk13 * invk,
                                               k2 * k3 * invsqrtk13 * invk));
            }
            _ktable.push_back(CoordinateVector<>(k1, k2, k3));
            const double gaussian_spectra = std::exp(-kdiff * kdiff * cinv);
            spectra_sum += gaussian_spectra;
            _kforce.push_back(gaussian_spectra * inv2pi * invkk);
          }
        }
      }
    }

    // Initialize the amplitude vectors to the right size
    _amplitudes_real.resize(_ktable.size());
    _amplitudes_imaginary.resize(_ktable.size());

    // Obtain full expression for the forcing amplitude
    const double norm = power_forcing / (spectra_sum * dtfor);
    for (uint_fast32_t i = 0; i < _kforce.size(); ++i) {
      _kforce[i] *= norm;
      _kforce[i] = std::sqrt(_kforce[i]);
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  AlveliusTurbulenceForcing(ParameterFile &params)
      : AlveliusTurbulenceForcing(
            params.get_physical_value< QUANTITY_INVERSE_LENGTH >(
                "TurbulenceForcing:minimum wavenumber", "1. m^-1"),
            params.get_physical_value< QUANTITY_INVERSE_LENGTH >(
                "TurbulenceForcing:maximum wavenumber", "3. m^-1"),
            params.get_physical_value< QUANTITY_INVERSE_LENGTH >(
                "TurbulenceForing:peak forcing wavenumber", "2.5 m^-1"),
            params.get_physical_value< QUANTITY_INVERSE_SURFACE_AREA >(
                "TurbulenceForcing:concentration factor", "0.2 m^-2"),
            params.get_physical_value< QUANTITY_FORCING_POWER >(
                "TurbulenceForcing:power forcing", "1. m^2 s^-3"),
            params.get_value< int_fast32_t >("TurbulenceForcing:random seed",
                                             42),
            params.get_physical_value< QUANTITY_TIME >(
                "TurbulenceForcing:time step", "1.e-6 s")) {}

  /**
   * @brief Update the turbulent amplitudes for the next time step.
   *
   * @param end_of_timestep End of the current hydro time step (in s).
   * @return True if the amplitudes were updated and we need to add the
   * turbulent forcing to the cells in the grid.
   */
  inline bool update_turbulence(const double end_of_timestep) {

    if (_number_of_driving_steps * _time_step < end_of_timestep) {
      for (uint_fast32_t i = 0; i < _ktable.size(); ++i) {
        double RealRand[2];
        double ImRand[2];
        get_random_factors(_random_generator, RealRand, ImRand);
        _amplitudes_real[i] = _kforce[i] * _e1[i] * RealRand[0] +
                              _kforce[i] * _e2[i] * RealRand[1];
        _amplitudes_imaginary[i] =
            _kforce[i] * _e1[i] * ImRand[0] + _kforce[i] * _e2[i] * ImRand[1];
      }
      ++_number_of_driving_steps;
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief Add the turbulent forcing for the given subgrid.
   *
   * @param subgrid HydroDensitySubGrid to operate on.
   */
  inline void add_turbulent_forcing(HydroDensitySubGrid &subgrid) const {

    for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
         ++cellit) {

      const CoordinateVector<> x = cellit.get_cell_midpoint();
      CoordinateVector<> force;
      for (uint_fast32_t ik = 0; ik < _ktable.size(); ++ik) {
        const CoordinateVector<> k = _ktable[ik];
        const CoordinateVector<> fr = _amplitudes_real[ik];
        const CoordinateVector<> fi = _amplitudes_imaginary[ik];

        const double kdotx = CoordinateVector<>::dot_product(k, x);
        force += fr * std::cos(kdotx) - fi * std::sin(kdotx);
      }

      const double mdt =
          cellit.get_hydro_variables().get_conserved_mass() * _time_step;
      cellit.get_hydro_variables().conserved(1) += mdt * force.x();
      cellit.get_hydro_variables().conserved(2) += mdt * force.y();
      cellit.get_hydro_variables().conserved(3) += mdt * force.z();
    }
  }

  /**
   * @brief Dump the forcing object to the given restart file.
   *
   * @param restart_writer RestartWriter to write to.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    _random_generator.write_restart_file(restart_writer);
    restart_writer.write(_time_step);
    restart_writer.write(_number_of_driving_steps);

    const size_t number_of_modes = _ktable.size();
    restart_writer.write(number_of_modes);
    for (size_t i = 0; i < number_of_modes; ++i) {
      _ktable[i].write_restart_file(restart_writer);
      _e1[i].write_restart_file(restart_writer);
      _e2[i].write_restart_file(restart_writer);
      restart_writer.write(_kforce[i]);
    }
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline AlveliusTurbulenceForcing(RestartReader &restart_reader)
      : _random_generator(restart_reader),
        _time_step(restart_reader.read< double >()),
        _number_of_driving_steps(restart_reader.read< uint_fast32_t >()) {

    const size_t number_of_modes = restart_reader.read< size_t >();
    _amplitudes_real.resize(number_of_modes);
    _amplitudes_imaginary.resize(number_of_modes);
    _ktable.resize(number_of_modes);
    _e1.resize(number_of_modes);
    _e2.resize(number_of_modes);
    _kforce.resize(number_of_modes);
    for (size_t i = 0; i < number_of_modes; ++i) {
      _ktable[i] = CoordinateVector<>(restart_reader);
      _e1[i] = CoordinateVector<>(restart_reader);
      _e2[i] = CoordinateVector<>(restart_reader);
      _kforce[i] = restart_reader.read< double >();
    }
  }
};

#endif // ALVELIUSTURBULENCEFORCING_HPP
