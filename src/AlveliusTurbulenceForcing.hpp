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
 * @brief Turbulence forcing using the method of Alvelius (1998).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 * @author Nina Sartorio (sartorio.nina@gmail.com)
 */
#ifndef ALVELIUSTURBULENCEFORCING_HPP
#define ALVELIUSTURBULENCEFORCING_HPP

#include "HydroDensitySubGrid.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief Photon package.
 */
class AlveliusTurbulenceForcing {
private:
  /*! @brief Real amplitudes. */
  std::vector< CoordinateVector<> > _amplitudes_real;

  /*! @brief Imaginary amplitudes. */
  std::vector< CoordinateVector<> > _amplitudes_imaginary;

  /**
   * @brief A table that has the number of the mode and the
   * x,y and z wavenumber components of that mode
   */
  std::vector< CoordinateVector< double > > _ktable;
  /**
   * @brief Direction unit vectors describing the
   * direction of the two force terms.
   */
  std::vector< CoordinateVector< double > > _e1;
  /**
   * @brief Direction unit vectors describing the
   * direction of the two force terms.
   */
  std::vector< CoordinateVector< double > > _e2;
  /**
   * @brief The forcing for each k mode between
   * kmin and kmax
   */
  std::vector< double > _kforce;

  /**
   * @brief Random Generator for turbulence
   * */
  RandomGenerator _RandomGenerator;

  /*! @brief Driving time step (in s). */
  double _time_step;

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
   * @param kmin Minimum wave number.
   * @param kmax Maximum wave number.
   * @param kforcing Wave number of highest forcing.
   * @param concentration_factor Width of the spectral function.
   * @param power_forcing Input power (probably has units).
   * @param seed Seed for the random generator.
   * @param dtfor Forcing time step (in s).
   */
  AlveliusTurbulenceForcing(double kmin, double kmax, double kforcing,
                            double concentration_factor, double power_forcing,
                            int_fast32_t seed, double dtfor)
      : _RandomGenerator(seed), _time_step(dtfor) {
    /**
     * @brief The force spectrum here prescribed is  Gaussian in shape:
     *F(k) = amplitude*exp^((k-kforcing)^2/concentration_factor)
     *       amplitude*gaussian_exponential
     */
    double spectra_sum = 0.;
    double cinv = 1. / concentration_factor;
    double inv2pi = 0.5 * (1 / M_PI);
    /**
     * @brief Iterate over all possible wavenumbers for
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
          double kk = std::sqrt(pwrk1 + pwrk2 + pwrk3);
          if (kk <= kmax && kk >= kmin) {
            const double kdiff = (kk - kforcing);
            const double sqrtk12 = std::sqrt(pwrk1 + pwrk2);
            const double invkk = 1. / kk;
            if (sqrtk12 > 0.) {
              const double invsqrtk12 = 1. / sqrtk12;
              _e1.push_back(
                  CoordinateVector<>(k2 * invsqrtk12, -k1 * invsqrtk12, 0.));
              _e2.push_back(CoordinateVector<>(k1 * k3 * invsqrtk12 * invkk,
                                               k2 * k3 * invsqrtk12 * invkk,
                                               -sqrtk12 * invkk));
            } else {
              const double sqrtk13 = std::sqrt(pwrk1 + pwrk3);
              const double invsqrtk13 = 1. / sqrtk13;
              _e1.push_back(
                  CoordinateVector<>(-k3 * invsqrtk13, 0., k1 * invsqrtk13));
              _e2.push_back(CoordinateVector<>(k1 * k2 * invsqrtk13 * invkk,
                                               -sqrtk13 * invkk,
                                               k2 * k3 * invsqrtk13 * invkk));
            }
            _ktable.push_back(CoordinateVector<>(k1, k2, k3));
            const double gaussian_spectra = std::exp(-kdiff * kdiff * cinv);
            spectra_sum += gaussian_spectra;
            _kforce.push_back(gaussian_spectra * inv2pi * invkk);
          }
        }
      }
    }
    /**
 @brief Initialize _amplitudes_imaginary _amplitudes_real
  */
    _amplitudes_real.resize(_ktable.size());
    _amplitudes_imaginary.resize(_ktable.size());
    /**
     * @brief Obtaining full expression for the forcing amplitude
     */
    double invSpectralSum = 1 / spectra_sum;
    double invdt = 1 / dtfor;
    for (uint_fast32_t i = 0; i < _kforce.size(); ++i) {
      _kforce[i] = _kforce[i] * power_forcing * invdt * invSpectralSum;
      _kforce[i] = std::sqrt(_kforce[i]);
    }
  }

  /**
   * @brief Update the turbulent amplitudes for the next time step.
   */
  inline void update_turbulence() {
    for (uint_fast32_t i = 0; i < _ktable.size(); ++i) {
      double RealRand[2];
      double ImRand[2];
      get_random_factors(_RandomGenerator, RealRand, ImRand);
      _amplitudes_real[i] =
          _kforce[i] * _e1[i] * RealRand[0] + _kforce[i] * _e2[i] * RealRand[1];
      _amplitudes_imaginary[i] =
          _kforce[i] * _e1[i] * ImRand[0] + _kforce[i] * _e2[i] * ImRand[1];
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

      cellit.get_hydro_variables().conserved(1) += _time_step * force.x();
      cellit.get_hydro_variables().conserved(2) += _time_step * force.y();
      cellit.get_hydro_variables().conserved(3) += _time_step * force.z();
    }
  }
};

#endif // ALVELIUSTURBULENCEFORCING_HPP
