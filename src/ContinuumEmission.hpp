/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ContinuumEmission.hpp
 *
 * @brief Routines to compute continuum emission from HII regions.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Routines to compute continuum emission from HII regions.
 *
 * The routines in this class are based on section 4.3 of Osterbrock & Ferland
 * (2006) (https://ui.adsabs.harvard.edu/abs/2006agna.book.....O). They use
 * expressions given by Brown & Mathews (1970)
 * (https://ui.adsabs.harvard.edu/abs/1970ApJ...160..939B) and the asymptotic
 * fits to the Kramers-Gaunt and free-free Gaunt factor from Seaton (1960)
 * (https://ui.adsabs.harvard.edu/abs/1960RPPh...23..313S).
 *
 * The results were tested against the values in Table 1 of Brown & Mathews
 * (1970).
 */
class ContinuumEmission {
private:
  /*! @brief Ionization energy of hydrogen (in J). */
  const double _IH;

  /*! @brief Constant factor in equation 2.2 of Seaton (1960) (in m^2). */
  const double _Seaton_factor;

  /*! @brief Constant factor in the first equation of Brown & Mathews (1970). */
  const double _gamma_n_factor;

public:
  /**
   * @brief Constructor.
   */
  ContinuumEmission()
      : _IH(13.6 * PhysicalConstants::get_physical_constant(
                       PHYSICALCONSTANT_ELECTRONVOLT)),
        _Seaton_factor(64. *
                       PhysicalConstants::get_physical_constant(
                           PHYSICALCONSTANT_FINE_STRUCTURE_CONSTANT) *
                       M_PI *
                       PhysicalConstants::get_physical_constant(
                           PHYSICALCONSTANT_BOHR_RADIUS) *
                       PhysicalConstants::get_physical_constant(
                           PHYSICALCONSTANT_BOHR_RADIUS) /
                       (3. * std::sqrt(3.))),
        _gamma_n_factor(
            std::sqrt(2. / M_PI) * 2. *
            PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK) /
            (PhysicalConstants::get_physical_constant(
                 PHYSICALCONSTANT_LIGHTSPEED) *
             PhysicalConstants::get_physical_constant(
                 PHYSICALCONSTANT_LIGHTSPEED) *
             PhysicalConstants::get_physical_constant(
                 PHYSICALCONSTANT_ELECTRON_MASS) *
             std::sqrt(PhysicalConstants::get_physical_constant(
                 PHYSICALCONSTANT_ELECTRON_MASS)))) {}

  /**
   * @brief Asymptotic fit to the Kramers-Gaunt factor from Seaton (1960),
   * equation 2.4.
   *
   * @param n Starting level of the electron.
   * @param epsilon Reduced photon energy, @f$\frac{h\nu{}}{Z^2I_{\rm{}H}} -
   * \frac{1}{n^2}@f$ (this is the relative difference between the photon energy
   * and the ground state of the atom).
   * @return Corresponding Kramers-Gaunt factor.
   */
  static inline double gauntII(const uint_fast32_t n, const double epsilon) {

    cmac_assert(n > 0);
    cmac_assert(epsilon >= 0.);

    const double u = n * n * epsilon;
    const double n_up1 = n * (u + 1.);
    const double factor = 1. / std::cbrt(n_up1 * n_up1);
    return 1. + 0.1728 * factor * (u - 1.) -
           0.0496 * factor * factor * (u * u + 4. * u / 3. + 1.);
  }

  /**
   * @brief Asymptotic fit to the velocity-averaged mean free-free Gaunt factor
   * from Seaton (1960), equation 2.38.
   *
   * @param Z Charge number of the ion.
   * @param hnu Energy of the photon (in J).
   * @param kT Energy of the electron gas (in J).
   * @return Corresponding velocity-averaged mean free-free Gaunt factor.
   */
  inline double gauntIII(const uint_fast32_t Z, const double hnu,
                         const double kT) const {

    cmac_assert(Z > 0);
    cmac_assert(hnu > 0.);
    cmac_assert(kT > 0.);

    const double u = kT / hnu;
    const double epsilon = hnu / (Z * Z * _IH);
    const double cbrtepsilon = std::cbrt(epsilon);
    return 1. + 0.1728 * cbrtepsilon * (1. + 2. * u) -
           0.0496 * cbrtepsilon * cbrtepsilon *
               (1. + 2. * u / 3. + 4. * u * u / 3.);
  }

  /**
   * @brief Get the photoionization cross section for photoionization of a
   * hydrogenic ion with charge @f$Z@f$ with an electron in the @f$n@f$th
   * excited state.
   *
   * This function computes equation 2.1 in Seaton (1960).
   *
   * @param Z Charge number of the ion.
   * @param n Starting level of the electron.
   * @param hnu Energy of the photon (in J).
   * @return Photoionization cross section for photoionization of the ion (in
   * m^2).
   */
  inline double level_photoionization_cross_section(const uint_fast32_t Z,
                                                    const uint_fast32_t n,
                                                    const double hnu) const {

    cmac_assert(Z > 0);
    cmac_assert(n > 0);
    cmac_assert(hnu >= 0.);

    const double epsilon = hnu / (Z * Z * _IH) - 1.0 / (n * n);
    return (epsilon >= 0.) ? _Seaton_factor * n /
                                 (Z * Z * std::cbrt(1. + n * n * epsilon)) *
                                 gauntII(n, epsilon)
                           : 0.;
  }

  /**
   * @brief Get the emission coefficient for continuous emission at photon
   * energy @f$h\nu{}@f$ due to recombinations to the @f$n@f$th excited state
   * of the hydrogenic ion with charge @f$Z@f$.
   *
   * This corresponds to the first equation in Brown & Mathews (1970).
   *
   * @param Z Charge number of the ion.
   * @param n Level to which the electron recombines.
   * @param hnu Photon energy (in J).
   * @param kT Energy of the electron gas (in J).
   * @return Corresponding continuous emission coefficient
   * (in J m^3 s^-1 Hz^-1).
   */
  inline double gamma_n(const uint_fast32_t Z, const uint_fast32_t n,
                        const double hnu, const double kT) const {

    cmac_assert(Z > 0);
    cmac_assert(n > 0);
    cmac_assert(hnu > 0.);
    cmac_assert(kT > 0.);

    const double n2 = n * n;
    const double In = Z * Z * _IH / n2;
    const double hnu_sqrtkT = hnu / std::sqrt(kT);
    return _gamma_n_factor * n2 * std::exp((In - hnu) / kT) * hnu_sqrtkT *
           hnu_sqrtkT * hnu_sqrtkT *
           level_photoionization_cross_section(Z, n, hnu);
  }
};
