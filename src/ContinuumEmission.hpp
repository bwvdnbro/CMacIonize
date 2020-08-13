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

  /*! @brief Constant factor in the second equation of Brown & Mathews
   *  (1970). */
  const double _free_free_factor;

  /**
   * @brief Compute the constant prefactor in the Seaton (1960) cross section
   * expression.
   *
   * @return Constant prefactor.
   */
  static inline double compute_Seaton_factor() {

    const double alpha = PhysicalConstants::get_physical_constant(
        PHYSICALCONSTANT_FINE_STRUCTURE_CONSTANT);
    const double a0 =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOHR_RADIUS);

    return 64. * alpha * M_PI * a0 * a0 / (3. * std::sqrt(3.));
  }

  /**
   * @brief Compute the constant prefactor in the expression for bound-free
   * emission.
   *
   * @return Constant prefactor.
   */
  static inline double compute_gamma_n_factor() {

    const double h =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
    const double c =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED);
    const double m_e = PhysicalConstants::get_physical_constant(
        PHYSICALCONSTANT_ELECTRON_MASS);
    return std::sqrt(2. / M_PI) * 2. * h / (c * c * m_e * std::sqrt(m_e));
  }

  /**
   * @brief Compute the constant prefactor in the expression for free-free
   * emission.
   *
   * @param IH Ionization energy of hydrogen (in J).
   * @return Constant prefactor.
   */
  static inline double compute_free_free_factor(const double IH) {

    // note that the Brown & Mathews (1970) and Osterbrock & Ferland (2006)
    // equations use Gaussian CGS units. In these units, the electric
    // constant epsilon0 is set to 1, and a 4*pi factor is missing from the
    // Maxwell equations. As a result, the electron charge e has units of
    // statC, with 1 statC = 1 g^(1/2) cm^(3/2) s^-1
    // We convert the electron charge in SI units (C = A s) to the equivalent
    // "Gaussian SI" units kg^(1/2) m^(3/2) s^-1 and use the same expression
    // as in the works mentioned above.

    // electron charge in SI units: A s
    const double eSI = PhysicalConstants::get_physical_constant(
        PHYSICALCONSTANT_ELECTRON_CHARGE);
    // electron charge in "Gaussian SI" units: kg^(1/2) m^(3/2) s^-1
    const double e = eSI / std::sqrt(4. * M_PI *
                                     PhysicalConstants::get_physical_constant(
                                         PHYSICALCONSTANT_ELECTRIC_CONSTANT));
    const double e2 = e * e;
    const double m_e = PhysicalConstants::get_physical_constant(
        PHYSICALCONSTANT_ELECTRON_MASS);
    const double c =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED);
    const double h =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);

    return 32. * e2 * e2 * h / (3. * m_e * m_e * c * c * c) *
           std::sqrt(M_PI * IH / 3.);
  }

public:
  /**
   * @brief Constructor.
   */
  ContinuumEmission()
      : _IH(13.6 * PhysicalConstants::get_physical_constant(
                       PHYSICALCONSTANT_ELECTRONVOLT)),
        _Seaton_factor(compute_Seaton_factor()),
        _gamma_n_factor(compute_gamma_n_factor()),
        _free_free_factor(compute_free_free_factor(_IH)) {}

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
    const double nnepsp1 = 1. + n * n * epsilon;
    return (epsilon >= 0.)
               ? _Seaton_factor * n / (Z * Z * nnepsp1 * nnepsp1 * nnepsp1) *
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

  /**
   * @brief Get the total emission coefficient for continuous emission at photon
   * energy @f$h\nu{}@f$ due to recombination to all excited states of the
   * hydrogenic ion with charge @f$Z@f$.
   *
   * This is simply a sum over gamma_n() values for a sufficient number of
   * levels.
   *
   * @param Z Charge number of the ion.
   * @param hnu Photon energy (in J).
   * @param kT Energy of the electron gas (in J).
   * @return Corresponding total bound-free continuous emission coefficient
   * (in J m^3 s^-1 Hz^-1).
   */
  inline double gamma_bound_free(const uint_fast32_t Z, const double hnu,
                                 const double kT) const {

    cmac_assert(Z > 0);
    cmac_assert(hnu > 0.);
    cmac_assert(kT > 0.);

    double result = 0.;
    for (uint_fast32_t n = 1; n < 1000u; ++n) {
      result += gamma_n(Z, n, hnu, kT);
    }
    return result;
  }

  /**
   * @brief Get the total emission coefficient for continuous emission at photon
   * energy @f$h\nu{}@f$ due to free-free Brehmsstrahlung.
   *
   * This corresponds to the second equation in Brown & Mathews (1970).
   *
   * @param Z Charge number of the ion.
   * @param hnu Photon energy (in J).
   * @param kT Energy of the electron gas (in J).
   * @return Corresponding total free-free continuous emission coefficient
   * (in J m^3 s^-1 Hz^-1).
   */
  inline double gamma_free_free(const uint_fast32_t Z, const double hnu,
                                const double kT) const {

    cmac_assert(Z > 0);
    cmac_assert(hnu > 0.);
    cmac_assert(kT > 0.);

    return _free_free_factor * Z * Z / std::sqrt(kT) * std::exp(-hnu / kT) *
           gauntIII(Z, hnu, kT);
  }

  /**
   * @brief Get the continuous emission coefficient for emission by neutral
   * hydrogen at the given wavelength and the given electron gas temperature.
   *
   * This is simply the sum of gamma_bound_free() and gamma_free_free() for
   * @f$Z=1@f$.
   *
   * @param lambda Wavelength of emission (in m).
   * @param T Temperature of the electron gas (in K).
   * @return Continuous emission coefficient (in J m^3 s^-1 Hz^-1).
   */
  inline double gamma_HI(const double lambda, const double T) const {

    cmac_assert(lambda > 0.);
    cmac_assert(T > 0.);

    const double nu =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED) /
        lambda;
    const double hnu =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK) * nu;
    const double kT =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN) *
        T;

    return gamma_bound_free(1, hnu, kT) + gamma_free_free(1, hnu, kT);
  }
};
