/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file IonizationStateCalculator.cpp
 *
 * @brief IonizationStateCalculator implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "IonizationStateCalculator.hpp"
#include "DensityGrid.hpp"
#include "DensityValues.hpp"
#include "Error.hpp"
#include "RecombinationRates.hpp"
#include "UnitConverter.hpp"
#include <algorithm>
#include <cmath>

/**
 * @brief Constructor.
 *
 * @param luminosity Total ionizing luminosity of all photon sources (in s^-1).
 * @param helium_abundance Helium abundance.
 * @param recombination_rates RecombinationRates used in ionization balance
 * calculation.
 */
IonizationStateCalculator::IonizationStateCalculator(
    double luminosity, double helium_abundance,
    RecombinationRates &recombination_rates)
    : _luminosity(luminosity), _helium_abundance(helium_abundance),
      _recombination_rates(recombination_rates) {}

/**
 * @brief Solves the ionization and temperature equations based on the values of
 * the mean intensity integrals in each cell.
 *
 * @param nphoton Number of ionizing photons emitted by all sources.
 * @param grid DensityGrid for which the calculation is done.
 */
void IonizationStateCalculator::calculate_ionization_state(unsigned int nphoton,
                                                           DensityGrid &grid) {
  // Kenny's jfac contains a lot of unit conversion factors. These drop out
  // since we work in SI units.
  double jfac = _luminosity / nphoton;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    double cellvolume = it.get_volume();
    DensityValues &cell = it.get_values();
    double jfaccell = jfac / cellvolume;
    cell.set_old_neutral_fraction_H(cell.get_ionic_fraction(ELEMENT_H));
    double jH = jfaccell * cell.get_mean_intensity(ELEMENT_H);
    double jHe = jfaccell * cell.get_mean_intensity(ELEMENT_He);
    double ntot = cell.get_total_density();
    if (jH > 0. && ntot > 0.) {
      double T = cell.get_temperature();
      double alphaH = _recombination_rates.get_recombination_rate(ELEMENT_H, T);
      double alphaHe =
          _recombination_rates.get_recombination_rate(ELEMENT_He, T);
      // h0find
      double h0, he0;
      if (_helium_abundance) {
        find_H0(alphaH, alphaHe, jH, jHe, ntot, _helium_abundance, T, h0, he0);
      } else {
        find_H0_simple(alphaH, jH, ntot, T, h0);
        he0 = 0.;
      }

      cell.set_ionic_fraction(ELEMENT_H, h0);
      cell.set_ionic_fraction(ELEMENT_He, he0);

      // coolants
      double ne = ntot * (1. - h0 + _helium_abundance * (1. - he0));
      double t4 = T * 1.e-4;

      // carbon
      double C21 = jfac * cell.get_mean_intensity(ELEMENT_Cp1) / ne /
                   _recombination_rates.get_recombination_rate(ELEMENT_Cp1, T);
      double CTHerecom = 1.e-9 * 0.046 * t4 * t4; // units?
      double C32 =
          jfac * cell.get_mean_intensity(ELEMENT_Cp2) /
          (ne * _recombination_rates.get_recombination_rate(ELEMENT_Cp2, T)
           /*+ ntot * h0 * HCTRecom(4,6,T)*/
           + ntot * he0 * _helium_abundance * CTHerecom);
      double C31 = C32 * C21;
      double sumC = C21 + C31;
      cell.set_ionic_fraction(ELEMENT_Cp1, C21 / (1. + sumC));
      cell.set_ionic_fraction(ELEMENT_Cp2, C31 / (1. + sumC));
    } else {
      if (ntot > 0.) {
        cell.set_ionic_fraction(ELEMENT_H, 1.);
        cell.set_ionic_fraction(ELEMENT_He, 1.);
      } else {
        cell.set_ionic_fraction(ELEMENT_H, 0.);
        cell.set_ionic_fraction(ELEMENT_He, 0.);
      }
    }
  }
}

/**
 * @brief Iteratively find the neutral fractions of hydrogen and helium based on
 * the given value current values, the values of the intensity integrals and the
 * recombination rate.
 *
 * The equation for the ionization balance of hydrogen is
 * \f[
 *   n({\rm{}H}^0)\int_{\nu{}_i}^\infty{} \frac{4\pi{}J_\nu{}}{h\nu{}}
 *   a_\nu{}({\rm{}H}^0) {\rm{}d}\nu{} = n_e n({\rm{}H}^+)
 *   \alpha{}({\rm{}H}^0, T_e) - n_e P({\rm{}H}_{\rm{}OTS})n({\rm{}He}^+)
 *   \alpha{}_{2^1{\rm{}P}}^{\rm{}eff},
 * \f]
 * and that of helium
 * \f[
 *   n({\rm{}He}^0) \int_{\nu{}_i}^\infty{} \frac{4\pi{}J_\nu{}}{h\nu{}}
 *   a_\nu{}({\rm{}He}^0) {\rm{}d}\nu{} = n({\rm{}He}^0) n_e
 *   \alpha{}({\rm{}He}^0, T_e),
 * \f]
 * where the value of the integral on the left hand side of the equations, the
 * temperature \f$T_e\f$, the recombination rates \f$\alpha{}\f$, and the
 * current values of \f$n({\rm{}H}^0)\f$ and \f$n({\rm{}He}^0)\f$ (and hence all
 * other densities) are given. We want to determine the new values for the
 * densities.
 *
 * We start with helium. First, we derive the following expression for the
 * electron density \f$n_e\f$:
 * \f[
 *   n_e = (1-n({\rm{}H}^0)) n_{\rm{}H} + (1-n({\rm{}He}^0)) n_{\rm{}He},
 * \f]
 * which follows from charge conservation (assuming other elements do not
 * contribute free electrons due to their low abundances). Since the total
 * density we store in every cell is the hydrogen density, and abundances are
 * expressed relative to the hydrogen density, we can rewrite this as
 * \f[
 *   n_e = [(1-n({\rm{}H}^0)) + (1-n({\rm{}He}^0)) A_{\rm{}He}] n_{\rm{}tot}.
 * \f]
 * Using this, we can rewrite the helium ionization balance as
 * \f[
 *   n({\rm{}He}^0) = C_{\rm{}He} (1-n({\rm{}He}^0)) \frac{n_e}{n_{\rm{}tot}},
 * \f]
 * with \f$C_{\rm{}He} = \frac{\alpha{}({\rm{}He}^0, T_e) n_{\rm{}tot} }
 * {J_{\rm{}He}}\f$, a constant for a given temperature. Due to the \f$n_e\f$
 * factor, this equation is coupled to the ionization balance of hydrogen.
 * However, if we assume the hydrogen neutral fraction to be known, we can find
 * a closed expression for the helium neutral fraction:
 * \f[
 *   n({\rm{}He}^0) = \frac{-D - \sqrt{D^2 - 4A_{\rm{}He}C_{\rm{}He}B}}
 *   {2A_{\rm{}He}C_{\rm{}He}},
 * \f]
 * where
 * \f[
 *   D = -1 - 2A_{\rm{}He}C_{\rm{}He} - C_{\rm{}He} + C_{\rm{}He}n({\rm{}H}^0),
 * \f]
 * and
 * \f[
 *   B = C_{\rm{}He} + C_{\rm{}He}n({\rm{}H}^0) - A_{\rm{}He} C_{\rm{}He}.
 * \f]
 * This expression can be found by solving the quadratic equation in
 * \f$n({\rm{}He}^0)\f$. We choose the minus sign based on the fact that
 * \f$D < 0\f$ and the requirement that the helium neutral fraction be positive.
 *
 * To find the hydrogen neutral fraction, we have to address the fact that the
 * probability of a  Ly\f$\alpha{}\f$
 * photon being absorbed on the spot (\f$P({\rm{}H}_{\rm{}OTS})\f$) also depends
 * on the neutral fractions. We therefore use an iterative scheme to find the
 * hydrogen neutral fraction. The equation we want to solve is
 * \f[
 *   n({\rm{}H}^0) = C_{\rm{}H} (1-n({\rm{}H}^0)) \frac{n_e}{n_{\rm{}tot}},
 * \f]
 * which (conveniently) has the same form as the equation for helium. But know
 * we have
 * \f[
 *   C_{\rm{}H} = C_{{\rm{}H},1} + C_{{\rm{}H},2} P({\rm{}H}_{\rm{}OTS})
 *   \frac{1-n({\rm{}He}^0)}{1-n({\rm{}H}^0)},
 * \f]
 * with \f$C_{{\rm{}H},1} = \frac{\alpha{}({\rm{}H}^0, T_e) n_{\rm{}tot} }
 * {J_{\rm{}H}}\f$ and \f$C_{{\rm{}H},2} = \frac{A_{\rm{}He}
 * \alpha{}_{2^1{\rm{}P}}^{\rm{}eff} n_{\rm{}tot} }
 * {J_{\rm{}H}}\f$ two constants.
 *
 * For every iteration of the scheme, we calculate an approximate value for
 * \f$C_{\rm{}H}\f$, based on the neutral fractions obtained during the previous
 * iteration. We also calculate the helium neutral fraction, based on the
 * hydrogen neutral fraction from the previous iteration. Using these values,
 * we can then solve for the hydrogen neutral fraction. We repeat the process
 * until the relative difference between the obtained neutral fractions is
 * below some tolerance value.
 *
 * @param alphaH Hydrogen recombination rate (in m^3s^-1).
 * @param alphaHe Helium recombination rate (in m^3s^-1).
 * @param jH Hydrogen intensity integral (in s^-1).
 * @param jHe Helium intensity integral (in s^-1).
 * @param nH Hydrogen number density (in m^-3).
 * @param AHe Helium abundance \f$A_{\rm{}He}\f$ (relative w.r.t. hydrogen).
 * @param T Temperature (in K).
 * @param h0 Variable to store resulting hydrogen neutral fraction in.
 * @param he0 Variable to store resulting helium neutral fraction in.
 */
void IonizationStateCalculator::find_H0(double alphaH, double alphaHe,
                                        double jH, double jHe, double nH,
                                        double AHe, double T, double &h0,
                                        double &he0) {
  // unit conversions
  // we use Kenny's units inside this function, but all input and output units
  // are SI units
  alphaH = UnitConverter< QUANTITY_REACTION_RATE >::to_unit(alphaH, "cm^3s^-1");
  alphaHe =
      UnitConverter< QUANTITY_REACTION_RATE >::to_unit(alphaHe, "cm^3s^-1");
  jH = UnitConverter< QUANTITY_FREQUENCY >::to_unit(jH, "s^-1");
  jHe = UnitConverter< QUANTITY_FREQUENCY >::to_unit(jHe, "s^-1");
  nH = UnitConverter< QUANTITY_NUMBER_DENSITY >::to_unit(nH, "cm^-3");

  double alpha_e_2sP = 4.27e-14 * std::pow(T * 1.e-4, -0.695);
  double ch1 = alphaH * nH / jH;
  double ch2 = AHe * alpha_e_2sP * nH / jH;
  double che = 0.;
  if (jHe) {
    che = alphaHe * nH / jHe;
  }
  // initial guesses for the neutral fractions
  double h0old = 0.99 * (1. - std::exp(-0.5 / ch1));
  // by enforcing a relative difference of 10%, we make sure we have at least
  // one iteration
  h0 = 0.9 * h0old;
  double he0old = 1.;
  // we make sure che is 0 if the helium intensity integral is 0
  if (che) {
    he0old = 0.5 / che;
    // make sure the neutral fraction is at most 100%
    he0old = std::min(he0old, 1.);
  }
  // again, by using this value we make sure we have at least one iteration
  he0 = 0.;
  unsigned int niter = 0;
  while (std::abs(h0 - h0old) > 1.e-4 * h0old &&
         std::abs(he0 - he0old) > 1.e-4 * he0old) {
    ++niter;
    h0old = h0;
    if (he0 > 0.) {
      he0old = he0;
    } else {
      he0old = 0.;
    }
    // calculate a new guess for C_H
    double pHots = 1. / (1. + 77. * he0old / std::sqrt(T) / h0old);
    double ch = ch1 - ch2 * AHe * (1. - he0old) * pHots / (1. - h0old);

    // find the helium neutral fraction
    he0 = 1.;
    if (che) {
      double bhe = (1. + 2. * AHe - h0) * che + 1.;
      double t1he = 4. * AHe * (1. + AHe - h0) * che * che / bhe / bhe;
      if (t1he < 1.e-3) {
        // first order expansion of the square root in the exact solution of the
        // quadratic equation
        he0 = (1. + AHe - h0) * che / bhe;
      } else {
        // exact solution of the quadratic equation
        he0 = (bhe -
               std::sqrt(bhe * bhe - 4. * AHe * (1. + AHe - h0) * che * che)) /
              (2. * AHe * che);
      }
    }
    // find the hydrogen neutral fraction
    double b = ch * (2. + AHe - he0 * AHe) + 1.;
    double t1 = 4. * ch * ch * (1. + AHe - he0 * AHe) / b / b;
    if (t1 < 1.e-3) {
      h0 = ch * (1. + AHe - he0 * AHe) / b;
    } else {
      h0 = (b - std::sqrt(b * b - 4. * ch * ch * (1. + AHe - he0 * AHe))) /
           (2. * ch);
    }
    if (niter > 10) {
      // if we have a lot of iterations: use the mean value to speed up
      // convergence
      h0 = 0.5 * (h0 + h0old);
      he0 = 0.5 * (he0 + he0old);
    }
    if (niter > 20) {
      error("Too many iterations in ionization loop!");
    }
  }
}

/**
 * @brief find_H0() for a system without helium.
 *
 * We do not need to iterate in this case: the solution is simply given by the
 * solution of a quadratic equation.
 *
 * @param alphaH Hydrogen recombination rate (in m^3s^-1).
 * @param jH Hydrogen intensity integral (in s^-1).
 * @param nH Hydrogen number density (in m^-3).
 * @param T Temperature (in K).
 * @param h0 Variable to store resulting hydrogen neutral fraction in.
 */
void IonizationStateCalculator::find_H0_simple(double alphaH, double jH,
                                               double nH, double T,
                                               double &h0) {
  // unit conversions
  // we use Kenny's units inside this function, but all input and output units
  // are SI units
  alphaH = UnitConverter< QUANTITY_REACTION_RATE >::to_unit(alphaH, "cm^3s^-1");
  jH = UnitConverter< QUANTITY_FREQUENCY >::to_unit(jH, "s^-1");
  nH = UnitConverter< QUANTITY_NUMBER_DENSITY >::to_unit(nH, "cm^-3");

  if (jH > 0. && nH > 0.) {
    double aa = 0.5 * jH / nH / alphaH;
    double bb = 2. / aa;
    double cc = std::sqrt(bb + 1.);
    h0 = 1. + aa * (1. - cc);
  } else {
    h0 = 1.;
  }
}
