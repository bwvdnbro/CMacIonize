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
 * @file TemperatureCalculator.cpp
 *
 * @brief TemperatureCalculator implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "TemperatureCalculator.hpp"
#include "Abundances.hpp"
#include "ChargeTransferRates.hpp"
#include "DensityGrid.hpp"
#include "DensityGridTraversalJobMarket.hpp"
#include "DensityValues.hpp"
#include "IonizationStateCalculator.hpp"
#include "LineCoolingData.hpp"
#include "RecombinationRates.hpp"
#include "WorkDistributor.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * @param luminosity Total ionizing luminosity of all photon sources (in s^-1).
 * @param abundances Abundances.
 * @param pahfac PAH heating factor.
 * @param crfac Cosmic ray heating factor.
 * @param line_cooling_data LineCoolingData use to calculate cooling due to line
 * emission.
 * @param recombination_rates RecombinationRates used to calculate ionic
 * fractions.
 * @param charge_transfer_rates ChargeTransferRates used to calculate ionic
 * fractions.
 */
TemperatureCalculator::TemperatureCalculator(
    double luminosity, Abundances &abundances, double pahfac, double crfac,
    LineCoolingData &line_cooling_data, RecombinationRates &recombination_rates,
    ChargeTransferRates &charge_transfer_rates)
    : _luminosity(luminosity), _abundances(abundances), _pahfac(pahfac),
      _crfac(crfac), _line_cooling_data(line_cooling_data),
      _recombination_rates(recombination_rates),
      _charge_transfer_rates(charge_transfer_rates) {}

/**
 * @brief Function that calculates the cooling and heating rate for a given
 * cell.
 *
 * @param h0 Variable to store the hydrogen neutral fraction in.
 * @param he0 Variable to store the helium neutral fraction in.
 * @param gain Total energy gain due to heating.
 * @param loss Total energy loss due to cooling.
 * @param T Temperature (in K).
 * @param cell DensityValues of the cell.
 * @param jfac Normalization factor for the mean intensities.
 * @param abundances Abundances.
 * @param hfac Normalization factor for the heating integrals.
 * @param pahfac Normalization factor for PAH heating.
 * @param crfac Normalization factor for cosmic ray heating.
 * @param data LineCoolingData used to calculate line cooling.
 * @param rates RecombinationRates used to calculate ionic fractions.
 * @param ctr ChargeTransferRates used to calculate ionic fractions.
 */
void TemperatureCalculator::ioneng(
    double &h0, double &he0, double &gain, double &loss, double T,
    DensityGrid::iterator &cell, double jfac, Abundances &abundances,
    double hfac, double pahfac, double crfac, LineCoolingData &data,
    RecombinationRates &rates, ChargeTransferRates &ctr) {
  double alphaH = rates.get_recombination_rate(ION_H_n, T);
  double alphaHe = rates.get_recombination_rate(ION_He_n, T);
  double alphaC[2];
  alphaC[0] = rates.get_recombination_rate(ION_C_p1, T);
  alphaC[1] = rates.get_recombination_rate(ION_C_p2, T);
  double alphaN[3];
  alphaN[0] = rates.get_recombination_rate(ION_N_n, T);
  alphaN[1] = rates.get_recombination_rate(ION_N_p1, T);
  alphaN[2] = rates.get_recombination_rate(ION_N_p2, T);
  double alphaO[2];
  alphaO[0] = rates.get_recombination_rate(ION_O_n, T);
  alphaO[1] = rates.get_recombination_rate(ION_O_p1, T);
  double alphaNe[2];
  alphaNe[0] = rates.get_recombination_rate(ION_Ne_n, T);
  alphaNe[1] = rates.get_recombination_rate(ION_Ne_p1, T);
  double alphaS[3];
  alphaS[0] = rates.get_recombination_rate(ION_S_p1, T);
  alphaS[1] = rates.get_recombination_rate(ION_S_p2, T);
  alphaS[2] = rates.get_recombination_rate(ION_S_p3, T);

  double t4 = T * 1.e-4;
  double alpha_e_2sP = 4.27e-14 * std::pow(t4, -0.695);
  double n = cell.get_number_density();

  double jH = jfac * cell.get_mean_intensity(ION_H_n);
  double jHe = jfac * cell.get_mean_intensity(ION_He_n);
  double AHe = abundances.get_abundance(ELEMENT_He);
  IonizationStateCalculator::find_H0(alphaH, alphaHe, jH, jHe, n, AHe, T, h0,
                                     he0);

  double ne = n * (1. - h0 + AHe * (1. - he0));
  cmac_assert(ne == ne);
  double nhp = n * (1. - h0);
  double nhep = (1. - he0) * n * AHe;
  double pHots = 1. / (1. + 77. / std::sqrt(T) * he0 / h0);

  // we multiplied Kenny's value with 1.e-12 to convert densities to m^-3
  // we then multiplied with 0.1 to convert to J m^-3s^-1
  double heatHeLa = pHots * 1.2196e-12 * alpha_e_2sP * ne * nhep * 1.e-12;
  gain = hfac * n *
         (cell.get_heating_H() * h0 + cell.get_heating_He() * AHe * he0);
  // pahs
  // we multiplied Kenny's value with 1.e-12 to convert densities to m^-3
  // we then multiplied with 0.1 to convert to J m^-3s^-1
  double heatpah = 3.e-38 * 5. * n * ne * pahfac;

  // cosmic rays
  // erg/cm^(9/2)/s --> J/m^(9/2)/s ==> 1.2e-27 --> 1.2e-25
  // value comes from equation (53) in Wiener, Zweibel & Oh, 2013, ApJ, 767, 87
  double heatcr = crfac * (1. - h0) * 1.2e-25 / std::sqrt(ne);

  gain += heatpah;
  gain += heatHeLa;
  gain += heatcr;

  // coolants

  // carbon
  double C21 = jfac * cell.get_mean_intensity(ION_C_p1) / ne / alphaC[0];
  // as can be seen below, CTHerecom has the same units as a recombination
  // rate: m^3s^-1
  // in Kenny's code, recombination rates are in cm^3s^-1
  // to put them in m^3s^-1 as well, we hence need to multiply Kenny's
  // original factor 1.e-9 with 1.e-6
  double CTHerecom = 1.e-15 * 0.046 * t4 * t4;
  double C32 = jfac * cell.get_mean_intensity(ION_C_p2) /
               (ne * alphaC[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(4, 6, T) +
                n * he0 * AHe * CTHerecom);
  double C31 = C32 * C21;
  double sumC = C21 + C31;
  cell.set_ionic_fraction(ION_C_p1, C21 / (1. + sumC));
  cell.set_ionic_fraction(ION_C_p2, C31 / (1. + sumC));

  // nitrogen
  double N21 = (jfac * cell.get_mean_intensity(ION_N_n) +
                nhp * ctr.get_charge_transfer_ionization_rate(1, 7, T)) /
               (ne * alphaN[0] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(2, 7, T));
  // multiplied Kenny's value with 1.e-6
  CTHerecom =
      1.e-15 * 0.33 * std::pow(t4, 0.29) * (1. + 1.3 * std::exp(-4.5 / t4));
  double N32 = jfac * cell.get_mean_intensity(ION_N_p1) /
               (ne * alphaN[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(3, 7, T) +
                n * he0 * AHe * CTHerecom);
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 1.e-15 * 0.15;
  double N43 = jfac * cell.get_mean_intensity(ION_N_p2) /
               (ne * alphaN[2] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(4, 7, T) +
                n * he0 * AHe * CTHerecom);
  double N31 = N32 * N21;
  double N41 = N43 * N31;
  double sumN = N21 + N31 + N41;
  cell.set_ionic_fraction(ION_N_n, N21 / (1. + sumN));
  cell.set_ionic_fraction(ION_N_p1, N31 / (1. + sumN));
  cell.set_ionic_fraction(ION_N_p2, N41 / (1. + sumN));

  // Sulphur
  double S21 = jfac * cell.get_mean_intensity(ION_S_p1) /
               (ne * alphaS[0] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(3, 16, T));
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 1.e-15 * 1.1 * std::pow(t4, 0.56);
  double S32 = jfac * cell.get_mean_intensity(ION_S_p2) /
               (ne * alphaS[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(4, 16, T) +
                n * he0 * AHe * CTHerecom);
  // multiplied Kenny's value with 1.e-6
  CTHerecom =
      1.e-15 * 7.6e-4 * std::pow(t4, 0.32) * (1. + 3.4 * std::exp(-5.25 * t4));
  double S43 = jfac * cell.get_mean_intensity(ION_S_p3) /
               (ne * alphaS[2] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(5, 16, T) +
                n * he0 * AHe * CTHerecom);
  double S31 = S32 * S21;
  double S41 = S43 * S31;
  double sumS = S21 + S31 + S41;
  cell.set_ionic_fraction(ION_S_p1, S21 / (1. + sumS));
  cell.set_ionic_fraction(ION_S_p2, S31 / (1. + sumS));
  cell.set_ionic_fraction(ION_S_p3, S41 / (1. + sumS));

  // Neon
  double Ne21 = jfac * cell.get_mean_intensity(ION_Ne_n) / (ne * alphaNe[0]);
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 1.e-15 * 1.e-5;
  double Ne32 = jfac * cell.get_mean_intensity(ION_Ne_p1) /
                (ne * alphaNe[1] +
                 n * h0 * ctr.get_charge_transfer_recombination_rate(3, 10, T) +
                 n * he0 * AHe * CTHerecom);
  double Ne31 = Ne32 * Ne21;
  double sumNe = Ne21 + Ne31;
  cell.set_ionic_fraction(ION_Ne_n, Ne21 / (1. + sumNe));
  cell.set_ionic_fraction(ION_Ne_p1, Ne31 / (1. + sumNe));

  // Oxygen
  double O21 = (jfac * cell.get_mean_intensity(ION_O_n) +
                nhp * ctr.get_charge_transfer_ionization_rate(1, 8, T)) /
               (ne * alphaO[0] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(2, 8, T));
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 0.2e-15 * std::pow(t4, 0.95);
  double O32 = jfac * cell.get_mean_intensity(ION_O_p1) /
               (ne * alphaO[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(3, 8, T) +
                n * he0 * AHe * CTHerecom);
  double O31 = O32 * O21;
  double sumO = O21 + O31;
  cell.set_ionic_fraction(ION_O_n, O21 / (1. + sumO));
  cell.set_ionic_fraction(ION_O_p1, O31 / (1. + sumO));

  double abund[12];
  abund[0] =
      abundances.get_abundance(ELEMENT_N) *
      (1. - cell.get_ionic_fraction(ION_N_n) -
       cell.get_ionic_fraction(ION_N_p1) - cell.get_ionic_fraction(ION_N_p2));
  abund[1] =
      abundances.get_abundance(ELEMENT_N) * cell.get_ionic_fraction(ION_N_n);
  abund[2] = abundances.get_abundance(ELEMENT_O) *
             (1. - cell.get_ionic_fraction(ION_O_n) -
              cell.get_ionic_fraction(ION_O_p1));
  abund[3] =
      abundances.get_abundance(ELEMENT_O) * cell.get_ionic_fraction(ION_O_n);
  abund[4] =
      abundances.get_abundance(ELEMENT_O) * cell.get_ionic_fraction(ION_O_p1);
  abund[5] =
      abundances.get_abundance(ELEMENT_Ne) * cell.get_ionic_fraction(ION_Ne_p1);
  abund[6] =
      abundances.get_abundance(ELEMENT_S) *
      (1. - cell.get_ionic_fraction(ION_S_p1) -
       cell.get_ionic_fraction(ION_S_p2) - cell.get_ionic_fraction(ION_S_p3));
  abund[7] =
      abundances.get_abundance(ELEMENT_S) * cell.get_ionic_fraction(ION_S_p1);
  abund[8] = abundances.get_abundance(ELEMENT_C) *
             (1. - cell.get_ionic_fraction(ION_C_p1) -
              cell.get_ionic_fraction(ION_C_p2));
  abund[9] =
      abundances.get_abundance(ELEMENT_C) * cell.get_ionic_fraction(ION_C_p1);
  abund[10] =
      abundances.get_abundance(ELEMENT_N) * cell.get_ionic_fraction(ION_N_p1);
  abund[11] =
      abundances.get_abundance(ELEMENT_Ne) * cell.get_ionic_fraction(ION_Ne_n);

  // FFcool
  double c = 5.5 - std::log(T);
  double gff = 1.1 + 0.34 * std::exp(-c * c / 3.);
  // we multiplied Kenny's value with 1.e-12 to convert the densities into m^-3
  // we then multiplied with 0.1 to convert them to J m^-3s^-1
  double Lff = 1.42e-40 * gff * std::sqrt(T) * (nhp + nhep) * ne;

  // RECcool
  // we multiplied Kenny's value with 1.e-12 to convert the densities into m^-3
  // we then multiplied with 0.1 to convert them to J m^-3s^-1
  double Lhp = 2.85e-14 * ne * nhp * std::sqrt(T) *
               (5.914 - 0.5 * std::log(T) + 0.01184 * std::pow(T, 0.33333));
  double Lhep = 2.6e-13 * ne * nhep * std::pow(T, 0.32);
  double LRec = 1.e-26 * (Lhp + Lhep);

  double Lc = data.get_cooling(T, ne, abund) * n;

  loss = Lc + Lff + LRec;
}

/**
 * @brief Calculate a new temperature for the given cell.
 *
 * @param jfac Normalization factor for the mean intensity integrals.
 * @param hfac Normalization factor for the heating integrals.
 * @param cell DensityGrid::iterator pointing to a cell.
 */
void TemperatureCalculator::calculate_temperature(
    double jfac, double hfac, DensityGrid::iterator &cell) const {
  const double eps = 1.e-3;
  const unsigned int max_iterations = 100;

  if ((cell.get_heating_H() == 0. && cell.get_mean_intensity(ION_He_n) == 0.) ||
      cell.get_number_density() == 0.) {
    cell.set_temperature(500.);

    cell.set_ionic_fraction(ION_H_n, 1.);

    cell.set_ionic_fraction(ION_He_n, 1.);

    cell.set_ionic_fraction(ION_C_p1, 0.);
    cell.set_ionic_fraction(ION_C_p2, 0.);

    cell.set_ionic_fraction(ION_N_n, 0.);
    cell.set_ionic_fraction(ION_N_p1, 0.);
    cell.set_ionic_fraction(ION_N_p2, 0.);

    cell.set_ionic_fraction(ION_O_n, 0.);
    cell.set_ionic_fraction(ION_O_p1, 0.);

    cell.set_ionic_fraction(ION_Ne_n, 0.);
    cell.set_ionic_fraction(ION_Ne_p1, 0.);

    cell.set_ionic_fraction(ION_S_p1, 0.);
    cell.set_ionic_fraction(ION_S_p2, 0.);
    cell.set_ionic_fraction(ION_S_p3, 0.);

    return;
  }

  double T0;
  if (cell.get_temperature() > 4000.) {
    T0 = cell.get_temperature();
  } else {
    T0 = 8000.;
  }

  unsigned int niter = 0;
  double gain0 = 1.;
  double loss0 = 0.;
  double h0 = 0.;
  double he0 = 0.;
  while (std::abs(gain0 - loss0) > eps * gain0 && niter < max_iterations) {
    ++niter;
    double T1 = 1.1 * T0;
    // ioneng
    double h01, he01, gain1, loss1;
    ioneng(h01, he01, gain1, loss1, T1, cell, jfac, _abundances, hfac, _pahfac,
           _crfac, _line_cooling_data, _recombination_rates,
           _charge_transfer_rates);

    double T2 = 0.9 * T0;
    // ioneng
    double h02, he02, gain2, loss2;
    ioneng(h02, he02, gain2, loss2, T2, cell, jfac, _abundances, hfac, _pahfac,
           _crfac, _line_cooling_data, _recombination_rates,
           _charge_transfer_rates);

    // ioneng - this one sets h0, he0, gain0 and loss0
    ioneng(h0, he0, gain0, loss0, T0, cell, jfac, _abundances, hfac, _pahfac,
           _crfac, _line_cooling_data, _recombination_rates,
           _charge_transfer_rates);

    double logtt = std::log(T1 / T2);
    double expgain = std::log(gain1 / gain2) / logtt;
    double exploss = std::log(loss1 / loss2) / logtt;
    T0 *= std::pow(loss0 / gain0, 1. / (expgain - exploss));

    if (T0 < 4000.) {
      // gas is neutral, temperature is 500 K
      T0 = 500.;
      h0 = 1.;
      he0 = 1.;
      // force exit out of loop
      gain0 = 1.;
      loss0 = 1.;
    }

    if (T0 > 1.e10) {
      // gas is ionized, temperature is 10^10 K (should probably be a lower
      // value)
      T0 = 1.e10;
      h0 = 0.;
      he0 = 0.;
      // force exit out of loop
      gain0 = 1.;
      loss0 = 1.;
    }
  }
  if (niter == max_iterations) {
    cmac_warning("Maximum number of iterations reached (temperature: %g, "
                 "relative difference cooling/heating: %g, aim: %g)!",
                 T0, std::abs(loss0 - gain0) / gain0, eps);
  }

  cell.set_temperature(T0);
  cell.set_ionic_fraction(ION_H_n, h0);
  cell.set_ionic_fraction(ION_He_n, he0);

  if (cell.get_mean_intensity(ION_H_n) == 0.) {
    cell.set_ionic_fraction(ION_H_n, 1.);
  }
  if (cell.get_mean_intensity(ION_He_n) == 0.) {
    cell.set_ionic_fraction(ION_He_n, 1.);
  }

  if (h0 == 1.) {
    cell.set_ionic_fraction(ION_C_p1, 0.);
    cell.set_ionic_fraction(ION_C_p2, 0.);

    cell.set_ionic_fraction(ION_N_n, 0.);
    cell.set_ionic_fraction(ION_N_p1, 0.);
    cell.set_ionic_fraction(ION_N_p2, 0.);

    cell.set_ionic_fraction(ION_O_n, 0.);
    cell.set_ionic_fraction(ION_O_p1, 0.);

    cell.set_ionic_fraction(ION_Ne_n, 0.);
    cell.set_ionic_fraction(ION_Ne_p1, 0.);

    cell.set_ionic_fraction(ION_S_p1, 0.);
    cell.set_ionic_fraction(ION_S_p2, 0.);
    cell.set_ionic_fraction(ION_S_p3, 0.);
  }
}

/**
 * @brief Calculate a new temperature for each cell after shooting the given
 * number of photons.
 *
 * @param totweight Total weight of all photons that were used.
 * @param grid DensityGrid on which to operate.
 * @param block Block that should be traversed by the local MPI process.
 */
void TemperatureCalculator::calculate_temperature(
    double totweight, DensityGrid &grid,
    std::pair< unsigned long, unsigned long > &block) const {
  double jfac = _luminosity / totweight;
  // the integral calculation uses the photon frequency (in Hz)
  // we want to convert this to the photon energy (in Joule)
  // we do this by multiplying with the Planck constant (in Js)
  double hfac = jfac * 6.626070040e-34;

  WorkDistributor<
      DensityGridTraversalJobMarket< TemperatureCalculatorFunction >,
      DensityGridTraversalJob< TemperatureCalculatorFunction > >
      workers;
  TemperatureCalculatorFunction do_calculation(*this, jfac, hfac);
  DensityGridTraversalJobMarket< TemperatureCalculatorFunction > jobs(
      grid, do_calculation, block);
  workers.do_in_parallel(jobs);
}
