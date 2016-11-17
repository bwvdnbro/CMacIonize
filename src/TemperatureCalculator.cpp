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
#include "ChargeTransferRates.hpp"
#include "DensityValues.hpp"
#include "IonizationStateCalculator.hpp"
#include "LineCoolingData.hpp"
#include "RecombinationRates.hpp"
#include <cmath>

void TemperatureCalculator::ioneng(
    double &h0, double &he0, double &gain, double &loss, DensityValues &cell,
    double jfac, double AHe, double AC, double AN, double AO, double AS,
    double ANe, double hfac, double pahfac, LineCoolingData &data,
    RecombinationRates &rates, ChargeTransferRates &ctr) {
  double T = cell.get_temperature();

  double alphaH = rates.get_recombination_rate(ELEMENT_H, T);
  double alphaHe = rates.get_recombination_rate(ELEMENT_He, T);
  double alphaC[2];
  alphaC[0] = rates.get_recombination_rate(ELEMENT_Cp1, T);
  alphaC[1] = rates.get_recombination_rate(ELEMENT_Cp2, T);
  double alphaN[3];
  alphaN[0] = rates.get_recombination_rate(ELEMENT_N, T);
  alphaN[1] = rates.get_recombination_rate(ELEMENT_Np1, T);
  alphaN[2] = rates.get_recombination_rate(ELEMENT_Np2, T);
  double alphaO[2];
  alphaO[0] = rates.get_recombination_rate(ELEMENT_O, T);
  alphaO[1] = rates.get_recombination_rate(ELEMENT_Op1, T);
  double alphaNe[2];
  alphaNe[0] = rates.get_recombination_rate(ELEMENT_Ne, T);
  alphaNe[1] = rates.get_recombination_rate(ELEMENT_Nep1, T);
  double alphaS[3];
  alphaS[0] = rates.get_recombination_rate(ELEMENT_Sp1, T);
  alphaS[1] = rates.get_recombination_rate(ELEMENT_Sp2, T);
  alphaS[2] = rates.get_recombination_rate(ELEMENT_Sp3, T);

  double t4 = T * 1.e-4;
  double alpha_e_2sP = 4.27e-14 * std::pow(t4, -0.695);
  double n = cell.get_total_density();

  double jH = jfac * cell.get_mean_intensity(ELEMENT_H);
  double jHe = jfac * cell.get_mean_intensity(ELEMENT_He);
  IonizationStateCalculator::find_H0(alphaH, alphaHe, jH, jHe, n, AHe, T, h0,
                                     he0);

  double ne = n * (1. - h0 + AHe * (1. - he0));
  double nhp = n * (1. - h0);
  double nhep = (1. - he0) * n * AHe;
  double pHots = 1. / (1. + 77. / std::sqrt(T) * he0 / h0);

  // UNITS???????
  double heatHeLa = pHots * 1.2196e-11 * alpha_e_2sP * ne * nhep * 1.e20;
  gain = hfac * n *
         (cell.get_heating_H() * h0 + cell.get_heating_He() * AHe * he0);
  // pahs
  // the value of G1 is 0 in Kenny's code...
  double G1 = 0.;
  double heatpah = 3.e-25 * 5. * n * ne * pahfac;
  gain += heatpah;
  gain *= 1.e20;
  gain += heatHeLa + G1;

  // coolants

  // carbon
  double C21 = jfac * cell.get_mean_intensity(ELEMENT_Cp1) / ne / alphaC[0];
  // as can be seen below, CTHerecom has the same units as a recombination
  // rate: m^3s^-1
  // in Kenny's code, recombination rates are in cm^3s^-1
  // to put them in m^3s^-1 as well, we hence need to multiply Kenny's
  // original factor 1.e-9 with 1.e-6
  double CTHerecom = 1.e-15 * 0.046 * t4 * t4;
  double C32 = jfac * cell.get_mean_intensity(ELEMENT_Cp2) /
               (ne * alphaC[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(4, 6, T) +
                n * he0 * AHe * CTHerecom);
  double C31 = C32 * C21;
  double sumC = C21 + C31;
  cell.set_ionic_fraction(ELEMENT_Cp1, C21 / (1. + sumC));
  cell.set_ionic_fraction(ELEMENT_Cp2, C31 / (1. + sumC));

  // nitrogen
  double N21 = (jfac * cell.get_mean_intensity(ELEMENT_N) +
                nhp * ctr.get_charge_transfer_ionization_rate(1, 7, T)) /
               (ne * alphaN[0] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(2, 7, T));
  // multiplied Kenny's value with 1.e-6
  CTHerecom =
      1.e-15 * 0.33 * std::pow(t4, 0.29) * (1. + 1.3 * std::exp(-4.5 / t4));
  double N32 = jfac * cell.get_mean_intensity(ELEMENT_Np1) /
               (ne * alphaN[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(3, 7, T) +
                n * he0 * AHe * CTHerecom);
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 1.e-15 * 0.15;
  double N43 = jfac * cell.get_mean_intensity(ELEMENT_Np2) /
               (ne * alphaN[2] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(4, 7, T) +
                n * he0 * AHe * CTHerecom);
  double N31 = N32 * N21;
  double N41 = N43 * N31;
  double sumN = N21 + N31 + N41;
  cell.set_ionic_fraction(ELEMENT_N, N21 / (1. + sumN));
  cell.set_ionic_fraction(ELEMENT_Np1, N31 / (1. + sumN));
  cell.set_ionic_fraction(ELEMENT_Np2, N41 / (1. + sumN));

  // Sulfur
  double S21 = jfac * cell.get_mean_intensity(ELEMENT_Sp1) /
               (ne * alphaS[0] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(3, 16, T));
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 1.e-15 * 1.1 * std::pow(t4, 0.56);
  double S32 = jfac * cell.get_mean_intensity(ELEMENT_Sp2) /
               (ne * alphaS[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(4, 16, T) +
                n * he0 * AHe * CTHerecom);
  // multiplied Kenny's value with 1.e-6
  CTHerecom =
      1.e-15 * 7.6e-4 * std::pow(t4, 0.32) * (1. + 3.4 * std::exp(-5.25 * t4));
  double S43 = jfac * cell.get_mean_intensity(ELEMENT_Sp3) /
               (ne * alphaS[2] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(5, 16, T) +
                n * he0 * AHe * CTHerecom);
  double S31 = S32 * S21;
  double S41 = S43 * S31;
  double sumS = S21 + S31 + S41;
  cell.set_ionic_fraction(ELEMENT_Sp1, S21 / (1. + sumS));
  cell.set_ionic_fraction(ELEMENT_Sp2, S31 / (1. + sumS));
  cell.set_ionic_fraction(ELEMENT_Sp3, S41 / (1. + sumS));

  // Neon
  double Ne21 = jfac * cell.get_mean_intensity(ELEMENT_Ne) / (ne * alphaNe[0]);
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 1.e-15 * 1.e-5;
  double Ne32 = jfac * cell.get_mean_intensity(ELEMENT_Nep1) /
                (ne * alphaNe[1] +
                 n * h0 * ctr.get_charge_transfer_recombination_rate(3, 10, T) +
                 n * he0 * AHe * CTHerecom);
  double Ne31 = Ne32 * Ne21;
  double sumNe = Ne21 + Ne31;
  cell.set_ionic_fraction(ELEMENT_Ne, Ne21 / (1. + sumNe));
  cell.set_ionic_fraction(ELEMENT_Nep1, Ne31 / (1. + sumNe));

  // Oxygen
  double O21 = (jfac * cell.get_mean_intensity(ELEMENT_O) +
                nhp * ctr.get_charge_transfer_ionization_rate(1, 8, T)) /
               (ne * alphaO[0] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(2, 8, T));
  // multiplied Kenny's value with 1.e-6
  CTHerecom = 0.2e-15 * std::pow(t4, 0.95);
  double O32 = jfac * cell.get_mean_intensity(ELEMENT_Op1) /
               (ne * alphaO[1] +
                n * h0 * ctr.get_charge_transfer_recombination_rate(3, 8, T) +
                n * he0 * AHe * CTHerecom);
  double O31 = O32 * O21;
  double sumO = O21 + O31;
  cell.set_ionic_fraction(ELEMENT_O, O21 / (1. + sumO));
  cell.set_ionic_fraction(ELEMENT_Op1, O31 / (1. + sumO));
}

/**
 * @brief Calculate a new temperature for the given cell.
 *
 * @param jfac Normalization factor for the mean intensity integrals.
 * @param cell DensityValues of the cell.
 */
void TemperatureCalculator::calculate_temperature(double jfac,
                                                  DensityValues &cell) {
  const double eps = 1.e-3;

  if ((cell.get_heating_H() == 0. &&
       cell.get_mean_intensity(ELEMENT_He) == 0.) ||
      cell.get_total_density() == 0.) {
    cell.set_temperature(500.);

    cell.set_ionic_fraction(ELEMENT_H, 1.);

    cell.set_ionic_fraction(ELEMENT_He, 1.);

    cell.set_ionic_fraction(ELEMENT_Cp1, 0.);
    cell.set_ionic_fraction(ELEMENT_Cp2, 0.);

    cell.set_ionic_fraction(ELEMENT_N, 0.);
    cell.set_ionic_fraction(ELEMENT_Np1, 0.);
    cell.set_ionic_fraction(ELEMENT_Np2, 0.);

    cell.set_ionic_fraction(ELEMENT_O, 0.);
    cell.set_ionic_fraction(ELEMENT_Op1, 0.);

    cell.set_ionic_fraction(ELEMENT_Ne, 0.);
    cell.set_ionic_fraction(ELEMENT_Nep1, 0.);

    cell.set_ionic_fraction(ELEMENT_Sp1, 0.);
    cell.set_ionic_fraction(ELEMENT_Sp2, 0.);
    cell.set_ionic_fraction(ELEMENT_Sp3, 0.);

    return;
  }

  double T0;
  if (cell.get_temperature() > 4000.) {
    T0 = cell.get_temperature();
  } else {
    T0 = 8000.;
  }

  unsigned int niter = 0;
  double gain0 = 0.;
  double loss0 = 0.;
  double h0 = 0.;
  double he0 = 0.;
  while (std::abs(gain0 - loss0) > eps * gain0 && niter <= 10) {
    ++niter;
    double T1 = 1.1 * T0;
    // ioneng
    double gain1 = 0.;
    double loss1 = 0.;

    double T2 = 0.9 * T0;
    // ioneng
    double gain2 = 0.;
    double loss2 = 0.;

    double expgain = std::log(gain1 / gain2) / std::log(T1 / T2);
    double exploss = std::log(loss1 / loss2) / std::log(T1 / T2);
    double Tnew = T0 * std::pow(loss0 / gain0, 1. / (expgain - exploss));
    T0 = Tnew;
    if (T0 < 4000.) {
      T0 = 500.;
      h0 = 1.;
      he0 = 1.;
      gain0 = 1.;
      loss0 = 1.;
    }
  }

  if (cell.get_mean_intensity(ELEMENT_H) == 0.) {
    h0 = 1.;
  }
  if (cell.get_mean_intensity(ELEMENT_He) == 0.) {
    he0 = 1.;
  }

  cell.set_temperature(T0);
  cell.set_ionic_fraction(ELEMENT_H, h0);
  cell.set_ionic_fraction(ELEMENT_He, he0);

  if (h0 == 1.) {
    cell.set_ionic_fraction(ELEMENT_Cp1, 0.);
    cell.set_ionic_fraction(ELEMENT_Cp2, 0.);

    cell.set_ionic_fraction(ELEMENT_N, 0.);
    cell.set_ionic_fraction(ELEMENT_Np1, 0.);
    cell.set_ionic_fraction(ELEMENT_Np2, 0.);

    cell.set_ionic_fraction(ELEMENT_O, 0.);
    cell.set_ionic_fraction(ELEMENT_Op1, 0.);

    cell.set_ionic_fraction(ELEMENT_Ne, 0.);
    cell.set_ionic_fraction(ELEMENT_Nep1, 0.);

    cell.set_ionic_fraction(ELEMENT_Sp1, 0.);
    cell.set_ionic_fraction(ELEMENT_Sp2, 0.);
    cell.set_ionic_fraction(ELEMENT_Sp3, 0.);
  }
}

/**
 * @brief Calculate a new temperature for each cell after shooting the given
 * number of photons.
 *
 * @param nphoton Number of ionizing photons emitted by all sources.
 * @param grid DensityGrid on which to operate.
 */
void TemperatureCalculator::calculate_temperature(unsigned int nphoton,
                                                  DensityGrid &grid) {
  // do stuff
}
