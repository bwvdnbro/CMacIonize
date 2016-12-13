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
 * @file EmissivityCalculator.cpp
 *
 * @brief EmissivityCalculator implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "EmissivityCalculator.hpp"
#include "Abundances.hpp"
#include "DensityGrid.hpp"
#include "DensityValues.hpp"
#include "LineCoolingData.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Initializes hydrogen and helium continuous emission coefficient tables.
 *
 * @param abundances Abundances.
 */
EmissivityCalculator::EmissivityCalculator(Abundances &abundances)
    : _abundances(abundances) {
  // these values come from Brown & Mathew, 1970, ApJ, 160, 939
  // the hmit and heplt tables correspond to wavelength 3646 in table 1 and
  // wavelength 3680 in table 5 respectively
  // the hplt and hemit tables correspond to the values you get when you
  // linearly interpolate in table 1 to wavelength 3681 and in table 5 to
  // wavelength 3646 respectively
  // all values are in 1.e-40 erg cm^3s^-1Hz^-1 (except for the ttab values,
  // which are in K).
  double ttab[8] = {4.e3, 6.e3, 8.e3, 1.e4, 1.2e4, 1.4e4, 1.6e4, 1.8e4};
  double hplt[8] = {0.162, 0.584, 1.046, 1.437, 1.742, 1.977, 2.159, 2.297};
  double hmit[8] = {92.6, 50.9, 33.8, 24.8, 19.53, 16.09, 13.7, 11.96};
  double heplt[8] = {0.189, 0.622, 1.076, 1.45, 1.74, 1.963, 2.14, 2.27};
  double hemit[8] = {15.7, 9.23, 6.71, 5.49, 4.83, 4.41, 4.135, 3.94};

  for (unsigned int i = 0; i < 8; ++i) {
    _logttab[i] = std::log(ttab[i]);
    _loghplt[i] = std::log(hplt[i]);
    _loghmit[i] = std::log(hmit[i]);
    _logheplt[i] = std::log(heplt[i]);
    _loghemit[i] = std::log(hemit[i]);
  }
}

/**
 * @brief Calculate the hydrogen and helium continuous emission coefficients at
 * the given temperature.
 *
 * @param T Temperature (in K).
 * @param emhpl Hydrogen coefficient 1 (in J m^3s^-1angstrom^-1).
 * @param emhmi Hydrogen coefficient 2 (in J m^3s^-1angstrom^-1).
 * @param emhepl Helium coefficient 1 (in J m^3s^-1angstrom^-1).
 * @param emhemi Helium coefficient 2 (in J m^3s^-1angstrom^-1).
 */
void EmissivityCalculator::bjump(double T, double &emhpl, double &emhmi,
                                 double &emhepl, double &emhemi) const {
  double logt = std::log(T);

  int i = Utilities::locate(logt, _logttab, 8);
  i = std::max(i, 0);
  i = std::min(i, 6);

  emhpl = _loghplt[i] +
          (logt - _logttab[i]) * (_loghplt[i + 1] - _loghplt[i]) /
              (_logttab[i + 1] - _logttab[i]);
  emhmi = _loghmit[i] +
          (logt - _logttab[i]) * (_loghmit[i + 1] - _loghmit[i]) /
              (_logttab[i + 1] - _logttab[i]);
  emhepl = _logheplt[i] +
           (logt - _logttab[i]) * (_logheplt[i + 1] - _logheplt[i]) /
               (_logttab[i + 1] - _logttab[i]);
  emhemi = _loghemit[i] +
           (logt - _logttab[i]) * (_loghemit[i + 1] - _loghemit[i]) /
               (_logttab[i + 1] - _logttab[i]);

  emhpl = std::exp(emhpl);
  emhmi = std::exp(emhmi);
  emhepl = std::exp(emhepl);
  emhemi = std::exp(emhemi);
  // the values above are in 1.e-40 erg cm^3s^-1Hz^-1
  // convert to J m^3s^-1angstrom^-1
  emhpl *= 1.e-53 * 2.99792458e18 / 3681. / 3681.;
  emhmi *= 1.e-53 * 2.99792458e18 / 3643. / 3643.;
  emhepl *= 1.e-53 * 2.99792458e18 / 3681. / 3681.;
  emhemi *= 1.e-53 * 2.99792458e18 / 3643. / 3643.;
}

/**
 * @brief Calculate the emissivity values for a single cell.
 *
 * @param cell DensityValues of the cell.
 * @param abundances Abundances.
 * @param lines LineCoolingData used to calculate emission line strengths.
 * @return EmissivityValues in the cell.
 */
EmissivityValues EmissivityCalculator::calculate_emissivities(
    DensityValues &cell, Abundances &abundances,
    const LineCoolingData &lines) const {
  const double h0max = 0.2;

  EmissivityValues eval;

  if (cell.get_ionic_fraction(ION_H_n) < h0max &&
      cell.get_temperature() > 3000.) {
    double ntot = cell.get_total_density();
    double nhp = ntot * (1. - cell.get_ionic_fraction(ION_H_n));
    double nhep = ntot * (1. - cell.get_ionic_fraction(ION_He_n)) *
                  abundances.get_abundance(ELEMENT_He);
    double ne = nhp + nhep;

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
    abund[5] = abundances.get_abundance(ELEMENT_Ne) *
               cell.get_ionic_fraction(ION_Ne_p1);
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
    abund[11] = abundances.get_abundance(ELEMENT_Ne) *
                cell.get_ionic_fraction(ION_Ne_n);

    double c6300 = 0.;
    double c9405 = 0.;
    double c6312 = 0.;
    double c33mu = 0.;
    double c19mu = 0.;
    double c3729 = 0.;
    double c3727 = 0.;
    double c7330 = 0.;
    double c4363 = 0.;
    double c5007 = 0.;
    double c52mu = 0.;
    double c5755 = 0.;
    double c6584 = 0.;
    double c4072 = 0.;
    double c6717 = 0.;
    double c6725 = 0.;
    double c3869 = 0.;
    double cniii57 = 0.;
    double cneii12 = 0.;
    double cneiii15 = 0.;
    double cnii122 = 0.;
    double cii2325 = 0.;
    double ciii1908 = 0.;
    double coii7325 = 0.;
    double csiv10 = 0.;
    double c88mu = 0.;
    lines.linestr(cell.get_temperature(), ne, abund, c6300, c9405, c6312, c33mu,
                  c19mu, c3729, c3727, c7330, c4363, c5007, c52mu, c88mu, c5755,
                  c6584, c4072, c6717, c6725, c3869, cniii57, cneii12, cneiii15,
                  cnii122, cii2325, ciii1908, coii7325, csiv10);

    double t4 = cell.get_temperature() * 1.e-4;
    // we added correction factors 1.e-12 to convert densities to cm^-3
    // and an extra factor 1.e-1 to convert to J m^-3s^-1
    eval.set_emissivity(EMISSIONLINE_HAlpha,
                        ne * nhp * 2.87 * 1.24e-38 * std::pow(t4, -0.938));
    eval.set_emissivity(EMISSIONLINE_HBeta,
                        ne * nhp * 1.24e-38 * std::pow(t4, -0.878));
    eval.set_emissivity(EMISSIONLINE_HII,
                        nhp * ne * 4.9e-40 * std::pow(t4, -0.848));

    double emhpl = 0.;
    double emhmi = 0.;
    double emhepl = 0.;
    double emhemi = 0.;
    bjump(cell.get_temperature(), emhpl, emhmi, emhepl, emhemi);
    eval.set_emissivity(EMISSIONLINE_BALMER_JUMP_LOW,
                        ne * (nhp * emhmi + nhep * emhemi));
    eval.set_emissivity(EMISSIONLINE_BALMER_JUMP_HIGH,
                        ne * (nhp * emhpl + nhep * emhepl));

    eval.set_emissivity(EMISSIONLINE_OI_6300, ntot * c6300);
    eval.set_emissivity(EMISSIONLINE_OII_3727, ntot * c3727);
    eval.set_emissivity(EMISSIONLINE_OIII_5007, ntot * c5007);
    eval.set_emissivity(EMISSIONLINE_OIII_4363, ntot * c4363);
    eval.set_emissivity(EMISSIONLINE_OIII_88mu, ntot * c88mu);
    eval.set_emissivity(EMISSIONLINE_NII_5755, ntot * c5755);
    eval.set_emissivity(EMISSIONLINE_NII_6584, ntot * c6584);
    eval.set_emissivity(EMISSIONLINE_NeIII_3869, ntot * c3869);
    eval.set_emissivity(EMISSIONLINE_SII_6725, ntot * c6725);
    eval.set_emissivity(EMISSIONLINE_SII_4072, ntot * c4072);
    eval.set_emissivity(EMISSIONLINE_SIII_9405, ntot * c9405);
    eval.set_emissivity(EMISSIONLINE_SIII_6312, ntot * c6312);
    eval.set_emissivity(EMISSIONLINE_SIII_19mu, ntot * c19mu);
    double T = cell.get_temperature();
    eval.set_emissivity(EMISSIONLINE_avg_T, ne * nhp * T);
    eval.set_emissivity(EMISSIONLINE_avg_T_count, ne * nhp);
    eval.set_emissivity(EMISSIONLINE_avg_nH_nHe,
                        (1. - cell.get_ionic_fraction(ION_H_n)) *
                            (1. - cell.get_ionic_fraction(ION_He_n)));
    eval.set_emissivity(EMISSIONLINE_avg_nH_nHe_count, 1.);
    eval.set_emissivity(EMISSIONLINE_NeII_12mu, ntot * cneii12);
    eval.set_emissivity(EMISSIONLINE_NIII_57mu, ntot * cniii57);
    eval.set_emissivity(EMISSIONLINE_NeIII_15mu, ntot * cneiii15);
    eval.set_emissivity(EMISSIONLINE_NII_122mu, ntot * cnii122);
    eval.set_emissivity(EMISSIONLINE_CII_2325, ntot * cii2325);
    eval.set_emissivity(EMISSIONLINE_CIII_1908, ntot * ciii1908);
    eval.set_emissivity(EMISSIONLINE_OII_7325, ntot * coii7325);
    eval.set_emissivity(EMISSIONLINE_SIV_10mu, ntot * csiv10);
    // we converted Kenny's constant from 1.e20 erg/cm^6/s to J/m^6/s
    eval.set_emissivity(EMISSIONLINE_HeI_5876,
                        ne * nhep * 1.69e-38 * std::pow(t4, -1.065));
    eval.set_emissivity(EMISSIONLINE_Hrec_s,
                        ne * nhp * 7.982e-23 /
                            (std::sqrt(T / 3.148) *
                             std::pow(1. + std::sqrt(T / 3.148), 0.252) *
                             std::pow(1. + std::sqrt(T / 7.036e5), 1.748)));
  }

  return eval;
}

/**
 * @brief Calculate the emissivities for all cells in the given DensityGrid.
 *
 * @param grid DensityGrid to operate on.
 */
void EmissivityCalculator::calculate_emissivities(DensityGrid &grid) const {
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    DensityValues &cell = it.get_values();
    EmissivityValues *emissivities =
        new EmissivityValues(calculate_emissivities(cell, _abundances, _lines));
    cell.set_emissivities(emissivities);
  }
}

/**
 * @brief Get the emissivities for all cells in the given DensityGrid.
 *
 * @param grid DensityGrid to operate on.
 * @return std::vector containing EmissivityValues, one for each cell in the
 * grid (in the same order the grid is traversed).
 */
std::vector< EmissivityValues >
EmissivityCalculator::get_emissivities(DensityGrid &grid) const {
  std::vector< EmissivityValues > result(grid.get_number_of_cells());
  unsigned int index = 0;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    DensityValues &cell = it.get_values();
    result[index] = calculate_emissivities(cell, _abundances, _lines);
    ++index;
  }
  return result;
}
