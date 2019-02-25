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
#include "PhysicalConstants.hpp"
#include "Utilities.hpp"
#include <cmath>

/**
 * @brief Constructor.
 *
 * Initializes hydrogen and helium continuous emission coefficient tables.
 *
 * @param abundances Abundances.
 */
EmissivityCalculator::EmissivityCalculator(const Abundances &abundances)
    : _abundances(abundances) {
  // these values come from Brown & Mathews (1970)
  // the hmit and heplt tables correspond to wavelength 3646 in table 1 and
  // wavelength 3680 in table 5 respectively
  // the hplt and hemit tables correspond to the values you get when you
  // linearly interpolate in table 1 to wavelength 3681 and in table 5 to
  // wavelength 3646 respectively
  // all values are in 1.e-40 erg cm^3 s^-1 Hz^-1 (except for the ttab values,
  // which are in K).
  const double ttab[8] = {4.e3, 6.e3, 8.e3, 1.e4, 1.2e4, 1.4e4, 1.6e4, 1.8e4};
  const double hplt[8] = {0.162, 0.584, 1.046, 1.437,
                          1.742, 1.977, 2.159, 2.297};
  const double hmit[8] = {92.6, 50.9, 33.8, 24.8, 19.53, 16.09, 13.7, 11.96};
  const double heplt[8] = {0.189, 0.622, 1.076, 1.45, 1.74, 1.963, 2.14, 2.27};
  const double hemit[8] = {15.7, 9.23, 6.71, 5.49, 4.83, 4.41, 4.135, 3.94};

  for (uint_fast8_t i = 0; i < 8; ++i) {
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
 * @param emission_hydrogen_high Hydrogen coefficient 1 (in J m^3 s^-1
 * angstrom^-1).
 * @param emission_hydrogen_low Hydrogen coefficient 2 (in J m^3 s^-1
 * angstrom^-1).
 * @param emission_helium_high Helium coefficient 1 (in J m^3 s^-1 angstrom^-1).
 * @param emission_helium_low Helium coefficient 2 (in J m^3 s^-1 angstrom^-1).
 */
void EmissivityCalculator::get_balmer_jump_emission(
    double T, double &emission_hydrogen_high, double &emission_hydrogen_low,
    double &emission_helium_high, double &emission_helium_low) const {

  const double logt = std::log(T);

  int_fast32_t i = Utilities::locate(logt, _logttab, 8);
  i = std::max(i, int_fast32_t(0));
  i = std::min(i, int_fast32_t(6));

  emission_hydrogen_high = _loghplt[i] + (logt - _logttab[i]) *
                                             (_loghplt[i + 1] - _loghplt[i]) /
                                             (_logttab[i + 1] - _logttab[i]);
  emission_hydrogen_low = _loghmit[i] + (logt - _logttab[i]) *
                                            (_loghmit[i + 1] - _loghmit[i]) /
                                            (_logttab[i + 1] - _logttab[i]);
  emission_helium_high = _logheplt[i] + (logt - _logttab[i]) *
                                            (_logheplt[i + 1] - _logheplt[i]) /
                                            (_logttab[i + 1] - _logttab[i]);
  emission_helium_low = _loghemit[i] + (logt - _logttab[i]) *
                                           (_loghemit[i + 1] - _loghemit[i]) /
                                           (_logttab[i + 1] - _logttab[i]);

  emission_hydrogen_high = std::exp(emission_hydrogen_high);
  emission_hydrogen_low = std::exp(emission_hydrogen_low);
  emission_helium_high = std::exp(emission_helium_high);
  emission_helium_low = std::exp(emission_helium_low);
  // the values above are in 1.e-40 erg cm^3 s^-1 Hz^-1
  // convert to J m^3 s^-1 angstrom^-1
  const double lightspeed =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED);
  emission_hydrogen_high *= 1.e-43 * lightspeed / (3681. * 3681.);
  emission_hydrogen_low *= 1.e-43 * lightspeed / (3643. * 3643.);
  emission_helium_high *= 1.e-43 * lightspeed / (3681. * 3681.);
  emission_helium_low *= 1.e-43 * lightspeed / (3643. * 3643.);
}

/**
 * @brief Calculate the emissivity values for a single cell.
 *
 * @param ionization_variables IonizationVariables of the cell.
 * @param abundances Abundances.
 * @param line_cooling_data LineCoolingData used to calculate emission line
 * strengths.
 * @return EmissivityValues in the cell.
 */
EmissivityValues EmissivityCalculator::calculate_emissivities(
    const IonizationVariables &ionization_variables,
    const Abundances &abundances,
    const LineCoolingData &line_cooling_data) const {

  const double h0max = 0.2;

  EmissivityValues eval;

  if (ionization_variables.get_ionic_fraction(ION_H_n) < h0max &&
      ionization_variables.get_temperature() > 3000.) {
    const double ntot = ionization_variables.get_number_density();
    const double nhp =
        ntot * (1. - ionization_variables.get_ionic_fraction(ION_H_n));
#ifdef HAS_HELIUM
    const double nhep =
        ntot * (1. - ionization_variables.get_ionic_fraction(ION_He_n)) *
        abundances.get_abundance(ELEMENT_He);
#else
    const double nhep = 0.;
#endif
    const double ne = nhp + nhep;

    // get the abundances of the ions used by the line cooling computation
    double abund[LINECOOLINGDATA_NUMELEMENTS];

#ifdef HAS_CARBON
    // carbon
    // we assume that all carbon is either C+, C++, or C+++
    // we only use C+ and C++
    // note that the ionic fraction of C_p1 corresponds to the fraction of
    // ionized C+, i.e. the fraction of C++
    abund[CII] = abundances.get_abundance(ELEMENT_C) *
                 (1. - ionization_variables.get_ionic_fraction(ION_C_p1) -
                  ionization_variables.get_ionic_fraction(ION_C_p2));
    abund[CIII] = abundances.get_abundance(ELEMENT_C) *
                  ionization_variables.get_ionic_fraction(ION_C_p1);
#endif

#ifdef HAS_NITROGEN
    // nitrogen
    // we assume all nitrogen is either N0, N+, N++ or N+++
    // we only use N0, N+ and N++
    abund[NI] = abundances.get_abundance(ELEMENT_N) *
                (1. - ionization_variables.get_ionic_fraction(ION_N_n) -
                 ionization_variables.get_ionic_fraction(ION_N_p1) -
                 ionization_variables.get_ionic_fraction(ION_N_p2));
    abund[NII] = abundances.get_abundance(ELEMENT_N) *
                 ionization_variables.get_ionic_fraction(ION_N_n);
    abund[NIII] = abundances.get_abundance(ELEMENT_N) *
                  ionization_variables.get_ionic_fraction(ION_N_p1);
#endif

#ifdef HAS_OXYGEN
    // oxygen
    // we assume all oxygen is either O0, O+ or O++
    // we use all of them
    abund[OI] = abundances.get_abundance(ELEMENT_O) *
                (1. - ionization_variables.get_ionic_fraction(ION_O_n) -
                 ionization_variables.get_ionic_fraction(ION_O_p1));
    abund[OII] = abundances.get_abundance(ELEMENT_O) *
                 ionization_variables.get_ionic_fraction(ION_O_n);
    abund[OIII] = abundances.get_abundance(ELEMENT_O) *
                  ionization_variables.get_ionic_fraction(ION_O_p1);
#endif

#ifdef HAS_NEON
    // neon
    // we make no assumptions on the relative abundances of different neon ions
    // we only use Ne+ and Ne++
    abund[NeII] = abundances.get_abundance(ELEMENT_Ne) *
                  ionization_variables.get_ionic_fraction(ION_Ne_n);
    abund[NeIII] = abundances.get_abundance(ELEMENT_Ne) *
                   ionization_variables.get_ionic_fraction(ION_Ne_p1);
#endif

#ifdef HAS_SULPHUR
    // sulphur
    // we assume all sulphur is either S+, S++, S+++ or S++++
    // we only use S+, S++ and S+++
    abund[SII] = abundances.get_abundance(ELEMENT_S) *
                 (1. - ionization_variables.get_ionic_fraction(ION_S_p1) -
                  ionization_variables.get_ionic_fraction(ION_S_p2) -
                  ionization_variables.get_ionic_fraction(ION_S_p3));
    abund[SIII] = abundances.get_abundance(ELEMENT_S) *
                  ionization_variables.get_ionic_fraction(ION_S_p1);
    abund[SIV] = abundances.get_abundance(ELEMENT_S) *
                 ionization_variables.get_ionic_fraction(ION_S_p2);
#endif

    std::vector< std::vector< double > > line_strengths =
        line_cooling_data.get_line_strengths(
            ionization_variables.get_temperature(), ne, abund);

    const double T = ionization_variables.get_temperature();
    const double T4 = T * 1.e-4;
    // we added correction factors 1.e-12 to convert densities to cm^-3
    // and an extra factor 1.e-1 to convert to J m^-3s^-1
    // Osterbrock & Ferland (2006), table 4.1
    eval.set_emissivity(EMISSIONLINE_HAlpha,
                        ne * nhp * 2.87 * 1.24e-38 * std::pow(T4, -0.938));
    // fits to Storey & Hummer (1995) data...
    eval.set_emissivity(EMISSIONLINE_HBeta,
                        ne * nhp * 1.24e-38 * std::pow(T4, -0.878));
    eval.set_emissivity(EMISSIONLINE_HII,
                        nhp * ne * 4.9e-40 * std::pow(T4, -0.848));

    double emhpl = 0.;
    double emhmi = 0.;
    double emhepl = 0.;
    double emhemi = 0.;
    get_balmer_jump_emission(ionization_variables.get_temperature(), emhpl,
                             emhmi, emhepl, emhemi);
    eval.set_emissivity(EMISSIONLINE_BALMER_JUMP_LOW,
                        ne * (nhp * emhmi + nhep * emhemi));
    eval.set_emissivity(EMISSIONLINE_BALMER_JUMP_HIGH,
                        ne * (nhp * emhpl + nhep * emhepl));

    // NII
    // Osterbrock & Ferland (2006), table 3.12
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    eval.set_emissivity(EMISSIONLINE_NII_5755,
                        ntot * line_strengths[NII][TRANSITION_3_to_4]);
    eval.set_emissivity(EMISSIONLINE_NII_6548,
                        ntot * line_strengths[NII][TRANSITION_1_to_3]);
    eval.set_emissivity(EMISSIONLINE_NII_6584,
                        ntot * line_strengths[NII][TRANSITION_2_to_3]);
    eval.set_emissivity(EMISSIONLINE_NII_122mu,
                        ntot * line_strengths[NII][TRANSITION_1_to_2]);

    // OI
    // Osterbrock & Ferland (2006), table 3.14
    // ground state: 3P2
    // excited states: 3P1, 3P0, 1D2, 1S0
    // this is the sum of the 6300.3 and 6363.8 angstrom transitions
    eval.set_emissivity(EMISSIONLINE_OI_6300,
                        ntot * (line_strengths[OI][TRANSITION_0_to_3] +
                                line_strengths[OI][TRANSITION_1_to_3]));

    // OII
    // Osterbrock & Ferland (2006), table 3.13
    // ground state: 4S3/2
    // excited states: 2D5/2, 2D3/2, 2P3/2, 2P1/2
    // this is the sum of the 3726.0 and 3728.8 angstrom transitions
    // note that Kenny's version wrongly included the 497.1 um transition as
    // well...
    eval.set_emissivity(EMISSIONLINE_OII_3727,
                        ntot * (line_strengths[OII][TRANSITION_0_to_1] +
                                line_strengths[OII][TRANSITION_0_to_2]));
    // this is the sum of the transitions at 7319.9, 7330.7, 7318.8 and 7329.6
    // angstrom
    eval.set_emissivity(EMISSIONLINE_OII_7325,
                        ntot * (line_strengths[OII][TRANSITION_1_to_4] +
                                line_strengths[OII][TRANSITION_2_to_4] +
                                line_strengths[OII][TRANSITION_1_to_3] +
                                line_strengths[OII][TRANSITION_2_to_3]));

    // OIII
    // Osterbrock & Ferland (2006), table 3.12
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    eval.set_emissivity(EMISSIONLINE_OIII_4363,
                        ntot * line_strengths[OIII][TRANSITION_3_to_4]);
    eval.set_emissivity(EMISSIONLINE_OIII_4959,
                        ntot * line_strengths[OIII][TRANSITION_1_to_3]);
    eval.set_emissivity(EMISSIONLINE_OIII_5007,
                        ntot * line_strengths[OIII][TRANSITION_2_to_3]);
    eval.set_emissivity(EMISSIONLINE_OIII_52mu,
                        ntot * line_strengths[OIII][TRANSITION_1_to_2]);
    eval.set_emissivity(EMISSIONLINE_OIII_88mu,
                        ntot * line_strengths[OIII][TRANSITION_0_to_1]);

    // NeIII
    // Osterbrock & Ferland (2006), table 3.14
    // ground state: 3P2
    // excited states: 3P1, 3P0, 1D2, 1S0
    eval.set_emissivity(EMISSIONLINE_NeIII_3869,
                        ntot * line_strengths[NeIII][TRANSITION_0_to_3]);
    eval.set_emissivity(EMISSIONLINE_NeIII_3968,
                        ntot * line_strengths[NeIII][TRANSITION_1_to_3]);
    eval.set_emissivity(EMISSIONLINE_NeIII_15mu,
                        ntot * line_strengths[NeIII][TRANSITION_0_to_1]);

    // SII
    // Osterbrock & Ferland (2006), table 3.13
    // ground state: 4S3/2
    // excited states: 2D3/2, 2D5/2, 2P1/2, 2P3/2
    // this is the sum of the 4068.6 and 4076.4 angstrom transitions
    eval.set_emissivity(EMISSIONLINE_SII_4072,
                        ntot * (line_strengths[SII][TRANSITION_0_to_3] +
                                line_strengths[SII][TRANSITION_0_to_4]));
    // this is the sum of the 6716.5 and 6730.8 angstrom transitions
    // note that Kenny's version wrongly includes the 314.5 um transition...
    eval.set_emissivity(EMISSIONLINE_SII_6725,
                        ntot * (line_strengths[SII][TRANSITION_0_to_1] +
                                line_strengths[SII][TRANSITION_0_to_2]));

    // SIII
    // Osterbrock & Ferland (2006), table 3.12
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    // this is the sum of the 9531.0 and 9068.9 angstrom transitions
    eval.set_emissivity(EMISSIONLINE_SIII_9405,
                        ntot * (line_strengths[SIII][TRANSITION_1_to_3] +
                                line_strengths[SIII][TRANSITION_2_to_3]));
    eval.set_emissivity(EMISSIONLINE_SIII_6312,
                        ntot * line_strengths[SIII][TRANSITION_3_to_4]);
    eval.set_emissivity(EMISSIONLINE_SIII_19mu,
                        ntot * line_strengths[SIII][TRANSITION_1_to_2]);
    eval.set_emissivity(EMISSIONLINE_SIII_33mu,
                        ntot * line_strengths[SIII][TRANSITION_0_to_1]);

    // CII
    // Osterbrock & Ferland (2006), table 3.9
    // ground state: 2P1/2
    // excited states: 2P3/2, 4P1/2, 4P3/2, 4P5/2
    // this should be the sum of all 4P to 2P transitions
    // note that Kenny's code wrongly includes some 4P to 4P transitions...
    eval.set_emissivity(EMISSIONLINE_CII_2325,
                        ntot * (line_strengths[CII][TRANSITION_0_to_2] +
                                line_strengths[CII][TRANSITION_1_to_2] +
                                line_strengths[CII][TRANSITION_0_to_3] +
                                line_strengths[CII][TRANSITION_1_to_3] +
                                line_strengths[CII][TRANSITION_0_to_4] +
                                line_strengths[CII][TRANSITION_1_to_4]));

    // CIII
    // Osterbrock & Ferland (2006), table 3.8
    // ground state: 1S0
    // excited states: 3P0, 3P1, 3P2, 1P1
    // this is the sum of all 3P to 1S transitions
    // note that Kenny's code wrongly includes some 3P to 3P transitions...
    eval.set_emissivity(EMISSIONLINE_CIII_1908,
                        ntot * (line_strengths[CIII][TRANSITION_0_to_1] +
                                line_strengths[CIII][TRANSITION_0_to_2] +
                                line_strengths[CIII][TRANSITION_0_to_3]));

    // NIII
    // Osterbrock & Ferland (2006), table 3.9
    // ground state: 2P1/2
    // excited state: 2P3/2
    eval.set_emissivity(EMISSIONLINE_NIII_57mu, ntot * line_strengths[NIII][0]);

    // NeII
    // Osterbrock & Ferland (2006), table 3.11
    // ground state: 2P3/2
    // excited state: 2P1/2
    eval.set_emissivity(EMISSIONLINE_NeII_12mu, ntot * line_strengths[NeII][0]);

    // SIV
    // Osterbrock & Ferland (2006), table 3.10
    // ground state: 2P1/2
    // excited state: 2P3/2
    eval.set_emissivity(EMISSIONLINE_SIV_10mu, ntot * line_strengths[SIV][0]);

    // density weighted average temperature of ionized particles
    eval.set_emissivity(EMISSIONLINE_avg_T, ne * nhp * T);
    eval.set_emissivity(EMISSIONLINE_avg_T_count, ne * nhp);
#ifdef HAS_HELIUM
    // average ionized hydrogen and helium density product
    eval.set_emissivity(
        EMISSIONLINE_avg_nH_nHe,
        ne * (1. - ionization_variables.get_ionic_fraction(ION_He_n)));
#endif
    eval.set_emissivity(
        EMISSIONLINE_avg_nH_nHe_count,
        ne * (1. - ionization_variables.get_ionic_fraction(ION_H_n)));
    // we converted Kenny's constant from 1.e20 erg/cm^6/s to J/m^6/s
    // Osterbrock & Ferland (2006), table 4.6 (fit?)
    eval.set_emissivity(EMISSIONLINE_HeI_5876,
                        ne * nhep * 1.69e-38 * std::pow(T4, -1.065));
    // Verner & Ferland (1996), table 1
    eval.set_emissivity(
        EMISSIONLINE_Hrec_s,
        ne * nhp * 7.982e-23 /
            (std::sqrt(T / 3.148) * std::pow(1. + std::sqrt(T / 3.148), 0.252) *
             std::pow(1. + std::sqrt(T / 7.036e5), 1.748)));

    // HST WFC2 filters
    // wavelength ranges were based on data found on
    // http://www.stsci.edu/hst/wfpc2/documents/handbook/IHB_17.html
    // F439W: effective wavelength: 4283 A, width: 464.4 A
    //        range: [4051 A, 4515 A]
    // lines that contribute: OIII: 4 --> 3
    //                        SIII: 4 --> 0, 3 --> 0
    eval.set_emissivity(EMISSIONLINE_WFC2_F439W,
                        ntot * (line_strengths[OIII][TRANSITION_3_to_4] +
                                line_strengths[SIII][TRANSITION_0_to_3] +
                                line_strengths[SIII][TRANSITION_0_to_4]));
    // F555W: effective wavelength: 5202 A, width: 1222.6 A
    //        range: [4591 A, 5813 A]
    // lines that contribute: NI: 2 --> 0, 1 --> 0
    //                        NII: 4 --> 3
    //                        OI: 4 --> 3
    //                        OIII: 3 --> 2, 3 --> 1, 3 --> 0
    //                        HBeta
    eval.set_emissivity(EMISSIONLINE_WFC2_F555W,
                        eval.get_emissivity(EMISSIONLINE_HBeta) +
                            ntot * (line_strengths[NI][TRANSITION_0_to_1] +
                                    line_strengths[NI][TRANSITION_0_to_2] +
                                    line_strengths[NII][TRANSITION_3_to_4] +
                                    line_strengths[OI][TRANSITION_3_to_4] +
                                    line_strengths[OIII][TRANSITION_0_to_3] +
                                    line_strengths[OIII][TRANSITION_1_to_3] +
                                    line_strengths[OIII][TRANSITION_2_to_3]));
    // F675W: effective wavelength: 6714 A, width: 889.5 A
    //        range: [6269 A, 7159 A]
    // lines that contribute: NII: 3 --> 2, 3 --> 1, 3 --> 0
    //                        OI: 3 --> 2, 3 --> 1, 3 --> 0
    //                        SII: 2 --> 0, 1 --> 0
    //                        SIII: 4 --> 3
    //                        HAlpha
    eval.set_emissivity(EMISSIONLINE_WFC2_F675W,
                        eval.get_emissivity(EMISSIONLINE_HAlpha) +
                            ntot * (line_strengths[NII][TRANSITION_0_to_3] +
                                    line_strengths[NII][TRANSITION_1_to_3] +
                                    line_strengths[NII][TRANSITION_2_to_3] +
                                    line_strengths[OI][TRANSITION_0_to_3] +
                                    line_strengths[OI][TRANSITION_1_to_3] +
                                    line_strengths[OI][TRANSITION_2_to_3] +
                                    line_strengths[SII][TRANSITION_0_to_1] +
                                    line_strengths[SII][TRANSITION_0_to_2] +
                                    line_strengths[SIII][TRANSITION_3_to_4]));
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
    EmissivityValues *emissivities =
        new EmissivityValues(calculate_emissivities(
            it.get_ionization_variables(), _abundances, _lines));
    it.set_emissivities(emissivities);
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
  size_t index = 0;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    result[index] = calculate_emissivities(it.get_ionization_variables(),
                                           _abundances, _lines);
    ++index;
  }
  return result;
}

/**
 * @brief Get the emissivities for the given cell.
 *
 * @param ionization_variables Cell variables.
 * @param do_line Array telling which lines to output.
 * @param output Array to store the output in.
 */
void EmissivityCalculator::calculate_emissivities(
    const IonizationVariables &ionization_variables,
    const bool do_line[NUMBER_OF_EMISSIONLINES],
    double output[NUMBER_OF_EMISSIONLINES]) const {

  EmissivityValues values =
      calculate_emissivities(ionization_variables, _abundances, _lines);
  for (int_fast32_t line = 0; line < NUMBER_OF_EMISSIONLINES; ++line) {
    if (do_line[line]) {
      output[line] = values.get_emissivity(line);
    }
  }
}
