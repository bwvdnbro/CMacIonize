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
 * @file LineCoolingData.cpp
 *
 * @brief Internal representation of the external data used for line cooling:
 * implementation
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "LineCoolingData.hpp"
#include "Error.hpp"
#include "PhysicalConstants.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

/**
 * @brief Constructor.
 *
 * Initializes the data values.
 */
LineCoolingData::LineCoolingData() {

  /// some energy conversion constants (for energy differences)

  // conversion factor from cm^-1 to K
  const double hc_over_k =
      100. * PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK) *
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED) /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

  // conversion factor from eV to K
  const double eV_over_k =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_ELECTRONVOLT) /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

  // conversion factor from Ry to K
  const double Ry_over_k =
      PhysicalConstants::get_physical_constant(
          PHYSICALCONSTANT_RYDBERG_ENERGY) /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

  /// five level elements

  /// NI
  {
    // data from Froese Fischer & Tachiev (2004), table 4
    // ground state: 4S3/2
    // excited states: 2D5/2, 2D3/2, 2P1/2, 2P3/2
    // in cm^-1
    const double energy_levels[4] = {19224.37, 19233.09, 28839.05, 28839.07};
    _inverse_statistical_weight[NI][0] = 0.25;
    _inverse_statistical_weight[NI][1] = 1. / 6.;
    _inverse_statistical_weight[NI][2] = 0.25;
    _inverse_statistical_weight[NI][3] = 0.5;
    _inverse_statistical_weight[NI][4] = 0.25;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[NI][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[NI][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[NI][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[NI][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[NI][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[NI][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[NI][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[NI][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[NI][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[NI][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we take the sum of all transition probabilities (in s^-1)
    _transition_probability[NI][TRANSITION_0_to_1] = 7.566e-6;
    _transition_probability[NI][TRANSITION_0_to_2] = 2.029e-5;
    _transition_probability[NI][TRANSITION_0_to_3] = 2.605e-3;
    _transition_probability[NI][TRANSITION_0_to_4] = 6.504e-3;
    _transition_probability[NI][TRANSITION_1_to_2] = 1.071e-8;
    _transition_probability[NI][TRANSITION_1_to_3] = 3.447e-2;
    _transition_probability[NI][TRANSITION_1_to_4] = 6.135e-2;
    _transition_probability[NI][TRANSITION_2_to_3] = 5.272e-2;
    _transition_probability[NI][TRANSITION_2_to_4] = 2.745e-2;
    _transition_probability[NI][TRANSITION_3_to_4] = 9.504e-17;
    // our own fits to the data of Tayal (2000)
    // these fits were made with the script data/linecooling/gamma_NI.py
    // (this script also outputs the code below)
    _collision_strength[NI][TRANSITION_0_to_1][0] = -0.352170660298;
    _collision_strength[NI][TRANSITION_0_to_1][1] = 290.570335282;
    _collision_strength[NI][TRANSITION_0_to_1][2] = -88972.2171524;
    _collision_strength[NI][TRANSITION_0_to_1][3] = -5864.47392717;
    _collision_strength[NI][TRANSITION_0_to_1][4] = -34.6307594052;
    _collision_strength[NI][TRANSITION_0_to_1][5] = 5864.48664195;
    _collision_strength[NI][TRANSITION_0_to_1][6] = 0.999999834749;
    _collision_strength[NI][TRANSITION_0_to_2][0] = 1.04662088312;
    _collision_strength[NI][TRANSITION_0_to_2][1] = -0.000192854143976;
    _collision_strength[NI][TRANSITION_0_to_2][2] = 0.0933869336879;
    _collision_strength[NI][TRANSITION_0_to_2][3] = 62.4700577714;
    _collision_strength[NI][TRANSITION_0_to_2][4] = 2.12972369491e-05;
    _collision_strength[NI][TRANSITION_0_to_2][5] = -62.4700577752;
    _collision_strength[NI][TRANSITION_0_to_2][6] = 0.999999999995;
    _collision_strength[NI][TRANSITION_0_to_3][0] = -0.143686075261;
    _collision_strength[NI][TRANSITION_0_to_3][1] = 4.85458860278;
    _collision_strength[NI][TRANSITION_0_to_3][2] = -1635.9084968;
    _collision_strength[NI][TRANSITION_0_to_3][3] = -85.5223871776;
    _collision_strength[NI][TRANSITION_0_to_3][4] = -0.570897281076;
    _collision_strength[NI][TRANSITION_0_to_3][5] = 85.5225946457;
    _collision_strength[NI][TRANSITION_0_to_3][6] = 0.999999813378;
    _collision_strength[NI][TRANSITION_0_to_4][0] = -0.146192751421;
    _collision_strength[NI][TRANSITION_0_to_4][1] = 10.0505771787;
    _collision_strength[NI][TRANSITION_0_to_4][2] = -3382.89322868;
    _collision_strength[NI][TRANSITION_0_to_4][3] = 3474.60389784;
    _collision_strength[NI][TRANSITION_0_to_4][4] = -1.1821988914;
    _collision_strength[NI][TRANSITION_0_to_4][5] = -3474.60346722;
    _collision_strength[NI][TRANSITION_0_to_4][6] = 1.00000000954;
    _collision_strength[NI][TRANSITION_1_to_2][0] = 1.08380124933;
    _collision_strength[NI][TRANSITION_1_to_2][1] = -0.00824248994995;
    _collision_strength[NI][TRANSITION_1_to_2][2] = 6.09352690919;
    _collision_strength[NI][TRANSITION_1_to_2][3] = 5.12778710603e-05;
    _collision_strength[NI][TRANSITION_1_to_2][4] = 0.000897416464378;
    _collision_strength[NI][TRANSITION_1_to_2][5] = -5.14370296385e-05;
    _collision_strength[NI][TRANSITION_1_to_2][6] = 0.999765473624;
    _collision_strength[NI][TRANSITION_1_to_3][0] = 0.825614395297;
    _collision_strength[NI][TRANSITION_1_to_3][1] = -0.0452768551043;
    _collision_strength[NI][TRANSITION_1_to_3][2] = 19.6917639859;
    _collision_strength[NI][TRANSITION_1_to_3][3] = -38.3932561041;
    _collision_strength[NI][TRANSITION_1_to_3][4] = 0.00508048299418;
    _collision_strength[NI][TRANSITION_1_to_3][5] = 38.3932550214;
    _collision_strength[NI][TRANSITION_1_to_3][6] = 1.00000000216;
    _collision_strength[NI][TRANSITION_1_to_4][0] = 0.673408607924;
    _collision_strength[NI][TRANSITION_1_to_4][1] = -1.30037729225;
    _collision_strength[NI][TRANSITION_1_to_4][2] = 549.88707741;
    _collision_strength[NI][TRANSITION_1_to_4][3] = -0.0687375130436;
    _collision_strength[NI][TRANSITION_1_to_4][4] = 0.146264389311;
    _collision_strength[NI][TRANSITION_1_to_4][5] = 0.0687059692896;
    _collision_strength[NI][TRANSITION_1_to_4][6] = 1.00003514965;
    _collision_strength[NI][TRANSITION_2_to_3][0] = 0.639898658243;
    _collision_strength[NI][TRANSITION_2_to_3][1] = -0.913913781086;
    _collision_strength[NI][TRANSITION_2_to_3][2] = 385.259711068;
    _collision_strength[NI][TRANSITION_2_to_3][3] = -49.3897876926;
    _collision_strength[NI][TRANSITION_2_to_3][4] = 0.102833403715;
    _collision_strength[NI][TRANSITION_2_to_3][5] = 49.3897653939;
    _collision_strength[NI][TRANSITION_2_to_3][6] = 1.00000003459;
    _collision_strength[NI][TRANSITION_2_to_4][0] = 0.779047782547;
    _collision_strength[NI][TRANSITION_2_to_4][1] = -0.129384116204;
    _collision_strength[NI][TRANSITION_2_to_4][2] = 55.57362588;
    _collision_strength[NI][TRANSITION_2_to_4][3] = -0.244902789134;
    _collision_strength[NI][TRANSITION_2_to_4][4] = 0.0145342075843;
    _collision_strength[NI][TRANSITION_2_to_4][5] = 0.244899675501;
    _collision_strength[NI][TRANSITION_2_to_4][6] = 1.00000097316;
    _collision_strength[NI][TRANSITION_3_to_4][0] = 0.979889001396;
    _collision_strength[NI][TRANSITION_3_to_4][1] = -0.0250254133594;
    _collision_strength[NI][TRANSITION_3_to_4][2] = 13.3611860055;
    _collision_strength[NI][TRANSITION_3_to_4][3] = 120.384831511;
    _collision_strength[NI][TRANSITION_3_to_4][4] = 0.00278533555664;
    _collision_strength[NI][TRANSITION_3_to_4][5] = -120.384832094;
    _collision_strength[NI][TRANSITION_3_to_4][6] = 0.999999999629;
  }

  /// NII
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 4 and 5
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {48.7, 130.8, 15316.3, 32688.9};
    _inverse_statistical_weight[NII][0] = 1.;
    _inverse_statistical_weight[NII][1] = 1. / 3.;
    _inverse_statistical_weight[NII][2] = 0.2;
    _inverse_statistical_weight[NII][3] = 0.2;
    _inverse_statistical_weight[NII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[NII][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[NII][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[NII][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[NII][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[NII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[NII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[NII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[NII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[NII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[NII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _transition_probability[NII][TRANSITION_0_to_1] = 2.077e-6;
    _transition_probability[NII][TRANSITION_0_to_2] = 1.127e-12;
    _transition_probability[NII][TRANSITION_0_to_3] = 3.554e-7;
    _transition_probability[NII][TRANSITION_0_to_4] = 0.;
    _transition_probability[NII][TRANSITION_1_to_2] = 7.463e-6;
    _transition_probability[NII][TRANSITION_1_to_3] = 1.016e-3;
    _transition_probability[NII][TRANSITION_1_to_4] = 3.297e-2;
    _transition_probability[NII][TRANSITION_2_to_3] = 3.005e-3;
    _transition_probability[NII][TRANSITION_2_to_4] = 1.315e-4;
    _transition_probability[NII][TRANSITION_3_to_4] = 1.023;
    // our own fits to the data of Lennon & Burke (1994)
    // these fits were made with the script data/linecooling/gamma_NII.py
    // (this script also outputs the code below)
    _collision_strength[NII][TRANSITION_0_to_1][0] = 0.651709467718;
    _collision_strength[NII][TRANSITION_0_to_1][1] = 0.0043840698298;
    _collision_strength[NII][TRANSITION_0_to_1][2] = 2.05865457437;
    _collision_strength[NII][TRANSITION_0_to_1][3] = -0.000183596343291;
    _collision_strength[NII][TRANSITION_0_to_1][4] = -0.000403172016047;
    _collision_strength[NII][TRANSITION_0_to_1][5] = 0.000183636758783;
    _collision_strength[NII][TRANSITION_0_to_1][6] = 0.999983437636;
    _collision_strength[NII][TRANSITION_0_to_2][0] = 0.99721996664;
    _collision_strength[NII][TRANSITION_0_to_2][1] = -2.63996934679e-06;
    _collision_strength[NII][TRANSITION_0_to_2][2] = 0.212680283457;
    _collision_strength[NII][TRANSITION_0_to_2][3] = -0.000299568497993;
    _collision_strength[NII][TRANSITION_0_to_2][4] = 1.36812200274e-06;
    _collision_strength[NII][TRANSITION_0_to_2][5] = 0.000299567319572;
    _collision_strength[NII][TRANSITION_0_to_2][6] = 1.00000030796;
    _collision_strength[NII][TRANSITION_0_to_3][0] = 0.985274169209;
    _collision_strength[NII][TRANSITION_0_to_3][1] = 2.66502677175e-05;
    _collision_strength[NII][TRANSITION_0_to_3][2] = 0.296461402639;
    _collision_strength[NII][TRANSITION_0_to_3][3] = 0.00018880364217;
    _collision_strength[NII][TRANSITION_0_to_3][4] = -2.56029417716e-06;
    _collision_strength[NII][TRANSITION_0_to_3][5] = -0.000188803354544;
    _collision_strength[NII][TRANSITION_0_to_3][6] = 1.00000011499;
    _collision_strength[NII][TRANSITION_0_to_4][0] = 1.00242568631;
    _collision_strength[NII][TRANSITION_0_to_4][1] = -1.16975224507e-07;
    _collision_strength[NII][TRANSITION_0_to_4][2] = 0.0306575644987;
    _collision_strength[NII][TRANSITION_0_to_4][3] = -0.000191065656783;
    _collision_strength[NII][TRANSITION_0_to_4][4] = 2.80204262642e-08;
    _collision_strength[NII][TRANSITION_0_to_4][5] = 0.000191065650186;
    _collision_strength[NII][TRANSITION_0_to_4][6] = 1.0000000024;
    _collision_strength[NII][TRANSITION_1_to_2][0] = 0.898544649079;
    _collision_strength[NII][TRANSITION_1_to_2][1] = 0.000626320672746;
    _collision_strength[NII][TRANSITION_1_to_2][2] = 1.55901247212;
    _collision_strength[NII][TRANSITION_1_to_2][3] = -3.57105490335e-07;
    _collision_strength[NII][TRANSITION_1_to_2][4] = -5.44511032111e-05;
    _collision_strength[NII][TRANSITION_1_to_2][5] = 3.57952429654e-07;
    _collision_strength[NII][TRANSITION_1_to_2][6] = 0.99987288108;
    _collision_strength[NII][TRANSITION_1_to_3][0] = 0.98527416532;
    _collision_strength[NII][TRANSITION_1_to_3][1] = 7.99508259627e-05;
    _collision_strength[NII][TRANSITION_1_to_3][2] = 0.88938422702;
    _collision_strength[NII][TRANSITION_1_to_3][3] = -0.000187717979801;
    _collision_strength[NII][TRANSITION_1_to_3][4] = -7.68088510971e-06;
    _collision_strength[NII][TRANSITION_1_to_3][5] = 0.000187718842681;
    _collision_strength[NII][TRANSITION_1_to_3][6] = 0.999999653029;
    _collision_strength[NII][TRANSITION_1_to_4][0] = 1.00242568648;
    _collision_strength[NII][TRANSITION_1_to_4][1] = -3.50925761487e-07;
    _collision_strength[NII][TRANSITION_1_to_4][2] = 0.0919726934133;
    _collision_strength[NII][TRANSITION_1_to_4][3] = 0.000184955818388;
    _collision_strength[NII][TRANSITION_1_to_4][4] = 8.40612886609e-08;
    _collision_strength[NII][TRANSITION_1_to_4][5] = -0.000184955838179;
    _collision_strength[NII][TRANSITION_1_to_4][6] = 0.999999992557;
    _collision_strength[NII][TRANSITION_2_to_3][0] = 0.985274172887;
    _collision_strength[NII][TRANSITION_2_to_3][1] = 0.000133251302615;
    _collision_strength[NII][TRANSITION_2_to_3][2] = 1.48230698308;
    _collision_strength[NII][TRANSITION_2_to_3][3] = 0.000192066237152;
    _collision_strength[NII][TRANSITION_2_to_3][4] = -1.28014668146e-05;
    _collision_strength[NII][TRANSITION_2_to_3][5] = -0.000192064799027;
    _collision_strength[NII][TRANSITION_2_to_3][6] = 1.00000056519;
    _collision_strength[NII][TRANSITION_2_to_4][0] = 1.00242568736;
    _collision_strength[NII][TRANSITION_2_to_4][1] = -5.84876871304e-07;
    _collision_strength[NII][TRANSITION_2_to_4][2] = 0.153287821582;
    _collision_strength[NII][TRANSITION_2_to_4][3] = -0.000189310275131;
    _collision_strength[NII][TRANSITION_2_to_4][4] = 1.40102211883e-07;
    _collision_strength[NII][TRANSITION_2_to_4][5] = 0.000189310242146;
    _collision_strength[NII][TRANSITION_2_to_4][6] = 1.00000001212;
    _collision_strength[NII][TRANSITION_3_to_4][0] = 1.02226099745;
    _collision_strength[NII][TRANSITION_3_to_4][1] = -0.000391814027158;
    _collision_strength[NII][TRANSITION_3_to_4][2] = 1.12144974144;
    _collision_strength[NII][TRANSITION_3_to_4][3] = -5.21675381019e-05;
    _collision_strength[NII][TRANSITION_3_to_4][4] = 3.93432063927e-05;
    _collision_strength[NII][TRANSITION_3_to_4][5] = 5.21628733234e-05;
    _collision_strength[NII][TRANSITION_3_to_4][6] = 1.00000664715;
  }

  /// OI
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 6 and 7
    // ground state: 3P2
    // excited states: 3P1, 3P0, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {158., 227., 15868., 33793.};
    _inverse_statistical_weight[OI][0] = 0.2;
    _inverse_statistical_weight[OI][1] = 1. / 3.;
    _inverse_statistical_weight[OI][2] = 1.;
    _inverse_statistical_weight[OI][3] = 0.2;
    _inverse_statistical_weight[OI][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[OI][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[OI][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[OI][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[OI][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[OI][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[OI][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[OI][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[OI][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[OI][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[OI][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _transition_probability[OI][TRANSITION_0_to_1] = 8.865e-5;
    _transition_probability[OI][TRANSITION_0_to_2] = 1.275e-10;
    _transition_probability[OI][TRANSITION_0_to_3] = 6.535e-3;
    _transition_probability[OI][TRANSITION_0_to_4] = 2.945e-4;
    _transition_probability[OI][TRANSITION_1_to_2] = 1.772e-5;
    _transition_probability[OI][TRANSITION_1_to_3] = 2.111e-3;
    _transition_probability[OI][TRANSITION_1_to_4] = 7.909e-2;
    _transition_probability[OI][TRANSITION_2_to_3] = 6.388e-7;
    _transition_probability[OI][TRANSITION_2_to_4] = 0.;
    _transition_probability[OI][TRANSITION_3_to_4] = 1.124;
    // our own fits to the data of Zatsarinny & Tayal (2003) and Berrington
    // (1988)
    // these fits were made with the script data/linecooling/gamma_OI.py
    // (this script also outputs the code below)
    // there is only data available for fine structure transitions below
    // 10,000 K, so we need to extrapolate for higher temperatures
    _collision_strength[OI][TRANSITION_0_to_1][0] = 0.531891968668;
    _collision_strength[OI][TRANSITION_0_to_1][1] = 0.000174043443031;
    _collision_strength[OI][TRANSITION_0_to_1][2] = -0.00289046667195;
    _collision_strength[OI][TRANSITION_0_to_1][3] = 6.05503334723e-08;
    _collision_strength[OI][TRANSITION_0_to_1][4] = 1.19694301884e-06;
    _collision_strength[OI][TRANSITION_0_to_1][5] = 311.969799446;
    _collision_strength[OI][TRANSITION_0_to_1][6] = -430.811510152;
    _collision_strength[OI][TRANSITION_0_to_2][0] = 0.709634874012;
    _collision_strength[OI][TRANSITION_0_to_2][1] = 4.19344490648e-05;
    _collision_strength[OI][TRANSITION_0_to_2][2] = -9.88929290533e-05;
    _collision_strength[OI][TRANSITION_0_to_2][3] = 2.34148281629e-09;
    _collision_strength[OI][TRANSITION_0_to_2][4] = -2.03916544987e-06;
    _collision_strength[OI][TRANSITION_0_to_2][5] = 0.450007261425;
    _collision_strength[OI][TRANSITION_0_to_2][6] = -522.886828465;
    _collision_strength[OI][TRANSITION_0_to_3][0] = 0.577956016286;
    _collision_strength[OI][TRANSITION_0_to_3][1] = -0.00328947147455;
    _collision_strength[OI][TRANSITION_0_to_3][2] = 0.552970782258;
    _collision_strength[OI][TRANSITION_0_to_3][3] = -1.0868651959e-08;
    _collision_strength[OI][TRANSITION_0_to_3][4] = 0.000449847254588;
    _collision_strength[OI][TRANSITION_0_to_3][5] = -3112.83302323;
    _collision_strength[OI][TRANSITION_0_to_3][6] = -426.886279035;
    _collision_strength[OI][TRANSITION_0_to_4][0] = 0.447157138186;
    _collision_strength[OI][TRANSITION_0_to_4][1] = -0.00199786430506;
    _collision_strength[OI][TRANSITION_0_to_4][2] = 0.367624253613;
    _collision_strength[OI][TRANSITION_0_to_4][3] = -4.72921770139e-09;
    _collision_strength[OI][TRANSITION_0_to_4][4] = 0.000249849906874;
    _collision_strength[OI][TRANSITION_0_to_4][5] = 362.968742149;
    _collision_strength[OI][TRANSITION_0_to_4][6] = -2159.85840959;
    _collision_strength[OI][TRANSITION_1_to_2][0] = 0.445351360073;
    _collision_strength[OI][TRANSITION_1_to_2][1] = 0.00019090023886;
    _collision_strength[OI][TRANSITION_1_to_2][2] = -0.00345127570802;
    _collision_strength[OI][TRANSITION_1_to_2][3] = 5.71292820253e-08;
    _collision_strength[OI][TRANSITION_1_to_2][4] = -3.18737719028e-05;
    _collision_strength[OI][TRANSITION_1_to_2][5] = 394.278912444;
    _collision_strength[OI][TRANSITION_1_to_2][6] = -19.4605389548;
    _collision_strength[OI][TRANSITION_1_to_3][0] = 0.577956016302;
    _collision_strength[OI][TRANSITION_1_to_3][1] = -0.00197368288433;
    _collision_strength[OI][TRANSITION_1_to_3][2] = 0.331782469279;
    _collision_strength[OI][TRANSITION_1_to_3][3] = -6.5211911744e-09;
    _collision_strength[OI][TRANSITION_1_to_3][4] = 0.000269908352701;
    _collision_strength[OI][TRANSITION_1_to_3][5] = -3112.83302714;
    _collision_strength[OI][TRANSITION_1_to_3][6] = -426.886279331;
    _collision_strength[OI][TRANSITION_1_to_4][0] = 0.447157138318;
    _collision_strength[OI][TRANSITION_1_to_4][1] = -0.00119871858132;
    _collision_strength[OI][TRANSITION_1_to_4][2] = 0.220574551826;
    _collision_strength[OI][TRANSITION_1_to_4][3] = -2.83753061784e-09;
    _collision_strength[OI][TRANSITION_1_to_4][4] = 0.000149909943915;
    _collision_strength[OI][TRANSITION_1_to_4][5] = 364.603929362;
    _collision_strength[OI][TRANSITION_1_to_4][6] = -2170.49654928;
    _collision_strength[OI][TRANSITION_2_to_3][0] = 0.577956016277;
    _collision_strength[OI][TRANSITION_2_to_3][1] = -0.000657894294987;
    _collision_strength[OI][TRANSITION_2_to_3][2] = 0.110594156466;
    _collision_strength[OI][TRANSITION_2_to_3][3] = -2.17373039198e-09;
    _collision_strength[OI][TRANSITION_2_to_3][4] = 8.99694509274e-05;
    _collision_strength[OI][TRANSITION_2_to_3][5] = -3112.83303105;
    _collision_strength[OI][TRANSITION_2_to_3][6] = -426.886279628;
    _collision_strength[OI][TRANSITION_2_to_4][0] = 0.44715713653;
    _collision_strength[OI][TRANSITION_2_to_4][1] = -0.000399572868221;
    _collision_strength[OI][TRANSITION_2_to_4][2] = 0.0735248521538;
    _collision_strength[OI][TRANSITION_2_to_4][3] = -9.45843552833e-10;
    _collision_strength[OI][TRANSITION_2_to_4][4] = 4.99699822537e-05;
    _collision_strength[OI][TRANSITION_2_to_4][5] = 286.606313767;
    _collision_strength[OI][TRANSITION_2_to_4][6] = -1714.06509463;
    _collision_strength[OI][TRANSITION_3_to_4][0] = 0.77784636923;
    _collision_strength[OI][TRANSITION_3_to_4][1] = 0.000256767786556;
    _collision_strength[OI][TRANSITION_3_to_4][2] = -0.0770418161828;
    _collision_strength[OI][TRANSITION_3_to_4][3] = 1.62454220113e-10;
    _collision_strength[OI][TRANSITION_3_to_4][4] = -1.98264637506e-05;
    _collision_strength[OI][TRANSITION_3_to_4][5] = 70.2374763213;
    _collision_strength[OI][TRANSITION_3_to_4][6] = -1707.98593908;
  }

  /// OII
  {
    // data from Froese Fischer & Tachiev (2004), table 4
    // ground state: 4S3/2
    // excited states: 2D5/2, 2D3/2, 2P3/2, 2P1/2
    // in cm^-1
    const double energy_levels[4] = {26810.73, 26830.45, 40468.36, 40470.96};
    _inverse_statistical_weight[OII][0] = 0.25;
    _inverse_statistical_weight[OII][1] = 1. / 6.;
    _inverse_statistical_weight[OII][2] = 0.25;
    _inverse_statistical_weight[OII][3] = 0.25;
    _inverse_statistical_weight[OII][4] = 0.5;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[OII][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[OII][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[OII][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[OII][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[OII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[OII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[OII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[OII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[OII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[OII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we sum contributions of all types (in s^-1)
    _transition_probability[OII][TRANSITION_0_to_1] = 4.124e-5;
    _transition_probability[OII][TRANSITION_0_to_2] = 1.635e-4;
    _transition_probability[OII][TRANSITION_0_to_3] = 5.646e-2;
    _transition_probability[OII][TRANSITION_0_to_4] = 2.265e-2;
    _transition_probability[OII][TRANSITION_1_to_2] = 1.241e-7;
    _transition_probability[OII][TRANSITION_1_to_3] = 1.106e-1;
    _transition_probability[OII][TRANSITION_1_to_4] = 5.824e-2;
    _transition_probability[OII][TRANSITION_2_to_3] = 5.871e-2;
    _transition_probability[OII][TRANSITION_2_to_4] = 9.668e-2;
    _transition_probability[OII][TRANSITION_3_to_4] = 3.158e-10;
    // our own fits to the data of Kisielius et al. (2009)
    // these fits were made with the script data/linecooling/gamma_OII.py
    // (this script also outputs the code below)
    _collision_strength[OII][TRANSITION_0_to_1][0] = 0.997482865693;
    _collision_strength[OII][TRANSITION_0_to_1][1] = -6.25129361407e-06;
    _collision_strength[OII][TRANSITION_0_to_1][2] = 0.85013614096;
    _collision_strength[OII][TRANSITION_0_to_1][3] = -0.000192875793831;
    _collision_strength[OII][TRANSITION_0_to_1][4] = 7.24338376045e-07;
    _collision_strength[OII][TRANSITION_0_to_1][5] = 0.000192875807718;
    _collision_strength[OII][TRANSITION_0_to_1][6] = 0.99999998866;
    _collision_strength[OII][TRANSITION_0_to_2][0] = 0.988746408499;
    _collision_strength[OII][TRANSITION_0_to_2][1] = -8.42756442139e-07;
    _collision_strength[OII][TRANSITION_0_to_2][2] = 0.603631681835;
    _collision_strength[OII][TRANSITION_0_to_2][3] = -0.000105058004912;
    _collision_strength[OII][TRANSITION_0_to_2][4] = 2.17871753378e-07;
    _collision_strength[OII][TRANSITION_0_to_2][5] = 0.000105057997886;
    _collision_strength[OII][TRANSITION_0_to_2][6] = 0.999999999636;
    _collision_strength[OII][TRANSITION_0_to_3][0] = 1.01835073138;
    _collision_strength[OII][TRANSITION_0_to_3][1] = -8.27164528694e-06;
    _collision_strength[OII][TRANSITION_0_to_3][2] = 0.218331491341;
    _collision_strength[OII][TRANSITION_0_to_3][3] = -0.000190517152037;
    _collision_strength[OII][TRANSITION_0_to_3][4] = 9.09027104327e-07;
    _collision_strength[OII][TRANSITION_0_to_3][5] = 0.000190517050403;
    _collision_strength[OII][TRANSITION_0_to_3][6] = 1.0000000379;
    _collision_strength[OII][TRANSITION_0_to_4][0] = 1.07789473168;
    _collision_strength[OII][TRANSITION_0_to_4][1] = -1.02597843387e-05;
    _collision_strength[OII][TRANSITION_0_to_4][2] = 0.0746395510181;
    _collision_strength[OII][TRANSITION_0_to_4][3] = -0.000184748722689;
    _collision_strength[OII][TRANSITION_0_to_4][4] = 1.04586672193e-06;
    _collision_strength[OII][TRANSITION_0_to_4][5] = 0.000184748595148;
    _collision_strength[OII][TRANSITION_0_to_4][6] = 1.00000005113;
    _collision_strength[OII][TRANSITION_1_to_2][0] = 1.09674389746;
    _collision_strength[OII][TRANSITION_1_to_2][1] = -0.000170228626507;
    _collision_strength[OII][TRANSITION_1_to_2][2] = 0.687307627909;
    _collision_strength[OII][TRANSITION_1_to_2][3] = -0.000224682674283;
    _collision_strength[OII][TRANSITION_1_to_2][4] = 1.70504922064e-05;
    _collision_strength[OII][TRANSITION_1_to_2][5] = 0.000224680734991;
    _collision_strength[OII][TRANSITION_1_to_2][6] = 1.00000063846;
    _collision_strength[OII][TRANSITION_1_to_3][0] = 0.989713168687;
    _collision_strength[OII][TRANSITION_1_to_3][1] = -7.34069674831e-05;
    _collision_strength[OII][TRANSITION_1_to_3][2] = 0.884455714824;
    _collision_strength[OII][TRANSITION_1_to_3][3] = -5.55972760398e-05;
    _collision_strength[OII][TRANSITION_1_to_3][4] = 9.29108980844e-06;
    _collision_strength[OII][TRANSITION_1_to_3][5] = 5.55950017245e-05;
    _collision_strength[OII][TRANSITION_1_to_3][6] = 1.00000309792;
    _collision_strength[OII][TRANSITION_1_to_4][0] = 0.947603573757;
    _collision_strength[OII][TRANSITION_1_to_4][1] = -6.73041136806e-06;
    _collision_strength[OII][TRANSITION_1_to_4][2] = 0.466652575724;
    _collision_strength[OII][TRANSITION_1_to_4][3] = -0.000414534389992;
    _collision_strength[OII][TRANSITION_1_to_4][4] = 1.97530866409e-06;
    _collision_strength[OII][TRANSITION_1_to_4][5] = 0.000414533353556;
    _collision_strength[OII][TRANSITION_1_to_4][6] = 1.00000019139;
    _collision_strength[OII][TRANSITION_2_to_3][0] = 0.965029355403;
    _collision_strength[OII][TRANSITION_2_to_3][1] = -2.31141091327e-05;
    _collision_strength[OII][TRANSITION_2_to_3][2] = 0.57711971218;
    _collision_strength[OII][TRANSITION_2_to_3][3] = 1.72396758467e-05;
    _collision_strength[OII][TRANSITION_2_to_3][4] = 3.75814829055e-06;
    _collision_strength[OII][TRANSITION_2_to_3][5] = -1.72410045449e-05;
    _collision_strength[OII][TRANSITION_2_to_3][6] = 0.99999413037;
    _collision_strength[OII][TRANSITION_2_to_4][0] = 1.01037890101;
    _collision_strength[OII][TRANSITION_2_to_4][1] = -3.42480054524e-05;
    _collision_strength[OII][TRANSITION_2_to_4][2] = 0.299599079452;
    _collision_strength[OII][TRANSITION_2_to_4][3] = 8.65931572433e-05;
    _collision_strength[OII][TRANSITION_2_to_4][4] = 4.00278953831e-06;
    _collision_strength[OII][TRANSITION_2_to_4][5] = -8.65939756713e-05;
    _collision_strength[OII][TRANSITION_2_to_4][6] = 0.999999286977;
    _collision_strength[OII][TRANSITION_3_to_4][0] = 1.0646553085;
    _collision_strength[OII][TRANSITION_3_to_4][1] = -3.64347494013e-05;
    _collision_strength[OII][TRANSITION_3_to_4][2] = 0.183829080779;
    _collision_strength[OII][TRANSITION_3_to_4][3] = 0.000166420572534;
    _collision_strength[OII][TRANSITION_3_to_4][4] = 3.85636161175e-06;
    _collision_strength[OII][TRANSITION_3_to_4][5] = -0.000166421128586;
    _collision_strength[OII][TRANSITION_3_to_4][6] = 0.999999750898;
  }

  /// OIII
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 4 and 5
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {114., 307., 20274., 43186.};
    _inverse_statistical_weight[OIII][0] = 1.;
    _inverse_statistical_weight[OIII][1] = 1. / 3.;
    _inverse_statistical_weight[OIII][2] = 0.2;
    _inverse_statistical_weight[OIII][3] = 0.2;
    _inverse_statistical_weight[OIII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[OIII][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[OIII][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[OIII][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[OIII][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[OIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[OIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[OIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[OIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[OIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[OIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _transition_probability[OIII][TRANSITION_0_to_1] = 2.664e-5;
    _transition_probability[OIII][TRANSITION_0_to_2] = 3.094e-11;
    _transition_probability[OIII][TRANSITION_0_to_3] = 1.69e-6;
    _transition_probability[OIII][TRANSITION_0_to_4] = 0.;
    _transition_probability[OIII][TRANSITION_1_to_2] = 9.695e-5;
    _transition_probability[OIII][TRANSITION_1_to_3] = 6.995e-3;
    _transition_probability[OIII][TRANSITION_1_to_4] = 2.268e-1;
    _transition_probability[OIII][TRANSITION_2_to_3] = 2.041e-2;
    _transition_probability[OIII][TRANSITION_2_to_4] = 6.091e-4;
    _transition_probability[OIII][TRANSITION_3_to_4] = 1.561;
    // our own fits to the data of Lennon & Burke (1994)
    // these fits were made with the script data/linecooling/gamma_OIII.py
    // (this script also outputs the code below)
    _collision_strength[OIII][TRANSITION_0_to_1][0] = 0.748234291392;
    _collision_strength[OIII][TRANSITION_0_to_1][1] = 0.00232792212453;
    _collision_strength[OIII][TRANSITION_0_to_1][2] = 2.17448230368;
    _collision_strength[OIII][TRANSITION_0_to_1][3] = 0.115352927697;
    _collision_strength[OIII][TRANSITION_0_to_1][4] = -0.000224722336436;
    _collision_strength[OIII][TRANSITION_0_to_1][5] = -0.115352902446;
    _collision_strength[OIII][TRANSITION_0_to_1][6] = 1.00000001641;
    _collision_strength[OIII][TRANSITION_0_to_2][0] = 0.711355746813;
    _collision_strength[OIII][TRANSITION_0_to_2][1] = 0.00146371349738;
    _collision_strength[OIII][TRANSITION_0_to_2][2] = 1.37408055389;
    _collision_strength[OIII][TRANSITION_0_to_2][3] = 0.0723506999835;
    _collision_strength[OIII][TRANSITION_0_to_2][4] = -0.000136510086638;
    _collision_strength[OIII][TRANSITION_0_to_2][5] = -0.0723506857993;
    _collision_strength[OIII][TRANSITION_0_to_2][6] = 1.00000001479;
    _collision_strength[OIII][TRANSITION_0_to_3][0] = 0.212048796845;
    _collision_strength[OIII][TRANSITION_0_to_3][1] = -0.0308141863599;
    _collision_strength[OIII][TRANSITION_0_to_3][2] = 50.1955663469;
    _collision_strength[OIII][TRANSITION_0_to_3][3] = 0.115631697056;
    _collision_strength[OIII][TRANSITION_0_to_3][4] = 0.00745960829944;
    _collision_strength[OIII][TRANSITION_0_to_3][5] = -0.115633905463;
    _collision_strength[OIII][TRANSITION_0_to_3][6] = 0.999998565297;
    _collision_strength[OIII][TRANSITION_0_to_4][0] = 0.473875264246;
    _collision_strength[OIII][TRANSITION_0_to_4][1] = 0.000161690904972;
    _collision_strength[OIII][TRANSITION_0_to_4][2] = 1.04668052111;
    _collision_strength[OIII][TRANSITION_0_to_4][3] = 0.596870866366;
    _collision_strength[OIII][TRANSITION_0_to_4][4] = 2.35439474228e-05;
    _collision_strength[OIII][TRANSITION_0_to_4][5] = -0.596870889829;
    _collision_strength[OIII][TRANSITION_0_to_4][6] = 0.999999996984;
    _collision_strength[OIII][TRANSITION_1_to_2][0] = 0.728575980601;
    _collision_strength[OIII][TRANSITION_1_to_2][1] = 0.00630786658198;
    _collision_strength[OIII][TRANSITION_1_to_2][2] = 5.83509344037;
    _collision_strength[OIII][TRANSITION_1_to_2][3] = 0.115352923656;
    _collision_strength[OIII][TRANSITION_1_to_2][4] = -0.000599131478658;
    _collision_strength[OIII][TRANSITION_1_to_2][5] = -0.115352859303;
    _collision_strength[OIII][TRANSITION_1_to_2][6] = 1.00000004191;
    _collision_strength[OIII][TRANSITION_1_to_3][0] = 0.212049723828;
    _collision_strength[OIII][TRANSITION_1_to_3][1] = -0.0924411081199;
    _collision_strength[OIII][TRANSITION_1_to_3][2] = 150.585634379;
    _collision_strength[OIII][TRANSITION_1_to_3][3] = 0.271112958338;
    _collision_strength[OIII][TRANSITION_1_to_3][4] = 0.0223785612179;
    _collision_strength[OIII][TRANSITION_1_to_3][5] = -0.271119583518;
    _collision_strength[OIII][TRANSITION_1_to_3][6] = 0.99999816428;
    _collision_strength[OIII][TRANSITION_1_to_4][0] = 0.473874398853;
    _collision_strength[OIII][TRANSITION_1_to_4][1] = 0.000485071869136;
    _collision_strength[OIII][TRANSITION_1_to_4][2] = 3.14006004469;
    _collision_strength[OIII][TRANSITION_1_to_4][3] = 0.207986730789;
    _collision_strength[OIII][TRANSITION_1_to_4][4] = 7.06330121882e-05;
    _collision_strength[OIII][TRANSITION_1_to_4][5] = -0.207986801177;
    _collision_strength[OIII][TRANSITION_1_to_4][6] = 0.999999974034;
    _collision_strength[OIII][TRANSITION_2_to_3][0] = 0.212049682479;
    _collision_strength[OIII][TRANSITION_2_to_3][1] = -0.154067564151;
    _collision_strength[OIII][TRANSITION_2_to_3][2] = 250.975916712;
    _collision_strength[OIII][TRANSITION_2_to_3][3] = -0.346947242556;
    _collision_strength[OIII][TRANSITION_2_to_3][4] = 0.0372974838928;
    _collision_strength[OIII][TRANSITION_2_to_3][5] = 0.346936200866;
    _collision_strength[OIII][TRANSITION_2_to_3][6] = 1.0000023908;
    _collision_strength[OIII][TRANSITION_2_to_4][0] = 0.473773531719;
    _collision_strength[OIII][TRANSITION_2_to_4][1] = 0.000808267411352;
    _collision_strength[OIII][TRANSITION_2_to_4][2] = 5.23703593485;
    _collision_strength[OIII][TRANSITION_2_to_4][3] = 0.346058699298;
    _collision_strength[OIII][TRANSITION_2_to_4][4] = 0.000117951167668;
    _collision_strength[OIII][TRANSITION_2_to_4][5] = -0.346058816732;
    _collision_strength[OIII][TRANSITION_2_to_4][6] = 0.999999973963;
    _collision_strength[OIII][TRANSITION_3_to_4][0] = 1.06251713793;
    _collision_strength[OIII][TRANSITION_3_to_4][1] = 0.000170026677194;
    _collision_strength[OIII][TRANSITION_3_to_4][2] = 0.176555659107;
    _collision_strength[OIII][TRANSITION_3_to_4][3] = 0.115626390808;
    _collision_strength[OIII][TRANSITION_3_to_4][4] = -1.7736910794e-05;
    _collision_strength[OIII][TRANSITION_3_to_4][5] = -0.115626388216;
    _collision_strength[OIII][TRANSITION_3_to_4][6] = 1.00000000168;
  }

  /// NeIII
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 6 and 7
    // ground state: 3P2
    // excited states: 3P1, 3P0, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {643., 921., 25841., 55751.};
    _inverse_statistical_weight[NeIII][0] = 0.2;
    _inverse_statistical_weight[NeIII][1] = 1. / 3.;
    _inverse_statistical_weight[NeIII][2] = 1.;
    _inverse_statistical_weight[NeIII][3] = 0.2;
    _inverse_statistical_weight[NeIII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[NeIII][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[NeIII][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[NeIII][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[NeIII][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[NeIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[NeIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[NeIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[NeIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[NeIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[NeIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _transition_probability[NeIII][TRANSITION_0_to_1] = 5.974e-3;
    _transition_probability[NeIII][TRANSITION_0_to_2] = 2.081e-8;
    _transition_probability[NeIII][TRANSITION_0_to_3] = 0.173;
    _transition_probability[NeIII][TRANSITION_0_to_4] = 3.985e-3;
    _transition_probability[NeIII][TRANSITION_1_to_2] = 1.159e-3;
    _transition_probability[NeIII][TRANSITION_1_to_3] = 5.344e-2;
    _transition_probability[NeIII][TRANSITION_1_to_4] = 2.028;
    _transition_probability[NeIII][TRANSITION_2_to_3] = 8.269e-6;
    _transition_probability[NeIII][TRANSITION_2_to_4] = 0.;
    _transition_probability[NeIII][TRANSITION_3_to_4] = 2.563;
    // our own fits to the data of Butler & Zeippen (1994)
    // these fits were made with the script data/linecooling/gamma_NeIII.py
    // (this script also outputs the code below)
    _collision_strength[NeIII][TRANSITION_0_to_1][0] = 0.229212961855;
    _collision_strength[NeIII][TRANSITION_0_to_1][1] = 0.391241413653;
    _collision_strength[NeIII][TRANSITION_0_to_1][2] = -98.8699646747;
    _collision_strength[NeIII][TRANSITION_0_to_1][3] = 3.33490329495e-07;
    _collision_strength[NeIII][TRANSITION_0_to_1][4] = -0.0315920221406;
    _collision_strength[NeIII][TRANSITION_0_to_1][5] = -1.79860927224e-34;
    _collision_strength[NeIII][TRANSITION_0_to_1][6] = 6.24902182391;
    _collision_strength[NeIII][TRANSITION_0_to_2][0] = 0.242231485943;
    _collision_strength[NeIII][TRANSITION_0_to_2][1] = 0.118413569516;
    _collision_strength[NeIII][TRANSITION_0_to_2][2] = -28.631450153;
    _collision_strength[NeIII][TRANSITION_0_to_2][3] = 0.00020164783497;
    _collision_strength[NeIII][TRANSITION_0_to_2][4] = -0.0104670084568;
    _collision_strength[NeIII][TRANSITION_0_to_2][5] = -0.00020065502745;
    _collision_strength[NeIII][TRANSITION_0_to_2][6] = 1.0003625347;
    _collision_strength[NeIII][TRANSITION_0_to_3][0] = 0.994550005514;
    _collision_strength[NeIII][TRANSITION_0_to_3][1] = -9.74857730722e-06;
    _collision_strength[NeIII][TRANSITION_0_to_3][2] = 0.813407066624;
    _collision_strength[NeIII][TRANSITION_0_to_3][3] = 1.47512404499e-07;
    _collision_strength[NeIII][TRANSITION_0_to_3][4] = 7.93314790996e-07;
    _collision_strength[NeIII][TRANSITION_0_to_3][5] = -1.47337527572e-07;
    _collision_strength[NeIII][TRANSITION_0_to_3][6] = 1.0000982771;
    _collision_strength[NeIII][TRANSITION_0_to_4][0] = 1.10193674216;
    _collision_strength[NeIII][TRANSITION_0_to_4][1] = -8.06762952432e-06;
    _collision_strength[NeIII][TRANSITION_0_to_4][2] = 0.0418940659381;
    _collision_strength[NeIII][TRANSITION_0_to_4][3] = 2.09442315818e-08;
    _collision_strength[NeIII][TRANSITION_0_to_4][4] = 8.07352815912e-07;
    _collision_strength[NeIII][TRANSITION_0_to_4][5] = -2.10320915394e-08;
    _collision_strength[NeIII][TRANSITION_0_to_4][6] = 0.999691826017;
    _collision_strength[NeIII][TRANSITION_1_to_2][0] = 0.382687205434;
    _collision_strength[NeIII][TRANSITION_1_to_2][1] = 0.0454576087822;
    _collision_strength[NeIII][TRANSITION_1_to_2][2] = -8.08623902229;
    _collision_strength[NeIII][TRANSITION_1_to_2][3] = -7.04290896513e-05;
    _collision_strength[NeIII][TRANSITION_1_to_2][4] = -0.00421458041452;
    _collision_strength[NeIII][TRANSITION_1_to_2][5] = 7.08575087635e-05;
    _collision_strength[NeIII][TRANSITION_1_to_2][6] = 0.999550531942;
    _collision_strength[NeIII][TRANSITION_1_to_3][0] = 0.998068711543;
    _collision_strength[NeIII][TRANSITION_1_to_3][1] = -8.15415554133e-06;
    _collision_strength[NeIII][TRANSITION_1_to_3][2] = 0.475826627698;
    _collision_strength[NeIII][TRANSITION_1_to_3][3] = -8.84425323564e-08;
    _collision_strength[NeIII][TRANSITION_1_to_3][4] = 7.00594666413e-07;
    _collision_strength[NeIII][TRANSITION_1_to_3][5] = 8.85242619285e-08;
    _collision_strength[NeIII][TRANSITION_1_to_3][6] = 0.999921004822;
    _collision_strength[NeIII][TRANSITION_1_to_4][0] = 0.355318972635;
    _collision_strength[NeIII][TRANSITION_1_to_4][1] = 0.0169307413432;
    _collision_strength[NeIII][TRANSITION_1_to_4][2] = -0.333065014807;
    _collision_strength[NeIII][TRANSITION_1_to_4][3] = -1.00690213878e-08;
    _collision_strength[NeIII][TRANSITION_1_to_4][4] = -0.00296274662344;
    _collision_strength[NeIII][TRANSITION_1_to_4][5] = 0.00214275937489;
    _collision_strength[NeIII][TRANSITION_1_to_4][6] = 0.190431927829;
    _collision_strength[NeIII][TRANSITION_2_to_3][0] = 0.978821637591;
    _collision_strength[NeIII][TRANSITION_2_to_3][1] = 7.06633577309e-06;
    _collision_strength[NeIII][TRANSITION_2_to_3][2] = 0.178784485715;
    _collision_strength[NeIII][TRANSITION_2_to_3][3] = 3.77212505483e-08;
    _collision_strength[NeIII][TRANSITION_2_to_3][4] = -7.76303145206e-07;
    _collision_strength[NeIII][TRANSITION_2_to_3][5] = -3.75278757078e-08;
    _collision_strength[NeIII][TRANSITION_2_to_3][6] = 1.00039864915;
    _collision_strength[NeIII][TRANSITION_2_to_4][0] = 0.369732170405;
    _collision_strength[NeIII][TRANSITION_2_to_4][1] = 0.00310820609142;
    _collision_strength[NeIII][TRANSITION_2_to_4][2] = 0.209181124787;
    _collision_strength[NeIII][TRANSITION_2_to_4][3] = 1.05339202675e-05;
    _collision_strength[NeIII][TRANSITION_2_to_4][4] = -0.000293668434264;
    _collision_strength[NeIII][TRANSITION_2_to_4][5] = -1.04869546902e-05;
    _collision_strength[NeIII][TRANSITION_2_to_4][6] = 1.0003424372;
    _collision_strength[NeIII][TRANSITION_3_to_4][0] = 1.05810601642;
    _collision_strength[NeIII][TRANSITION_3_to_4][1] = -4.02053295696e-05;
    _collision_strength[NeIII][TRANSITION_3_to_4][2] = 0.187646793296;
    _collision_strength[NeIII][TRANSITION_3_to_4][3] = -1.18487150347e-07;
    _collision_strength[NeIII][TRANSITION_3_to_4][4] = 4.26737262075e-06;
    _collision_strength[NeIII][TRANSITION_3_to_4][5] = 1.17809039437e-07;
    _collision_strength[NeIII][TRANSITION_3_to_4][6] = 1.0004318279;
  }

  /// SII
  {
    // data from Tayal & Zatsarinny (2010), tables 1 and 3
    // ground state: 4S3/2
    // excited states: 2D3/2, 2D5/2, 2P1/2, 2P3/2
    // in eV
    const double energy_levels[4] = {1.842, 1.845, 3.041, 3.046};
    _inverse_statistical_weight[SII][0] = 0.25;
    _inverse_statistical_weight[SII][1] = 0.25;
    _inverse_statistical_weight[SII][2] = 1. / 6.;
    _inverse_statistical_weight[SII][3] = 0.5;
    _inverse_statistical_weight[SII][4] = 0.25;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[SII][TRANSITION_0_to_1] = energy_levels[0] * eV_over_k;
    _energy_difference[SII][TRANSITION_0_to_2] = energy_levels[1] * eV_over_k;
    _energy_difference[SII][TRANSITION_0_to_3] = energy_levels[2] * eV_over_k;
    _energy_difference[SII][TRANSITION_0_to_4] = energy_levels[3] * eV_over_k;
    _energy_difference[SII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * eV_over_k;
    _energy_difference[SII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * eV_over_k;
    _energy_difference[SII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * eV_over_k;
    _energy_difference[SII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * eV_over_k;
    _energy_difference[SII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * eV_over_k;
    _energy_difference[SII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * eV_over_k;
    // we take the sum of all contributions (in s^-1)
    _transition_probability[SII][TRANSITION_0_to_1] = 6.32e-4;
    _transition_probability[SII][TRANSITION_0_to_2] = 2.20e-4;
    _transition_probability[SII][TRANSITION_0_to_3] = 7.64e-2;
    _transition_probability[SII][TRANSITION_0_to_4] = 1.90e-1;
    _transition_probability[SII][TRANSITION_1_to_2] = 1.71e-7;
    _transition_probability[SII][TRANSITION_1_to_3] = 1.47e-1;
    _transition_probability[SII][TRANSITION_1_to_4] = 1.165e-1;
    _transition_probability[SII][TRANSITION_2_to_3] = 7.16e-2;
    _transition_probability[SII][TRANSITION_2_to_4] = 1.61e-1;
    _transition_probability[SII][TRANSITION_3_to_4] = 2.43e-7;
    // our own fits to the data of Tayal & Zatsarinny (2010)
    // these fits were made with the script data/linecooling/gamma_SII.py
    // (this script also outputs the code below)
    _collision_strength[SII][TRANSITION_0_to_1][0] = 1.03641663445;
    _collision_strength[SII][TRANSITION_0_to_1][1] = -0.000249794248385;
    _collision_strength[SII][TRANSITION_0_to_1][2] = 2.20349579749;
    _collision_strength[SII][TRANSITION_0_to_1][3] = -0.120817610429;
    _collision_strength[SII][TRANSITION_0_to_1][4] = 2.39631107845e-05;
    _collision_strength[SII][TRANSITION_0_to_1][5] = 0.120817607822;
    _collision_strength[SII][TRANSITION_0_to_1][6] = 1.00000000161;
    _collision_strength[SII][TRANSITION_0_to_2][0] = 1.03717217955;
    _collision_strength[SII][TRANSITION_0_to_2][1] = -0.00039155644199;
    _collision_strength[SII][TRANSITION_0_to_2][2] = 3.28566942345;
    _collision_strength[SII][TRANSITION_0_to_2][3] = -0.121147740734;
    _collision_strength[SII][TRANSITION_0_to_2][4] = 3.77925782749e-05;
    _collision_strength[SII][TRANSITION_0_to_2][5] = 0.121147736484;
    _collision_strength[SII][TRANSITION_0_to_2][6] = 1.00000000262;
    _collision_strength[SII][TRANSITION_0_to_3][0] = 1.00775072872;
    _collision_strength[SII][TRANSITION_0_to_3][1] = -2.72593102367e-05;
    _collision_strength[SII][TRANSITION_0_to_3][2] = 0.638728107106;
    _collision_strength[SII][TRANSITION_0_to_3][3] = -0.123901293284;
    _collision_strength[SII][TRANSITION_0_to_3][4] = 3.54901615379e-06;
    _collision_strength[SII][TRANSITION_0_to_3][5] = 0.123901291979;
    _collision_strength[SII][TRANSITION_0_to_3][6] = 1.00000000081;
    _collision_strength[SII][TRANSITION_0_to_4][0] = 1.01789525407;
    _collision_strength[SII][TRANSITION_0_to_4][1] = -7.35398931928e-05;
    _collision_strength[SII][TRANSITION_0_to_4][2] = 1.19790395953;
    _collision_strength[SII][TRANSITION_0_to_4][3] = -0.121254946562;
    _collision_strength[SII][TRANSITION_0_to_4][4] = 8.85704384397e-06;
    _collision_strength[SII][TRANSITION_0_to_4][5] = 0.12125494385;
    _collision_strength[SII][TRANSITION_0_to_4][6] = 1.00000000172;
    _collision_strength[SII][TRANSITION_1_to_2][0] = 1.05197469278;
    _collision_strength[SII][TRANSITION_1_to_2][1] = -0.00110099988803;
    _collision_strength[SII][TRANSITION_1_to_2][2] = 5.72662457909;
    _collision_strength[SII][TRANSITION_1_to_2][3] = -0.12184170716;
    _collision_strength[SII][TRANSITION_1_to_2][4] = 0.000107776630232;
    _collision_strength[SII][TRANSITION_1_to_2][5] = 0.121841695288;
    _collision_strength[SII][TRANSITION_1_to_2][6] = 1.00000000723;
    _collision_strength[SII][TRANSITION_1_to_3][0] = 0.964511857799;
    _collision_strength[SII][TRANSITION_1_to_3][1] = -1.30353220838e-05;
    _collision_strength[SII][TRANSITION_1_to_3][2] = 1.96892680842;
    _collision_strength[SII][TRANSITION_1_to_3][3] = -0.535678236601;
    _collision_strength[SII][TRANSITION_1_to_3][4] = 2.60395988847e-06;
    _collision_strength[SII][TRANSITION_1_to_3][5] = 0.535678235452;
    _collision_strength[SII][TRANSITION_1_to_3][6] = 1.00000000016;
    _collision_strength[SII][TRANSITION_1_to_4][0] = 1.05748558757;
    _collision_strength[SII][TRANSITION_1_to_4][1] = -0.000213382981323;
    _collision_strength[SII][TRANSITION_1_to_4][2] = 1.64461231866;
    _collision_strength[SII][TRANSITION_1_to_4][3] = -0.199564473316;
    _collision_strength[SII][TRANSITION_1_to_4][4] = 2.16292515466e-05;
    _collision_strength[SII][TRANSITION_1_to_4][5] = 0.199564470263;
    _collision_strength[SII][TRANSITION_1_to_4][6] = 1.00000000115;
    _collision_strength[SII][TRANSITION_2_to_3][0] = 1.02405773355;
    _collision_strength[SII][TRANSITION_2_to_3][1] = -0.000132943341914;
    _collision_strength[SII][TRANSITION_2_to_3][2] = 1.55838205223;
    _collision_strength[SII][TRANSITION_2_to_3][3] = -0.282405174182;
    _collision_strength[SII][TRANSITION_2_to_3][4] = 1.3709801875e-05;
    _collision_strength[SII][TRANSITION_2_to_3][5] = 0.282405171934;
    _collision_strength[SII][TRANSITION_2_to_3][6] = 1.0000000006;
    _collision_strength[SII][TRANSITION_2_to_4][0] = 1.05203411618;
    _collision_strength[SII][TRANSITION_2_to_4][1] = -0.000370323954089;
    _collision_strength[SII][TRANSITION_2_to_4][2] = 2.9052960213;
    _collision_strength[SII][TRANSITION_2_to_4][3] = -0.120622459895;
    _collision_strength[SII][TRANSITION_2_to_4][4] = 3.77620972842e-05;
    _collision_strength[SII][TRANSITION_2_to_4][5] = 0.120622454552;
    _collision_strength[SII][TRANSITION_2_to_4][6] = 1.00000000333;
    _collision_strength[SII][TRANSITION_3_to_4][0] = 1.04324836606;
    _collision_strength[SII][TRANSITION_3_to_4][1] = -0.000111829064916;
    _collision_strength[SII][TRANSITION_3_to_4][2] = 1.28208938549;
    _collision_strength[SII][TRANSITION_3_to_4][3] = -0.0254376510096;
    _collision_strength[SII][TRANSITION_3_to_4][4] = 1.21897286763e-05;
    _collision_strength[SII][TRANSITION_3_to_4][5] = 0.0254376483263;
    _collision_strength[SII][TRANSITION_3_to_4][6] = 1.00000000807;
  }

  /// SIII
  {
    // data from Mendoza & Zeippen (1982), table 2 and 3
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    // in Ry
    const double energy_levels[4] = {0.002708, 0.007586, 0.103157, 0.247532};
    _inverse_statistical_weight[SIII][0] = 1.;
    _inverse_statistical_weight[SIII][1] = 1. / 3.;
    _inverse_statistical_weight[SIII][2] = 0.2;
    _inverse_statistical_weight[SIII][3] = 0.2;
    _inverse_statistical_weight[SIII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[SIII][TRANSITION_0_to_1] = energy_levels[0] * Ry_over_k;
    _energy_difference[SIII][TRANSITION_0_to_2] = energy_levels[1] * Ry_over_k;
    _energy_difference[SIII][TRANSITION_0_to_3] = energy_levels[2] * Ry_over_k;
    _energy_difference[SIII][TRANSITION_0_to_4] = energy_levels[3] * Ry_over_k;
    _energy_difference[SIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * Ry_over_k;
    _energy_difference[SIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * Ry_over_k;
    _energy_difference[SIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * Ry_over_k;
    _energy_difference[SIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * Ry_over_k;
    _energy_difference[SIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * Ry_over_k;
    _energy_difference[SIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * Ry_over_k;
    // we take the sum of all contributions (in s^-1)
    _transition_probability[SIII][TRANSITION_0_to_1] = 4.72e-4;
    _transition_probability[SIII][TRANSITION_0_to_2] = 4.61e-8;
    _transition_probability[SIII][TRANSITION_0_to_3] = 5.82e-6;
    _transition_probability[SIII][TRANSITION_0_to_4] = 0.;
    _transition_probability[SIII][TRANSITION_1_to_2] = 2.07e-3;
    _transition_probability[SIII][TRANSITION_1_to_3] = 2.20e-2;
    _transition_probability[SIII][TRANSITION_1_to_4] = 7.96e-1;
    _transition_probability[SIII][TRANSITION_2_to_3] = 5.76e-2;
    _transition_probability[SIII][TRANSITION_2_to_4] = 1.05e-2;
    _transition_probability[SIII][TRANSITION_3_to_4] = 2.22;
    // our own fits to the data of Hudson, Ramsbottom & Scott (2012)
    // these fits were made with the script data/linecooling/gamma_SIII.py
    // (this script also outputs the code below)
    _collision_strength[SIII][TRANSITION_0_to_1][0] = 0.319834036674;
    _collision_strength[SIII][TRANSITION_0_to_1][1] = 0.473094494664;
    _collision_strength[SIII][TRANSITION_0_to_1][2] = 23.4611409315;
    _collision_strength[SIII][TRANSITION_0_to_1][3] = 1.81292268139e-07;
    _collision_strength[SIII][TRANSITION_0_to_1][4] = -0.0391016530041;
    _collision_strength[SIII][TRANSITION_0_to_1][5] = 1.02381727611e-08;
    _collision_strength[SIII][TRANSITION_0_to_1][6] = -5.38956433778;
    _collision_strength[SIII][TRANSITION_0_to_2][0] = 0.946895106493;
    _collision_strength[SIII][TRANSITION_0_to_2][1] = -0.000253931252088;
    _collision_strength[SIII][TRANSITION_0_to_2][2] = 1.42893199246;
    _collision_strength[SIII][TRANSITION_0_to_2][3] = -0.0143421985165;
    _collision_strength[SIII][TRANSITION_0_to_2][4] = 3.36498345396e-05;
    _collision_strength[SIII][TRANSITION_0_to_2][5] = 0.0143421876514;
    _collision_strength[SIII][TRANSITION_0_to_2][6] = 1.00000005821;
    _collision_strength[SIII][TRANSITION_0_to_3][0] = 1.09296096143;
    _collision_strength[SIII][TRANSITION_0_to_3][1] = -0.000130949694096;
    _collision_strength[SIII][TRANSITION_0_to_3][2] = 0.40981640082;
    _collision_strength[SIII][TRANSITION_0_to_3][3] = -0.00193270090945;
    _collision_strength[SIII][TRANSITION_0_to_3][4] = 1.40040822412e-05;
    _collision_strength[SIII][TRANSITION_0_to_3][5] = 0.00193269826308;
    _collision_strength[SIII][TRANSITION_0_to_3][6] = 1.0000001045;
    _collision_strength[SIII][TRANSITION_0_to_4][0] = 0.268430758822;
    _collision_strength[SIII][TRANSITION_0_to_4][1] = -0.161012416751;
    _collision_strength[SIII][TRANSITION_0_to_4][2] = -1.727282253;
    _collision_strength[SIII][TRANSITION_0_to_4][3] = -7.73439548718e-08;
    _collision_strength[SIII][TRANSITION_0_to_4][4] = -0.0228575861268;
    _collision_strength[SIII][TRANSITION_0_to_4][5] = 0.220495114378;
    _collision_strength[SIII][TRANSITION_0_to_4][6] = 0.0599587031067;
    _collision_strength[SIII][TRANSITION_1_to_2][0] = 0.945260316234;
    _collision_strength[SIII][TRANSITION_1_to_2][1] = 0.00126278546768;
    _collision_strength[SIII][TRANSITION_1_to_2][2] = 6.58127073323;
    _collision_strength[SIII][TRANSITION_1_to_2][3] = 0.00656821073797;
    _collision_strength[SIII][TRANSITION_1_to_2][4] = -0.000120335068449;
    _collision_strength[SIII][TRANSITION_1_to_2][5] = -0.00656820056244;
    _collision_strength[SIII][TRANSITION_1_to_2][6] = 1.00000011531;
    _collision_strength[SIII][TRANSITION_1_to_3][0] = 1.08018916158;
    _collision_strength[SIII][TRANSITION_1_to_3][1] = -0.000410120544537;
    _collision_strength[SIII][TRANSITION_1_to_3][2] = 1.33349052392;
    _collision_strength[SIII][TRANSITION_1_to_3][3] = 0.00139791118784;
    _collision_strength[SIII][TRANSITION_1_to_3][4] = 4.41269104115e-05;
    _collision_strength[SIII][TRANSITION_1_to_3][5] = -0.00139791975694;
    _collision_strength[SIII][TRANSITION_1_to_3][6] = 0.999999531792;
    _collision_strength[SIII][TRANSITION_1_to_4][0] = -0.00314229593184;
    _collision_strength[SIII][TRANSITION_1_to_4][1] = 0.0273739699311;
    _collision_strength[SIII][TRANSITION_1_to_4][2] = -9.08588844556;
    _collision_strength[SIII][TRANSITION_1_to_4][3] = -0.00532056355739;
    _collision_strength[SIII][TRANSITION_1_to_4][4] = 0.0278125010909;
    _collision_strength[SIII][TRANSITION_1_to_4][5] = 0.00534754896895;
    _collision_strength[SIII][TRANSITION_1_to_4][6] = 0.999565574858;
    _collision_strength[SIII][TRANSITION_2_to_3][0] = 1.08428172406;
    _collision_strength[SIII][TRANSITION_2_to_3][1] = -0.000727729853996;
    _collision_strength[SIII][TRANSITION_2_to_3][2] = 2.41274280077;
    _collision_strength[SIII][TRANSITION_2_to_3][3] = -0.00181282426352;
    _collision_strength[SIII][TRANSITION_2_to_3][4] = 7.77416341249e-05;
    _collision_strength[SIII][TRANSITION_2_to_3][5] = 0.00181280958775;
    _collision_strength[SIII][TRANSITION_2_to_3][6] = 1.0000006179;
    _collision_strength[SIII][TRANSITION_2_to_4][0] = -0.454650835236;
    _collision_strength[SIII][TRANSITION_2_to_4][1] = -32.2631830151;
    _collision_strength[SIII][TRANSITION_2_to_4][2] = 2378.98316062;
    _collision_strength[SIII][TRANSITION_2_to_4][3] = 0.00245376919637;
    _collision_strength[SIII][TRANSITION_2_to_4][4] = 5.1731551752;
    _collision_strength[SIII][TRANSITION_2_to_4][5] = -2.30897023902e-06;
    _collision_strength[SIII][TRANSITION_2_to_4][6] = 1.56838057651;
    _collision_strength[SIII][TRANSITION_3_to_4][0] = 0.383531726856;
    _collision_strength[SIII][TRANSITION_3_to_4][1] = -0.0883803729726;
    _collision_strength[SIII][TRANSITION_3_to_4][2] = 49.4449183698;
    _collision_strength[SIII][TRANSITION_3_to_4][3] = -0.0110688636912;
    _collision_strength[SIII][TRANSITION_3_to_4][4] = 0.0149332210578;
    _collision_strength[SIII][TRANSITION_3_to_4][5] = 0.0110640993744;
    _collision_strength[SIII][TRANSITION_3_to_4][6] = 1.00003260333;
  }

  /// CII
  {
    // data from Froese Fischer & Tachiev (2004), table 2
    // ground state: 2P1/2
    // excited states: 2P3/2, 4P1/2, 4P3/2, 4P5/2
    // in cm^-1
    const double energy_levels[4] = {63.67, 43057.99, 43080.21, 43108.66};
    _inverse_statistical_weight[CII][0] = 0.5;
    _inverse_statistical_weight[CII][1] = 0.25;
    _inverse_statistical_weight[CII][2] = 0.5;
    _inverse_statistical_weight[CII][3] = 0.25;
    _inverse_statistical_weight[CII][4] = 1. / 6.;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[CII][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[CII][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[CII][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[CII][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[CII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[CII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[CII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[CII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[CII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[CII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we sum contributions of all types (in s^-1)
    _transition_probability[CII][TRANSITION_0_to_1] = 2.321e-6;
    _transition_probability[CII][TRANSITION_0_to_2] = 6.136e1;
    _transition_probability[CII][TRANSITION_0_to_3] = 1.463;
    _transition_probability[CII][TRANSITION_0_to_4] = 8.177e-4;
    _transition_probability[CII][TRANSITION_1_to_2] = 6.929e1;
    _transition_probability[CII][TRANSITION_1_to_3] = 8.853;
    _transition_probability[CII][TRANSITION_1_to_4] = 4.477e1;
    _transition_probability[CII][TRANSITION_2_to_3] = 2.467e-7;
    _transition_probability[CII][TRANSITION_2_to_4] = 3.571e-14;
    _transition_probability[CII][TRANSITION_3_to_4] = 3.725e-7;
    // our own fits to the data of Tayal (2008)
    // these fits were made with the script data/linecooling/gamma_CII.py
    // (this script also outputs the code below)
    _collision_strength[CII][TRANSITION_0_to_1][0] = 0.0480505832358;
    _collision_strength[CII][TRANSITION_0_to_1][1] = -5.10087630511;
    _collision_strength[CII][TRANSITION_0_to_1][2] = 1436.89664127;
    _collision_strength[CII][TRANSITION_0_to_1][3] = 1561.94244352;
    _collision_strength[CII][TRANSITION_0_to_1][4] = 0.742924478556;
    _collision_strength[CII][TRANSITION_0_to_1][5] = -1561.94264165;
    _collision_strength[CII][TRANSITION_0_to_1][6] = 0.999999990311;
    _collision_strength[CII][TRANSITION_0_to_2][0] = 1.13322434366;
    _collision_strength[CII][TRANSITION_0_to_2][1] = -2.73646504879e-05;
    _collision_strength[CII][TRANSITION_0_to_2][2] = 0.104888308389;
    _collision_strength[CII][TRANSITION_0_to_2][3] = 0.49980900604;
    _collision_strength[CII][TRANSITION_0_to_2][4] = 2.81074173787e-06;
    _collision_strength[CII][TRANSITION_0_to_2][5] = -0.499809006492;
    _collision_strength[CII][TRANSITION_0_to_2][6] = 0.999999999931;
    _collision_strength[CII][TRANSITION_0_to_3][0] = 1.03908907832;
    _collision_strength[CII][TRANSITION_0_to_3][1] = -4.40567142559e-05;
    _collision_strength[CII][TRANSITION_0_to_3][2] = 0.309950968041;
    _collision_strength[CII][TRANSITION_0_to_3][3] = 0.411097707865;
    _collision_strength[CII][TRANSITION_0_to_3][4] = 4.69915372252e-06;
    _collision_strength[CII][TRANSITION_0_to_3][5] = -0.411097708841;
    _collision_strength[CII][TRANSITION_0_to_3][6] = 0.999999999817;
    _collision_strength[CII][TRANSITION_0_to_4][0] = 0.940734820256;
    _collision_strength[CII][TRANSITION_0_to_4][1] = 1.37370102592e-05;
    _collision_strength[CII][TRANSITION_0_to_4][2] = 0.362085513794;
    _collision_strength[CII][TRANSITION_0_to_4][3] = 0.366342804713;
    _collision_strength[CII][TRANSITION_0_to_4][4] = -5.16981296168e-07;
    _collision_strength[CII][TRANSITION_0_to_4][5] = -0.366342805653;
    _collision_strength[CII][TRANSITION_0_to_4][6] = 0.999999999794;
    _collision_strength[CII][TRANSITION_1_to_2][0] = 1.21127092292;
    _collision_strength[CII][TRANSITION_1_to_2][1] = -9.94067404214e-06;
    _collision_strength[CII][TRANSITION_1_to_2][2] = 0.0365338010715;
    _collision_strength[CII][TRANSITION_1_to_2][3] = 0.50013648256;
    _collision_strength[CII][TRANSITION_1_to_2][4] = 1.01657928264e-06;
    _collision_strength[CII][TRANSITION_1_to_2][5] = -0.500136482722;
    _collision_strength[CII][TRANSITION_1_to_2][6] = 0.999999999975;
    _collision_strength[CII][TRANSITION_1_to_3][0] = 1.15008027301;
    _collision_strength[CII][TRANSITION_1_to_3][1] = -3.72283249595e-05;
    _collision_strength[CII][TRANSITION_1_to_3][2] = 0.165574827638;
    _collision_strength[CII][TRANSITION_1_to_3][3] = 0.49993982796;
    _collision_strength[CII][TRANSITION_1_to_3][4] = 3.81426383305e-06;
    _collision_strength[CII][TRANSITION_1_to_3][5] = -0.499939828584;
    _collision_strength[CII][TRANSITION_1_to_3][6] = 0.999999999905;
    _collision_strength[CII][TRANSITION_1_to_4][0] = 0.993974665119;
    _collision_strength[CII][TRANSITION_1_to_4][1] = -0.000117455969969;
    _collision_strength[CII][TRANSITION_1_to_4][2] = 1.22975481187;
    _collision_strength[CII][TRANSITION_1_to_4][3] = 0.786974771901;
    _collision_strength[CII][TRANSITION_1_to_4][4] = 1.3263664083e-05;
    _collision_strength[CII][TRANSITION_1_to_4][5] = -0.786974775449;
    _collision_strength[CII][TRANSITION_1_to_4][6] = 0.99999999965;
    _collision_strength[CII][TRANSITION_2_to_3][0] = 0.8353074866;
    _collision_strength[CII][TRANSITION_2_to_3][1] = -0.000459275718981;
    _collision_strength[CII][TRANSITION_2_to_3][2] = 1.96939409603;
    _collision_strength[CII][TRANSITION_2_to_3][3] = 348.77077441;
    _collision_strength[CII][TRANSITION_2_to_3][4] = 8.00777479566e-05;
    _collision_strength[CII][TRANSITION_2_to_3][5] = -348.77077445;
    _collision_strength[CII][TRANSITION_2_to_3][6] = 0.999999999991;
    _collision_strength[CII][TRANSITION_2_to_4][0] = 0.516732024641;
    _collision_strength[CII][TRANSITION_2_to_4][1] = -0.0017944348363;
    _collision_strength[CII][TRANSITION_2_to_4][2] = 15.5044834402;
    _collision_strength[CII][TRANSITION_2_to_4][3] = 199.020768398;
    _collision_strength[CII][TRANSITION_2_to_4][4] = 0.00104646848633;
    _collision_strength[CII][TRANSITION_2_to_4][5] = -199.020769192;
    _collision_strength[CII][TRANSITION_2_to_4][6] = 0.999999999688;
    _collision_strength[CII][TRANSITION_3_to_4][0] = 0.660041604974;
    _collision_strength[CII][TRANSITION_3_to_4][1] = -0.002364929605;
    _collision_strength[CII][TRANSITION_3_to_4][2] = 15.6088769634;
    _collision_strength[CII][TRANSITION_3_to_4][3] = 13186.8831818;
    _collision_strength[CII][TRANSITION_3_to_4][4] = 0.000712740711772;
    _collision_strength[CII][TRANSITION_3_to_4][5] = -13186.8831823;
    _collision_strength[CII][TRANSITION_3_to_4][6] = 0.999999999997;
  }

  /// CIII
  {
    // data from Froese Fischer & Tachiev (2004), table 1
    // ground state: 1S0
    // excited states: 3P0, 3P1, 3P2, 1P1
    // in cm^-1
    const double energy_levels[4] = {52391.87, 52415.53, 52472.16, 102446.94};
    _inverse_statistical_weight[CIII][0] = 1.;
    _inverse_statistical_weight[CIII][1] = 1.;
    _inverse_statistical_weight[CIII][2] = 1. / 3.;
    _inverse_statistical_weight[CIII][3] = 0.2;
    _inverse_statistical_weight[CIII][4] = 1. / 3.;
    // convert energy levels to energy differences (and convert units)
    _energy_difference[CIII][TRANSITION_0_to_1] = energy_levels[0] * hc_over_k;
    _energy_difference[CIII][TRANSITION_0_to_2] = energy_levels[1] * hc_over_k;
    _energy_difference[CIII][TRANSITION_0_to_3] = energy_levels[2] * hc_over_k;
    _energy_difference[CIII][TRANSITION_0_to_4] = energy_levels[3] * hc_over_k;
    _energy_difference[CIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _energy_difference[CIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _energy_difference[CIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _energy_difference[CIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _energy_difference[CIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _energy_difference[CIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we sum contributions of all types (in s^-1)
    _transition_probability[CIII][TRANSITION_0_to_1] = 0.;
    _transition_probability[CIII][TRANSITION_0_to_2] = 1.040e2;
    _transition_probability[CIII][TRANSITION_0_to_3] = 5.216e-3;
    _transition_probability[CIII][TRANSITION_0_to_4] = 1.769e9;
    _transition_probability[CIII][TRANSITION_1_to_2] = 2.381e-7;
    _transition_probability[CIII][TRANSITION_1_to_3] = 1.486e-13;
    _transition_probability[CIII][TRANSITION_1_to_4] = 1.595e-3;
    _transition_probability[CIII][TRANSITION_2_to_3] = 2.450e-6;
    _transition_probability[CIII][TRANSITION_2_to_4] = 1.269e-3;
    _transition_probability[CIII][TRANSITION_3_to_4] = 2.036e-3;
    // our own fits to the data of Berrington et al. (1985)
    // these fits were made with the script data/linecooling/gamma_CIII.py
    // (this script also outputs the code below)
    _collision_strength[CIII][TRANSITION_0_to_1][0] = 0.680381936743;
    _collision_strength[CIII][TRANSITION_0_to_1][1] = 0.000595111454514;
    _collision_strength[CIII][TRANSITION_0_to_1][2] = 0.999774036843;
    _collision_strength[CIII][TRANSITION_0_to_1][3] = -0.868336006673;
    _collision_strength[CIII][TRANSITION_0_to_1][4] = -5.24291039737e-05;
    _collision_strength[CIII][TRANSITION_0_to_1][5] = 0.868336008271;
    _collision_strength[CIII][TRANSITION_0_to_1][6] = 0.999999999877;
    _collision_strength[CIII][TRANSITION_0_to_2][0] = 0.993110668289;
    _collision_strength[CIII][TRANSITION_0_to_2][1] = 2.00775854613e-05;
    _collision_strength[CIII][TRANSITION_0_to_2][2] = 0.341938652078;
    _collision_strength[CIII][TRANSITION_0_to_2][3] = -0.077525729879;
    _collision_strength[CIII][TRANSITION_0_to_2][4] = -1.94419933081e-06;
    _collision_strength[CIII][TRANSITION_0_to_2][5] = 0.0775257299714;
    _collision_strength[CIII][TRANSITION_0_to_2][6] = 0.999999999919;
    _collision_strength[CIII][TRANSITION_0_to_3][0] = 0.989348686947;
    _collision_strength[CIII][TRANSITION_0_to_3][1] = 3.62365756685e-05;
    _collision_strength[CIII][TRANSITION_0_to_3][2] = 0.586929942311;
    _collision_strength[CIII][TRANSITION_0_to_3][3] = -0.0245235186989;
    _collision_strength[CIII][TRANSITION_0_to_3][4] = -3.49480183379e-06;
    _collision_strength[CIII][TRANSITION_0_to_3][5] = 0.0245235188626;
    _collision_strength[CIII][TRANSITION_0_to_3][6] = 0.999999999548;
    _collision_strength[CIII][TRANSITION_0_to_4][0] = 1.06360193686;
    _collision_strength[CIII][TRANSITION_0_to_4][1] = -3.04483388297e-05;
    _collision_strength[CIII][TRANSITION_0_to_4][2] = 2.23929533724;
    _collision_strength[CIII][TRANSITION_0_to_4][3] = -0.00938462756431;
    _collision_strength[CIII][TRANSITION_0_to_4][4] = 2.88405160957e-06;
    _collision_strength[CIII][TRANSITION_0_to_4][5] = 0.00938462748276;
    _collision_strength[CIII][TRANSITION_0_to_4][6] = 1.00000000057;
    _collision_strength[CIII][TRANSITION_1_to_2][0] = 0.738971130118;
    _collision_strength[CIII][TRANSITION_1_to_2][1] = 0.00386367008165;
    _collision_strength[CIII][TRANSITION_1_to_2][2] = 2.49111221111;
    _collision_strength[CIII][TRANSITION_1_to_2][3] = 0.00315734011239;
    _collision_strength[CIII][TRANSITION_1_to_2][4] = -0.000340295038591;
    _collision_strength[CIII][TRANSITION_1_to_2][5] = -0.00315733000806;
    _collision_strength[CIII][TRANSITION_1_to_2][6] = 1.00000021346;
    _collision_strength[CIII][TRANSITION_1_to_3][0] = 0.867329367873;
    _collision_strength[CIII][TRANSITION_1_to_3][1] = 0.00134115573239;
    _collision_strength[CIII][TRANSITION_1_to_3][2] = -0.397414520326;
    _collision_strength[CIII][TRANSITION_1_to_3][3] = -0.114393956811;
    _collision_strength[CIII][TRANSITION_1_to_3][4] = -0.000119702687795;
    _collision_strength[CIII][TRANSITION_1_to_3][5] = 0.11439396073;
    _collision_strength[CIII][TRANSITION_1_to_3][6] = 0.999999997705;
    _collision_strength[CIII][TRANSITION_1_to_4][0] = 1.11299784711;
    _collision_strength[CIII][TRANSITION_1_to_4][1] = -1.38231571261e-05;
    _collision_strength[CIII][TRANSITION_1_to_4][2] = 0.17974425361;
    _collision_strength[CIII][TRANSITION_1_to_4][3] = -0.114824292495;
    _collision_strength[CIII][TRANSITION_1_to_4][4] = 1.15358738505e-06;
    _collision_strength[CIII][TRANSITION_1_to_4][5] = 0.114824292472;
    _collision_strength[CIII][TRANSITION_1_to_4][6] = 1.00000000001;
    _collision_strength[CIII][TRANSITION_2_to_3][0] = 0.83743274479;
    _collision_strength[CIII][TRANSITION_2_to_3][1] = 0.00568976058875;
    _collision_strength[CIII][TRANSITION_2_to_3][2] = 0.450792629302;
    _collision_strength[CIII][TRANSITION_2_to_3][3] = 0.00290023704999;
    _collision_strength[CIII][TRANSITION_2_to_3][4] = -0.000506695387193;
    _collision_strength[CIII][TRANSITION_2_to_3][5] = -0.00290022076109;
    _collision_strength[CIII][TRANSITION_2_to_3][6] = 1.00000037589;
    _collision_strength[CIII][TRANSITION_2_to_4][0] = 1.11290316385;
    _collision_strength[CIII][TRANSITION_2_to_4][1] = -4.1482086006e-05;
    _collision_strength[CIII][TRANSITION_2_to_4][2] = 0.53966645149;
    _collision_strength[CIII][TRANSITION_2_to_4][3] = 0.0301543484698;
    _collision_strength[CIII][TRANSITION_2_to_4][4] = 3.46165385743e-06;
    _collision_strength[CIII][TRANSITION_2_to_4][5] = -0.0301543485388;
    _collision_strength[CIII][TRANSITION_2_to_4][6] = 0.99999999985;
    _collision_strength[CIII][TRANSITION_3_to_4][0] = 1.11294049421;
    _collision_strength[CIII][TRANSITION_3_to_4][1] = -6.91266088737e-05;
    _collision_strength[CIII][TRANSITION_3_to_4][2] = 0.899155762243;
    _collision_strength[CIII][TRANSITION_3_to_4][3] = -0.0189378233178;
    _collision_strength[CIII][TRANSITION_3_to_4][4] = 5.76865382713e-06;
    _collision_strength[CIII][TRANSITION_3_to_4][5] = 0.0189378232029;
    _collision_strength[CIII][TRANSITION_3_to_4][6] = 1.0000000004;
  }

  /// two level elements
  // note that we need to remap the indices of these levels, as the enum values
  // start counting from LINECOOLINGDATA_NUMFIVELEVELELEMENTS
  // this is better than needing to do the remapping outside of the class, which
  // is what would happen otherwise...

  /// NIII
  {
    // Blum & Pradhan (1992), table 5, first energy level (in Ry)
    // ground state: 2P1/2
    // excited state: 2P3/2
    _two_level_element_data[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_ENERGY_DIFFERENCE] =
                               0.00159 * Ry_over_k;
    // Galavis, Mendoza & Zeippen (1998), table 4, 1 to 2 transition (in s^-1)
    _two_level_element_data[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_TRANSITION_PROBABILITY] = 4.736e-5;
    // Blum & Pradhan (1992), table 3, value for 10,000 K, 1 to 2 transition
    _two_level_element_data[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_COLLISION_STRENGTH] = 1.4454;
    // statistical weights: level 0 is a P_{1/2} level, while level 1 is a
    // P_{3/2}
    _two_level_element_data[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0] = 0.5;
    _two_level_element_data[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1] = 0.25;
  }

  /// NeII
  {
    // Saraph & Tully (1994), table 2, fine structure splitting energy for
    // Z = 10 (in Ry)
    // ground state: 2P3/2
    // excited state: 2P1/2
    _two_level_element_data[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_ENERGY_DIFFERENCE] =
                               0.0071 * Ry_over_k;
    // Kaufman & Sugar (1986), table 7 (in s^-1)
    _two_level_element_data[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_TRANSITION_PROBABILITY] = 8.55e-3;
    // Griffin, Mitnik & Badnell (2001), table 4, value for 10,000 K
    _two_level_element_data[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_COLLISION_STRENGTH] = 0.314;
    // statistical weights: level 0 is a P_{3/2} level, while level 1 is a
    // P_{1/2}
    _two_level_element_data[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0] = 0.25;
    _two_level_element_data[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1] = 0.5;
  }

  /// numerical factors that are precomputed

  // Boltzmann constant (in J s^-1)
  const double kb =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
  // Planck constant (in J s)
  const double h =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  // electron mass (in kg)
  const double m_e =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_ELECTRON_MASS);

  _collision_strength_prefactor =
      h * h / (std::sqrt(kb) * std::pow(2. * M_PI * m_e, 1.5));
}

/**
 * @brief Get the transition probability for deexcitation for the given
 * transition of the given element.
 *
 * @param element LineCoolingDataFiveLevelElement.
 * @param transition LineCoolingDataTransition.
 * @return Transition probability for deexcitation (in s^-1).
 */
double LineCoolingData::get_transition_probability(
    LineCoolingDataFiveLevelElement element,
    LineCoolingDataTransition transition) const {
  return _transition_probability[element][transition];
}

/**
 * @brief Get the energy difference for the given transition of the given
 * element.
 *
 * @param element LineCoolingDataFiveLevelElement.
 * @param transition LineCoolingDataTransition.
 * @return Energy difference (in K).
 */
double LineCoolingData::get_energy_difference(
    LineCoolingDataFiveLevelElement element,
    LineCoolingDataTransition transition) const {
  return _energy_difference[element][transition];
}

/**
 * @brief Get the statistical weight for the given level of the given element.
 *
 * @param element LineCoolingDataFiveLevelElement.
 * @param level Level index.
 * @return Statistical weight.
 */
double
LineCoolingData::get_statistical_weight(LineCoolingDataFiveLevelElement element,
                                        unsigned char level) const {
  return 1. / _inverse_statistical_weight[element][level];
}

/**
 * @brief Solve a system of 5 coupled linear equations.
 *
 * We assume a system of equations of the form
 * @f[
 *   A.X = B,
 * @f]
 * with
 * @f[
 * A = \begin{pmatrix}
 * A00 & A01 & A02 & A03 & A04 \\
 * A10 & A11 & A12 & A13 & A14 \\
 * A20 & A21 & A22 & A23 & A24 \\
 * A30 & A31 & A32 & A33 & A34 \\
 * A40 & A41 & A42 & A43 & A44
 * \end{pmatrix},
 * @f]
 * and
 * @f[
 * B = \begin{pmatrix}
 * B0 \\
 * B1 \\
 * B2 \\
 * B3 \\
 * B4
 * \end{pmatrix}.
 * @f]
 * We want to solve for the elements of the column matrix @f$X@f$.
 *
 * To do this, we first reduce the combined matrix @f$AB@f$ (the matrix @f$A@f$
 * with the matrix @f$B@f$ added as an extra column) to an upper triangular
 * matrix using Gaussian elimination. This automatically gives us the solution
 * for the last element of @f$X@f$. We then use this element to recursively
 * solve for the others by substitution.
 *
 * To avoid division by zero, we always use the next row with the largest
 * coefficient, and interchange rows if necessary. This means that to eliminate
 * e.g. row 2, we first check which of the rows 2-5 contains the largest value
 * in column 2. Suppose this is row 4. We then interchange rows 2 and 4, and
 * divide all columns of the new row 2 (the original row 4) by the value in
 * column 2 of row 2. Row 2 now contains a 1 at position 2, which is what we
 * want. We then use the values in row 2 to make sure all elements up to column
 * 2 are zero in all rows below row 2. Since row 2, column 2 is 1, this can be
 * achieved by simply subtracting row 2, column i multiplied by row j, column 2
 * from row j, column i for all i and j larger than 2. We do not actually do the
 * calculation for column 2 and smaller, since we know the result is zero
 * (column 1 was already zero after the elimination of row 1).
 *
 * Note that both matrix @f$A@f$ and matrix @f$B@f$ are modified in place. When
 * the method returns, @f$B@f$ contains the elements of the matrix @f$X@f$.
 *
 * @param A Elements of the matrix @f$A@f$.
 * @param B Elements of the matrix @f$B@f$, and elements of the solution on
 * exit.
 * @return Exit code: 0 on success. If a non zero value is returned, the values
 * stored in B on exit are meaningless and should not be used.
 */
int LineCoolingData::simq(double A[5][5], double B[5]) {

  for (unsigned char j = 0; j < 5; ++j) {
    // find the next row with the largest coefficient
    unsigned char imax = 0;
    double Amax = 0.;
    for (unsigned char i = j; i < 5; ++i) {
      if (std::abs(A[i][j]) > std::abs(Amax)) {
        Amax = A[i][j];
        imax = i;
      }
    }
    // check that the matrix is non-singular
    if (Amax == 0.) {
      return 1;
    }
    const double Amax_inv = 1. / Amax;
    // imax now contains the index of the row with the largest coefficient
    // interchange rows if necessary to make sure that the row with the largest
    // coefficient is the next row to eliminate with
    for (unsigned char k = 0; k < 5; ++k) {
      if (imax != j) {
        const double save = A[j][k];
        A[j][k] = A[imax][k];
        A[imax][k] = save;
      }
      A[j][k] *= Amax_inv;
    }
    if (imax != j) {
      const double save = B[j];
      B[j] = B[imax];
      B[imax] = save;
    }
    B[j] *= Amax_inv;
    // row j now contains a 1 at position j
    // all elements in columns with indices smaller than j are supposed to be
    // zero due to previous eliminations (we do not set them to zero however to
    // save computations)
    if (j < 4) {
      // we are not finished yet: use row j to eliminate all rows below row j
      for (unsigned char i = j + 1; i < 5; ++i) {
        for (unsigned char k = j + 1; k < 5; ++k) {
          A[i][k] -= A[i][j] * A[j][k];
        }
        B[i] -= A[i][j] * B[j];
      }
    }
  }
  // the matrix now has the form
  //  1  A01 A02 A03 A04 B0
  //  0   1  A12 A13 A14 B1
  //  0   0   1  A23 A24 B2
  //  0   0   0   1  A34 B3
  //  0   0   0   0   1  B4
  // In other words: B[4] contains the value of the last variable
  // use it to recursively solve for the others
  for (unsigned char i = 0; i < 4; ++i) {
    for (unsigned char j = 0; j < i + 1; ++j) {
      B[3 - i] -= B[4 - j] * A[3 - i][4 - j];
    }
  }
  return 0;
}

/**
 * @brief Find the level populations for the given element at the given
 * temperature.
 *
 * @param element LineCoolingDataFiveLevelElement.
 * @param collision_strength_prefactor Prefactor for the collision strengths
 * (in s^-1).
 * @param T Temperature (in K).
 * @param Tinv Inverse of the temperature (in K^-1).
 * @param logT Logarithm of the temperature in K.
 * @param level_populations Array to store the resulting level populations in.
 */
void LineCoolingData::compute_level_populations(
    LineCoolingDataFiveLevelElement element,
    double collision_strength_prefactor, double T, double Tinv, double logT,
    double level_populations[5]) const {

  double level_matrix[5][5];
  // initialize the level populations and the first row of the coefficient
  // matrix
  for (unsigned char i = 0; i < 5; ++i) {
    level_matrix[0][i] = 1.;
    level_populations[i] = 0.;
  }
  // the first row of the coefficient matrix expresses the constant number of
  // particles: the sum of all level populations is unity
  level_populations[0] = 1.;

  // precompute the collision rates for the given temperature
  double collision_rate_down[NUMBER_OF_TRANSITIONS];
  double collision_rate_up[NUMBER_OF_TRANSITIONS];
  for (int i = 0; i < NUMBER_OF_TRANSITIONS; ++i) {
    const double collision_strength =
        collision_strength_prefactor *
        std::pow(T, _collision_strength[element][i][0]) *
        (_collision_strength[element][i][1] +
         _collision_strength[element][i][2] * Tinv +
         _collision_strength[element][i][3] * T +
         _collision_strength[element][i][4] * logT +
         _collision_strength[element][i][5] *
             std::pow(T, _collision_strength[element][i][6]));
    collision_rate_down[i] = collision_strength;
    collision_rate_up[i] =
        collision_strength * std::exp(-_energy_difference[element][i] * Tinv);
  }

  level_matrix[1][0] = collision_rate_up[TRANSITION_0_to_1] *
                       _inverse_statistical_weight[element][0];
  level_matrix[1][1] = -(_transition_probability[element][TRANSITION_0_to_1] +
                         _inverse_statistical_weight[element][1] *
                             (collision_rate_down[TRANSITION_0_to_1] +
                              collision_rate_up[TRANSITION_1_to_2] +
                              collision_rate_up[TRANSITION_1_to_3] +
                              collision_rate_up[TRANSITION_1_to_4]));
  level_matrix[1][2] = _transition_probability[element][TRANSITION_1_to_2] +
                       _inverse_statistical_weight[element][2] *
                           collision_rate_down[TRANSITION_1_to_2];
  level_matrix[1][3] = _transition_probability[element][TRANSITION_1_to_3] +
                       _inverse_statistical_weight[element][3] *
                           collision_rate_down[TRANSITION_1_to_3];
  level_matrix[1][4] = _transition_probability[element][TRANSITION_1_to_4] +
                       _inverse_statistical_weight[element][4] *
                           collision_rate_down[TRANSITION_1_to_4];

  level_matrix[2][0] = collision_rate_up[TRANSITION_0_to_2] *
                       _inverse_statistical_weight[element][0];
  level_matrix[2][1] = collision_rate_up[TRANSITION_1_to_2] *
                       _inverse_statistical_weight[element][1];
  level_matrix[2][2] = -(_transition_probability[element][TRANSITION_0_to_2] +
                         _transition_probability[element][TRANSITION_1_to_2] +
                         _inverse_statistical_weight[element][2] *
                             (collision_rate_down[TRANSITION_0_to_2] +
                              collision_rate_down[TRANSITION_1_to_2] +
                              collision_rate_up[TRANSITION_2_to_3] +
                              collision_rate_up[TRANSITION_2_to_4]));
  level_matrix[2][3] = _transition_probability[element][TRANSITION_2_to_3] +
                       collision_rate_down[TRANSITION_2_to_3] *
                           _inverse_statistical_weight[element][3];
  level_matrix[2][4] = _transition_probability[element][TRANSITION_2_to_4] +
                       collision_rate_down[TRANSITION_2_to_4] *
                           _inverse_statistical_weight[element][4];

  level_matrix[3][0] = collision_rate_up[TRANSITION_0_to_3] *
                       _inverse_statistical_weight[element][0];
  level_matrix[3][1] = collision_rate_up[TRANSITION_1_to_3] *
                       _inverse_statistical_weight[element][1];
  level_matrix[3][2] = collision_rate_up[TRANSITION_2_to_3] *
                       _inverse_statistical_weight[element][2];
  level_matrix[3][3] = -(_transition_probability[element][TRANSITION_0_to_3] +
                         _transition_probability[element][TRANSITION_1_to_3] +
                         _transition_probability[element][TRANSITION_2_to_3] +
                         _inverse_statistical_weight[element][3] *
                             (collision_rate_down[TRANSITION_0_to_3] +
                              collision_rate_down[TRANSITION_1_to_3] +
                              collision_rate_down[TRANSITION_2_to_3] +
                              collision_rate_up[TRANSITION_3_to_4]));
  level_matrix[3][4] = _transition_probability[element][TRANSITION_3_to_4] +
                       collision_rate_down[TRANSITION_3_to_4] *
                           _inverse_statistical_weight[element][4];

  level_matrix[4][0] = collision_rate_up[TRANSITION_0_to_4] *
                       _inverse_statistical_weight[element][0];
  level_matrix[4][1] = collision_rate_up[TRANSITION_1_to_4] *
                       _inverse_statistical_weight[element][1];
  level_matrix[4][2] = collision_rate_up[TRANSITION_2_to_4] *
                       _inverse_statistical_weight[element][2];
  level_matrix[4][3] = collision_rate_up[TRANSITION_3_to_4] *
                       _inverse_statistical_weight[element][3];
  level_matrix[4][4] = -(_transition_probability[element][TRANSITION_0_to_4] +
                         _transition_probability[element][TRANSITION_1_to_4] +
                         _transition_probability[element][TRANSITION_2_to_4] +
                         _transition_probability[element][TRANSITION_3_to_4] +
                         _inverse_statistical_weight[element][4] *
                             (collision_rate_down[TRANSITION_0_to_4] +
                              collision_rate_down[TRANSITION_1_to_4] +
                              collision_rate_down[TRANSITION_2_to_4] +
                              collision_rate_down[TRANSITION_3_to_4]));

  // find level populations
  const int status = simq(level_matrix, level_populations);
  if (status != 0) {
    // something went wrong
    cmac_error("Singular matrix in level population computation (element: %i, "
               "T: %g, n_e: %g)!",
               element, 1. / Tinv,
               collision_strength_prefactor /
                   (_collision_strength_prefactor * std::sqrt(Tinv)));
  }
}

/**
 * @brief Find the level population of the second level for the given two level
 * element at the given temperature.
 *
 * @param element LineCoolingDataTwoLevelElement.
 * @param collision_strength_prefactor Prefactor for the collision strengths
 * (in s^-1).
 * @param Tinv Inverse temperature (in K^-1).
 * @return Level population of the second level.
 */
double LineCoolingData::compute_level_population(
    LineCoolingDataTwoLevelElement element, double collision_strength_prefactor,
    double Tinv) const {

  // note that we need to remap the element index
  const int i = element - LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  const double ksi =
      _two_level_element_data[i][TWOLEVELFIELD_ENERGY_DIFFERENCE];
  const double A =
      _two_level_element_data[i][TWOLEVELFIELD_TRANSITION_PROBABILITY];
  const double Gamma =
      _two_level_element_data[i][TWOLEVELFIELD_COLLISION_STRENGTH];
  const double inv_omega_1 =
      _two_level_element_data[i][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0];
  const double inv_omega_2 =
      _two_level_element_data[i][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1];
  const double Texp = std::exp(-ksi * Tinv);
  return collision_strength_prefactor * Gamma * Texp * inv_omega_1 /
         (A +
          collision_strength_prefactor * Gamma *
              (inv_omega_2 + Texp * inv_omega_1));
}

/**
 * @brief Get the radiative energy losses due to line cooling at the given
 * temperature, electron density and coolant abundances.
 *
 * We consider 12 ions; 10 with 5 low lying collisionally excited levels, and 2
 * with only 2 levels. For the former, we solve equations (3.27) and (3.28) in
 * Osterbrock & Ferland (2006), with the cooling rate given by equation (3.29).
 *
 * For the latter, we use equations (3.24) and (3.25), together with the total
 * number of ions:
 *\f[
 * n_1 + n_2 = n.
 * \f]
 * Equation (3.24) gives a relation between the level populations of ion 1 and
 * 2:
 * \f[
 * \frac{n_2}{n_1} = f,
 * \f]
 * which together with the total number of ions leads to
 * \f[
 * n_2 = \frac{f}{1+f}n.
 * \f]
 * We substitute this in equation (3.25).
 *
 * @param temperature Temperature (in K).
 * @param electron_density Electron density (in m^-3).
 * @param abundances Abdunances of coolants.
 * @return Radiative cooling per hydrogen atom (in kg m^2s^-3).
 */
double LineCoolingData::get_cooling(double temperature, double electron_density,
                                    const double *abundances) const {

  if (electron_density == 0.) {
    // we cannot return a 0 cooling rate, because that crashes our iterative
    // temperature finding scheme
    return 1.e-99;
  }

  /// initialize some variables

  // Boltzmann constant (in J s^-1)
  const double kb =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

  // _collision_strength_prefactor has units K^0.5 m^3 s^-1;
  // collision_strength_prefactor has units s^-1
  const double collision_strength_prefactor =
      _collision_strength_prefactor * electron_density / std::sqrt(temperature);
  const double Tinv = 1. / temperature;
  const double logT = std::log(temperature);

  /// five level elements

  double cooling = 0.;
  for (int j = 0; j < LINECOOLINGDATA_NUMFIVELEVELELEMENTS; ++j) {

    const LineCoolingDataFiveLevelElement element =
        static_cast< LineCoolingDataFiveLevelElement >(j);

    double level_populations[5];
    compute_level_populations(element, collision_strength_prefactor,
                              temperature, Tinv, logT, level_populations);

    // compute the cooling for each transition
    // this corresponds to equation (3.29) in Osterbrock & Ferland (2006)
    const double cl2 = level_populations[1] *
                       _transition_probability[j][TRANSITION_0_to_1] *
                       _energy_difference[j][TRANSITION_0_to_1];
    const double cl3 =
        level_populations[2] * (_transition_probability[j][TRANSITION_0_to_2] *
                                    _energy_difference[j][TRANSITION_0_to_2] +
                                _transition_probability[j][TRANSITION_1_to_2] *
                                    _energy_difference[j][TRANSITION_1_to_2]);
    const double cl4 =
        level_populations[3] * (_transition_probability[j][TRANSITION_0_to_3] *
                                    _energy_difference[j][TRANSITION_0_to_3] +
                                _transition_probability[j][TRANSITION_1_to_3] *
                                    _energy_difference[j][TRANSITION_1_to_3] +
                                _transition_probability[j][TRANSITION_2_to_3] *
                                    _energy_difference[j][TRANSITION_2_to_3]);
    const double cl5 =
        level_populations[4] * (_transition_probability[j][TRANSITION_0_to_4] *
                                    _energy_difference[j][TRANSITION_0_to_4] +
                                _transition_probability[j][TRANSITION_1_to_4] *
                                    _energy_difference[j][TRANSITION_1_to_4] +
                                _transition_probability[j][TRANSITION_2_to_4] *
                                    _energy_difference[j][TRANSITION_2_to_4] +
                                _transition_probability[j][TRANSITION_3_to_4] *
                                    _energy_difference[j][TRANSITION_3_to_4]);

    cooling += abundances[j] * kb * (cl2 + cl3 + cl4 + cl5);
  }

  /// 2 level atoms

  // offset of two level elements in the abundances array
  const int offset = LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  for (int i = offset; i < LINECOOLINGDATA_NUMTWOLEVELELEMENTS; ++i) {

    const int index = i - offset;
    const LineCoolingDataTwoLevelElement element =
        static_cast< LineCoolingDataTwoLevelElement >(i);
    const double level_population =
        compute_level_population(element, collision_strength_prefactor, Tinv);
    cooling +=
        abundances[i] * kb *
        _two_level_element_data[index][TWOLEVELFIELD_ENERGY_DIFFERENCE] *
        _two_level_element_data[index][TWOLEVELFIELD_TRANSITION_PROBABILITY] *
        level_population;
  }

  return cooling;
}

/**
 * @brief Calculate the strength of all emission lines for which we have data.
 *
 * @param temperature Temperature (in K).
 * @param electron_density Electron density (in m^-3).
 * @param abundances Ion abundances.
 * @return std::vector containing, for each ion, a std::vector of line strengths
 * for each transition (in J s^-1).
 */
std::vector< std::vector< double > >
LineCoolingData::get_line_strengths(double temperature, double electron_density,
                                    const double *abundances) const {

  /// initialize some variables

  // Boltzmann constant (in J s^-1)
  const double kb =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

  const double collision_strength_prefactor =
      _collision_strength_prefactor * electron_density / std::sqrt(temperature);
  const double Tinv = 1. / temperature;
  const double logT = std::log(temperature);

  // vector to store line strengths in
  std::vector< std::vector< double > > line_strengths(
      LINECOOLINGDATA_NUMFIVELEVELELEMENTS +
      LINECOOLINGDATA_NUMTWOLEVELELEMENTS);

  /// 5 level elements

  // there are no lines for element 0 (NI), so we skip the iteration for that
  // element
  for (int j = 1; j < 10; ++j) {

    line_strengths[j].resize(NUMBER_OF_TRANSITIONS);

    const LineCoolingDataFiveLevelElement element =
        static_cast< LineCoolingDataFiveLevelElement >(j);

    double level_populations[5];
    compute_level_populations(element, collision_strength_prefactor,
                              temperature, Tinv, logT, level_populations);

    const double prefactor = abundances[j] * kb;

    line_strengths[j][TRANSITION_0_to_1] =
        prefactor * level_populations[1] *
        _transition_probability[j][TRANSITION_0_to_1] *
        _energy_difference[j][TRANSITION_0_to_1];
    line_strengths[j][TRANSITION_0_to_2] =
        prefactor * level_populations[2] *
        _transition_probability[j][TRANSITION_0_to_2] *
        _energy_difference[j][TRANSITION_0_to_2];
    line_strengths[j][TRANSITION_1_to_2] =
        prefactor * level_populations[2] *
        _transition_probability[j][TRANSITION_1_to_2] *
        _energy_difference[j][TRANSITION_1_to_2];
    line_strengths[j][TRANSITION_0_to_3] =
        prefactor * level_populations[3] *
        _transition_probability[j][TRANSITION_0_to_3] *
        _energy_difference[j][TRANSITION_0_to_3];
    line_strengths[j][TRANSITION_1_to_3] =
        prefactor * level_populations[3] *
        _transition_probability[j][TRANSITION_1_to_3] *
        _energy_difference[j][TRANSITION_1_to_3];
    line_strengths[j][TRANSITION_2_to_3] =
        prefactor * level_populations[3] *
        _transition_probability[j][TRANSITION_2_to_3] *
        _energy_difference[j][TRANSITION_2_to_3];
    line_strengths[j][TRANSITION_0_to_4] =
        prefactor * level_populations[4] *
        _transition_probability[j][TRANSITION_0_to_4] *
        _energy_difference[j][TRANSITION_0_to_4];
    line_strengths[j][TRANSITION_1_to_4] =
        prefactor * level_populations[4] *
        _transition_probability[j][TRANSITION_1_to_4] *
        _energy_difference[j][TRANSITION_1_to_4];
    line_strengths[j][TRANSITION_2_to_4] =
        prefactor * level_populations[4] *
        _transition_probability[j][TRANSITION_2_to_4] *
        _energy_difference[j][TRANSITION_2_to_4];
    line_strengths[j][TRANSITION_3_to_4] =
        prefactor * level_populations[4] *
        _transition_probability[j][TRANSITION_3_to_4] *
        _energy_difference[j][TRANSITION_3_to_4];
  }

  /// 2 level elements

  // offset of two level elements in the abundances array
  const int offset = LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  for (int i = offset; i < LINECOOLINGDATA_NUMTWOLEVELELEMENTS; ++i) {

    const int index = i - offset;

    line_strengths[i].resize(1);

    const LineCoolingDataTwoLevelElement element =
        static_cast< LineCoolingDataTwoLevelElement >(i);

    const double level_population =
        compute_level_population(element, collision_strength_prefactor, Tinv);

    line_strengths[i][0] =
        abundances[i] * kb * level_population *
        _two_level_element_data[index][TWOLEVELFIELD_ENERGY_DIFFERENCE] *
        _two_level_element_data[index][TWOLEVELFIELD_TRANSITION_PROBABILITY];
  }

  return line_strengths;
}
