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
    // these data could do with better fits...
    _collision_strength[NI][TRANSITION_0_to_1] = 0.0262;
    _collision_strength_exponent[NI][TRANSITION_0_to_1] = 3.00361496607;
    _collision_strength[NI][TRANSITION_0_to_2] = 0.0175;
    _collision_strength_exponent[NI][TRANSITION_0_to_2] = 2.89422551909;
    _collision_strength[NI][TRANSITION_0_to_3] = 0.0105;
    _collision_strength_exponent[NI][TRANSITION_0_to_3] = 1.13420100466;
    _collision_strength[NI][TRANSITION_0_to_4] = 0.021;
    _collision_strength_exponent[NI][TRANSITION_0_to_4] = 1.07040488898;
    _collision_strength[NI][TRANSITION_1_to_2] = 3.24;
    _collision_strength_exponent[NI][TRANSITION_1_to_2] = 0.713317341817;
    _collision_strength[NI][TRANSITION_1_to_3] = 0.56;
    _collision_strength_exponent[NI][TRANSITION_1_to_3] = 2.93912918004;
    _collision_strength[NI][TRANSITION_1_to_4] = 4.15;
    _collision_strength_exponent[NI][TRANSITION_1_to_4] = 2.83603951829;
    _collision_strength[NI][TRANSITION_2_to_3] = 2.11;
    _collision_strength_exponent[NI][TRANSITION_2_to_3] = 2.83805595321;
    _collision_strength[NI][TRANSITION_2_to_4] = 1.07;
    _collision_strength_exponent[NI][TRANSITION_2_to_4] = 2.91341222593;
    _collision_strength[NI][TRANSITION_3_to_4] = 1.99;
    _collision_strength_exponent[NI][TRANSITION_3_to_4] = 1.63273153773;
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
    // validity estimated by eye, is approximately [5,000 K; 25,000 K] for all
    // curves (and definitely better than the values used in Kenny's code)
    _collision_strength[NII][TRANSITION_0_to_1] = 0.4076;
    _collision_strength_exponent[NII][TRANSITION_0_to_1] = 0.120890503443;
    _collision_strength[NII][TRANSITION_0_to_2] = 0.272;
    _collision_strength_exponent[NII][TRANSITION_0_to_2] = 0.206218494855;
    _collision_strength[NII][TRANSITION_0_to_3] = 0.293433333333;
    _collision_strength_exponent[NII][TRANSITION_0_to_3] = 0.0454386659465;
    _collision_strength[NII][TRANSITION_0_to_4] = 0.0325666666667;
    _collision_strength_exponent[NII][TRANSITION_0_to_4] = 0.0531252760254;
    _collision_strength[NII][TRANSITION_1_to_2] = 1.12;
    _collision_strength_exponent[NII][TRANSITION_1_to_2] = 0.167121267122;
    _collision_strength[NII][TRANSITION_1_to_3] = 0.8803;
    _collision_strength_exponent[NII][TRANSITION_1_to_3] = 0.0459250739371;
    _collision_strength[NII][TRANSITION_1_to_4] = 0.0977;
    _collision_strength_exponent[NII][TRANSITION_1_to_4] = 0.0455168342781;
    _collision_strength[NII][TRANSITION_2_to_3] = 1.46716666667;
    _collision_strength_exponent[NII][TRANSITION_2_to_3] = 0.0459250736643;
    _collision_strength[NII][TRANSITION_2_to_4] = 0.162833333333;
    _collision_strength_exponent[NII][TRANSITION_2_to_4] = 0.0427417408534;
    _collision_strength[NII][TRANSITION_3_to_4] = 0.8338;
    _collision_strength_exponent[NII][TRANSITION_3_to_4] = -0.204643342755;
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
    _collision_strength[OI][TRANSITION_0_to_1] = 0.106;
    _collision_strength_exponent[OI][TRANSITION_0_to_1] = 1.15115610038;
    _collision_strength[OI][TRANSITION_0_to_2] = 0.0321;
    _collision_strength_exponent[OI][TRANSITION_0_to_2] = 0.969736258456;
    _collision_strength[OI][TRANSITION_0_to_3] = 0.162777777778;
    _collision_strength_exponent[OI][TRANSITION_0_to_3] = 0.758309213819;
    _collision_strength[OI][TRANSITION_0_to_4] = 0.0179444444444;
    _collision_strength_exponent[OI][TRANSITION_0_to_4] = 0.764336719318;
    _collision_strength[OI][TRANSITION_1_to_2] = 0.0283;
    _collision_strength_exponent[OI][TRANSITION_1_to_2] = 1.65586135957;
    _collision_strength[OI][TRANSITION_1_to_3] = 0.0976666666667;
    _collision_strength_exponent[OI][TRANSITION_1_to_3] = 0.758309213743;
    _collision_strength[OI][TRANSITION_1_to_4] = 0.0107666666667;
    _collision_strength_exponent[OI][TRANSITION_1_to_4] = 0.764336719513;
    _collision_strength[OI][TRANSITION_2_to_3] = 0.0325555555556;
    _collision_strength_exponent[OI][TRANSITION_2_to_3] = 0.758309213775;
    _collision_strength[OI][TRANSITION_2_to_4] = 0.00358888888889;
    _collision_strength_exponent[OI][TRANSITION_2_to_4] = 0.764336719235;
    _collision_strength[OI][TRANSITION_3_to_4] = 0.0883;
    _collision_strength_exponent[OI][TRANSITION_3_to_4] = 0.577711932569;
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
    // much better fitting is necessary...
    _collision_strength[OII][TRANSITION_0_to_1] = 0.834;
    _collision_strength_exponent[OII][TRANSITION_0_to_1] = 0.00228239318471;
    _collision_strength[OII][TRANSITION_0_to_2] = 0.554;
    _collision_strength_exponent[OII][TRANSITION_0_to_2] = -0.000138326600127;
    _collision_strength[OII][TRANSITION_0_to_3] = 0.256;
    _collision_strength_exponent[OII][TRANSITION_0_to_3] = 0.0198710452891;
    _collision_strength[OII][TRANSITION_0_to_4] = 0.132;
    _collision_strength_exponent[OII][TRANSITION_0_to_4] = 0.0183736539042;
    _collision_strength[OII][TRANSITION_1_to_2] = 1.203;
    _collision_strength_exponent[OII][TRANSITION_1_to_2] = -0.0309381910145;
    _collision_strength[OII][TRANSITION_1_to_3] = 0.851;
    _collision_strength_exponent[OII][TRANSITION_1_to_3] = 0.0395198992;
    _collision_strength[OII][TRANSITION_1_to_4] = 0.339;
    _collision_strength_exponent[OII][TRANSITION_1_to_4] = 0.04037113288;
    _collision_strength[OII][TRANSITION_2_to_3] = 0.472;
    _collision_strength_exponent[OII][TRANSITION_2_to_3] = 0.040966740518;
    _collision_strength[OII][TRANSITION_2_to_4] = 0.331;
    _collision_strength_exponent[OII][TRANSITION_2_to_4] = 0.0388864515812;
    _collision_strength[OII][TRANSITION_3_to_4] = 0.285;
    _collision_strength_exponent[OII][TRANSITION_3_to_4] = 0.0223104241703;
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
    // validity estimated by eye, is approximately [5,000 K; 25,000 K] for all
    // curves (and definitely better than the values used in Kenny's code)
    _collision_strength[OIII][TRANSITION_0_to_1] = 0.5454;
    _collision_strength_exponent[OIII][TRANSITION_0_to_1] = 0.0583811855633;
    _collision_strength[OIII][TRANSITION_0_to_2] = 0.2713;
    _collision_strength_exponent[OIII][TRANSITION_0_to_2] = 0.0893900923473;
    _collision_strength[OIII][TRANSITION_0_to_3] = 0.254355555556;
    _collision_strength_exponent[OIII][TRANSITION_0_to_3] = 0.138947107346;
    _collision_strength[OIII][TRANSITION_0_to_4] = 0.0325;
    _collision_strength_exponent[OIII][TRANSITION_0_to_4] = 0.141146060784;
    _collision_strength[OIII][TRANSITION_1_to_2] = 1.291;
    _collision_strength_exponent[OIII][TRANSITION_1_to_2] = 0.0731403594743;
    _collision_strength[OIII][TRANSITION_1_to_3] = 0.763066666667;
    _collision_strength_exponent[OIII][TRANSITION_1_to_3] = 0.138947108288;
    _collision_strength[OIII][TRANSITION_1_to_4] = 0.0975;
    _collision_strength_exponent[OIII][TRANSITION_1_to_4] = 0.159858565366;
    _collision_strength[OIII][TRANSITION_2_to_3] = 1.27177777778;
    _collision_strength_exponent[OIII][TRANSITION_2_to_3] = 0.13894710933;
    _collision_strength[OIII][TRANSITION_2_to_4] = 0.1625;
    _collision_strength_exponent[OIII][TRANSITION_2_to_4] = 0.159858567059;
    _collision_strength[OIII][TRANSITION_3_to_4] = 0.5815;
    _collision_strength_exponent[OIII][TRANSITION_3_to_4] = 0.157486912474;
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
    // these data values are not fitted properly with the proposed curve, and a
    // different fitting strategy will be necessary to fit them...
    _collision_strength[NeIII][TRANSITION_0_to_1] = 0.774;
    _collision_strength_exponent[NeIII][TRANSITION_0_to_1] = 0.0360460427437;
    _collision_strength[NeIII][TRANSITION_0_to_2] = 0.208;
    _collision_strength_exponent[NeIII][TRANSITION_0_to_2] = 0.0206794792523;
    _collision_strength[NeIII][TRANSITION_0_to_3] = 0.754;
    _collision_strength_exponent[NeIII][TRANSITION_0_to_3] = -0.0158496002015;
    _collision_strength[NeIII][TRANSITION_0_to_4] = 0.084;
    _collision_strength_exponent[NeIII][TRANSITION_0_to_4] = 0.0155831445199;
    _collision_strength[NeIII][TRANSITION_1_to_2] = 0.244;
    _collision_strength_exponent[NeIII][TRANSITION_1_to_2] = 0.0527279729302;
    _collision_strength[NeIII][TRANSITION_1_to_3] = 0.452;
    _collision_strength_exponent[NeIII][TRANSITION_1_to_3] = -0.0149474243201;
    _collision_strength[NeIII][TRANSITION_1_to_4] = 0.05;
    _collision_strength_exponent[NeIII][TRANSITION_1_to_4] = 0.0223434756313;
    _collision_strength[NeIII][TRANSITION_2_to_3] = 0.151;
    _collision_strength_exponent[NeIII][TRANSITION_2_to_3] = -0.0113458238066;
    _collision_strength[NeIII][TRANSITION_2_to_4] = 0.017;
    _collision_strength_exponent[NeIII][TRANSITION_2_to_4] = 0.0321061815655;
    _collision_strength[NeIII][TRANSITION_3_to_4] = 0.269;
    _collision_strength_exponent[NeIII][TRANSITION_3_to_4] = 0.0408073833617;
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
    // the fits are reasonable in this case
    _collision_strength[SII][TRANSITION_0_to_1] = 2.56;
    _collision_strength_exponent[SII][TRANSITION_0_to_1] = -0.0940407351461;
    _collision_strength[SII][TRANSITION_0_to_2] = 3.83;
    _collision_strength_exponent[SII][TRANSITION_0_to_2] = -0.0933764679975;
    _collision_strength[SII][TRANSITION_0_to_3] = 0.704;
    _collision_strength_exponent[SII][TRANSITION_0_to_3] = 0.042753148509;
    _collision_strength[SII][TRANSITION_0_to_4] = 1.42;
    _collision_strength_exponent[SII][TRANSITION_0_to_4] = 0.0381510814487;
    _collision_strength[SII][TRANSITION_1_to_2] = 6.89;
    _collision_strength_exponent[SII][TRANSITION_1_to_2] = -0.121018378151;
    _collision_strength[SII][TRANSITION_1_to_3] = 1.47;
    _collision_strength_exponent[SII][TRANSITION_1_to_3] = 0.0155625859981;
    _collision_strength[SII][TRANSITION_1_to_4] = 2.39;
    _collision_strength_exponent[SII][TRANSITION_1_to_4] = -0.0044941415539;
    _collision_strength[SII][TRANSITION_2_to_3] = 1.78;
    _collision_strength_exponent[SII][TRANSITION_2_to_3] = -0.0145851038874;
    _collision_strength[SII][TRANSITION_2_to_4] = 4.06;
    _collision_strength_exponent[SII][TRANSITION_2_to_4] = 0.00459345784943;
    _collision_strength[SII][TRANSITION_3_to_4] = 1.8;
    _collision_strength_exponent[SII][TRANSITION_3_to_4] = 0.0302135311733;
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
    // the fits are reasonable in this case
    _collision_strength[SIII][TRANSITION_0_to_1] = 2.26;
    _collision_strength_exponent[SIII][TRANSITION_0_to_1] = -0.0980735715102;
    _collision_strength[SIII][TRANSITION_0_to_2] = 1.02;
    _collision_strength_exponent[SIII][TRANSITION_0_to_2] = 0.187197754794;
    _collision_strength[SIII][TRANSITION_0_to_3] = 0.729;
    _collision_strength_exponent[SIII][TRANSITION_0_to_3] = 0.0646941085028;
    _collision_strength[SIII][TRANSITION_0_to_4] = 0.125;
    _collision_strength_exponent[SIII][TRANSITION_0_to_4] = 0.22673505046;
    _collision_strength[SIII][TRANSITION_1_to_2] = 5.1;
    _collision_strength_exponent[SIII][TRANSITION_1_to_2] = 0.0383455125877;
    _collision_strength[SIII][TRANSITION_1_to_3] = 2.17;
    _collision_strength_exponent[SIII][TRANSITION_1_to_3] = 0.0676576085243;
    _collision_strength[SIII][TRANSITION_1_to_4] = 0.33;
    _collision_strength_exponent[SIII][TRANSITION_1_to_4] = 0.156636569012;
    _collision_strength[SIII][TRANSITION_2_to_3] = 4.02;
    _collision_strength_exponent[SIII][TRANSITION_2_to_3] = 0.0513491293562;
    _collision_strength[SIII][TRANSITION_2_to_4] = 0.545;
    _collision_strength_exponent[SIII][TRANSITION_2_to_4] = 0.155917865017;
    _collision_strength[SIII][TRANSITION_3_to_4] = 1.38;
    _collision_strength_exponent[SIII][TRANSITION_3_to_4] = 0.255995735614;
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
    // these data values are not fitted properly with the proposed curve, and a
    // different fitting strategy will be necessary to fit them...
    _collision_strength[CII][TRANSITION_0_to_1] = 2.028;
    _collision_strength_exponent[CII][TRANSITION_0_to_1] = 0.0839065230943;
    _collision_strength[CII][TRANSITION_0_to_2] = 0.261;
    _collision_strength_exponent[CII][TRANSITION_0_to_2] = -0.00341736535797;
    _collision_strength[CII][TRANSITION_0_to_3] = 0.392;
    _collision_strength_exponent[CII][TRANSITION_0_to_3] = -0.00288752888829;
    _collision_strength[CII][TRANSITION_0_to_4] = 0.247;
    _collision_strength_exponent[CII][TRANSITION_0_to_4] = 0.0286538848796;
    _collision_strength[CII][TRANSITION_1_to_2] = 0.181;
    _collision_strength_exponent[CII][TRANSITION_1_to_2] = 0.0401915985567;
    _collision_strength[CII][TRANSITION_1_to_3] = 0.501;
    _collision_strength_exponent[CII][TRANSITION_1_to_3] = 0.0200040684801;
    _collision_strength[CII][TRANSITION_1_to_4] = 1.111;
    _collision_strength_exponent[CII][TRANSITION_1_to_4] = -0.0065220578966;
    _collision_strength[CII][TRANSITION_2_to_3] = 0.792;
    _collision_strength_exponent[CII][TRANSITION_2_to_3] = 0.285318320282;
    _collision_strength[CII][TRANSITION_2_to_4] = 0.836;
    _collision_strength_exponent[CII][TRANSITION_2_to_4] = 0.179371908326;
    _collision_strength[CII][TRANSITION_3_to_4] = 1.926;
    _collision_strength_exponent[CII][TRANSITION_3_to_4] = 0.225302392269;
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
    // these data values are not fitted properly with the proposed curve, and a
    // different fitting strategy will be necessary to fit them...
    _collision_strength[CIII][TRANSITION_0_to_1] = 0.115555555556;
    _collision_strength_exponent[CIII][TRANSITION_0_to_1] = -0.12535308605;
    _collision_strength[CIII][TRANSITION_0_to_2] = 0.346666666667;
    _collision_strength_exponent[CIII][TRANSITION_0_to_2] = -0.12535308689;
    _collision_strength[CIII][TRANSITION_0_to_3] = 0.577777777778;
    _collision_strength_exponent[CIII][TRANSITION_0_to_3] = -0.125353087003;
    _collision_strength[CIII][TRANSITION_0_to_4] = 4.0;
    _collision_strength_exponent[CIII][TRANSITION_0_to_4] = 0.118185461614;
    _collision_strength[CIII][TRANSITION_1_to_2] = 0.962;
    _collision_strength_exponent[CIII][TRANSITION_1_to_2] = -0.051083603727;
    _collision_strength[CIII][TRANSITION_1_to_3] = 0.718;
    _collision_strength_exponent[CIII][TRANSITION_1_to_3] = 0.0535157399273;
    _collision_strength[CIII][TRANSITION_1_to_4] = 0.411111111111;
    _collision_strength_exponent[CIII][TRANSITION_1_to_4] = -0.227262416038;
    _collision_strength[CIII][TRANSITION_2_to_3] = 2.78;
    _collision_strength_exponent[CIII][TRANSITION_2_to_3] = 0.0138404451029;
    _collision_strength[CIII][TRANSITION_2_to_4] = 1.23333333333;
    _collision_strength_exponent[CIII][TRANSITION_2_to_4] = -0.227262416396;
    _collision_strength[CIII][TRANSITION_3_to_4] = 2.05555555556;
    _collision_strength_exponent[CIII][TRANSITION_3_to_4] = -0.227262416164;
  }

  /// two level elements

  /// NIII
  {
    // Blum & Pradhan (1992), table 5, first energy level (in Ry)
    // ground state: 2P1/2
    // excited state: 2P3/2
    _two_level_element_data[NIII][TWOLEVELFIELD_ENERGY_DIFFERENCE] =
        0.00159 * Ry_over_k;
    // Galavis, Mendoza & Zeippen (1998), table 4, 1 to 2 transition (in s^-1)
    _two_level_element_data[NIII][TWOLEVELFIELD_TRANSITION_PROBABILITY] =
        4.736e-5;
    // Blum & Pradhan (1992), table 3, value for 10,000 K, 1 to 2 transition
    _two_level_element_data[NIII][TWOLEVELFIELD_COLLISION_STRENGTH] = 1.4454;
    // statistical weights: level 0 is a P_{1/2} level, while level 1 is a
    // P_{3/2}
    _two_level_element_data[NIII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0] =
        0.5;
    _two_level_element_data[NIII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1] =
        0.25;
  }

  /// NeII
  {
    // Saraph & Tully (1994), table 2, fine structure splitting energy for
    // Z = 10 (in Ry)
    // ground state: 2P3/2
    // excited state: 2P1/2
    _two_level_element_data[NeII][TWOLEVELFIELD_ENERGY_DIFFERENCE] =
        0.0071 * Ry_over_k;
    // Kaufman & Sugar (1986), table 7 (in s^-1)
    _two_level_element_data[NeII][TWOLEVELFIELD_TRANSITION_PROBABILITY] =
        8.55e-3;
    // Griffin, Mitnik & Badnell (2001), table 4, value for 10,000 K
    _two_level_element_data[NeII][TWOLEVELFIELD_COLLISION_STRENGTH] = 0.314;
    // statistical weights: level 0 is a P_{3/2} level, while level 1 is a
    // P_{1/2}
    _two_level_element_data[NeII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0] =
        0.25;
    _two_level_element_data[NeII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1] =
        0.5;
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
 * @brief Get the velocity-averaged collision strength for the given transition
 * of the given element (at 10,000 K).
 *
 * @param element LineCoolingDataFiveLevelElement.
 * @param transition LineCoolingDataTransition.
 * @return Velocity-averaged collision strength (at 10,000 K).
 */
double LineCoolingData::get_collision_strength(
    LineCoolingDataFiveLevelElement element,
    LineCoolingDataTransition transition) const {
  return _collision_strength[element][transition];
}

/**
 * @brief Get the exponent for the temperature variation of the collision
 * strength for the given transition of the given element.
 *
 * @param element LineCoolingDataFiveLevelElement.
 * @param transition LineCoolingDataTransition.
 * @return Exponent for the temperature variation of the collision strength.
 */
double LineCoolingData::get_collision_strength_exponent(
    LineCoolingDataFiveLevelElement element,
    LineCoolingDataTransition transition) const {
  return _collision_strength_exponent[element][transition];
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
 * @param Tinv Inverse of the temperature (in K^-1).
 * @param T4 Temperature (in 10^4 K).
 * @param level_populations Array to store the resulting level populations in.
 */
void LineCoolingData::compute_level_populations(
    LineCoolingDataFiveLevelElement element,
    double collision_strength_prefactor, double Tinv, double T4,
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
        collision_strength_prefactor * _collision_strength[element][i] *
        std::pow(T4, _collision_strength_exponent[element][i]);
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

  const double ksi =
      _two_level_element_data[element][TWOLEVELFIELD_ENERGY_DIFFERENCE];
  const double A =
      _two_level_element_data[element][TWOLEVELFIELD_TRANSITION_PROBABILITY];
  const double Gamma =
      _two_level_element_data[element][TWOLEVELFIELD_COLLISION_STRENGTH];
  const double inv_omega_1 =
      _two_level_element_data[element]
                             [TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0];
  const double inv_omega_2 =
      _two_level_element_data[element]
                             [TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1];
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
  const double T4 = temperature * 1.e-4;
  const double Tinv = 1. / temperature;

  /// five level elements

  double cooling = 0.;
  for (int j = 0; j < LINECOOLINGDATA_NUMFIVELEVELELEMENTS; ++j) {

    const LineCoolingDataFiveLevelElement element =
        static_cast< LineCoolingDataFiveLevelElement >(j);

    double level_populations[5];
    compute_level_populations(element, collision_strength_prefactor, Tinv, T4,
                              level_populations);

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
  for (int i = 0; i < LINECOOLINGDATA_NUMTWOLEVELELEMENTS; ++i) {
    const LineCoolingDataTwoLevelElement element =
        static_cast< LineCoolingDataTwoLevelElement >(i);
    const double level_population =
        compute_level_population(element, collision_strength_prefactor, Tinv);
    cooling +=
        abundances[offset + i] * kb *
        _two_level_element_data[i][TWOLEVELFIELD_ENERGY_DIFFERENCE] *
        _two_level_element_data[i][TWOLEVELFIELD_TRANSITION_PROBABILITY] *
        level_population;
  }

  return cooling;
}

/**
 * @brief Calculate the strength of a number of emission lines for the given
 * temperature, electron density and ion abundances.
 *
 * The matching of lines to transitions is based on tables 3.8-14 in Osterbrock
 * & Ferland (2006) and similar data in the data papers used for transition
 * data.
 *
 * @param temperature Temperature (in K).
 * @param electron_density Electron density (in m^-3).
 * @param abundances Ion abundances.
 * @param c6300_6363 Variable to store the oxygen 6300 and 6363 angstrom
 * emission line strengths in (in J s^-1).
 * @param c9405 Variable to store the sulphur 9405 angstrom emission line
 * strength in (in J s^-1).
 * @param c6312 Variable to store the sulphur 6312 angstrom emission line
 * strength in (in J s^-1).
 * @param c33mu Variable that is not used in the output (in J s^-1).
 * @param c19mu Variable to store the sulphur 18.7 micrometre emission line
 * strength in (in J s^-1).
 * @param c3729 Variable to store the oxygen 3729 angstrom emission line
 * strength in (in J s^-1).
 * @param c3727 Variable to store the oxygen 3727 angstrom emission line
 * strength in (in J s^-1).
 * @param c7330 Variable to store the oxygen 7330 angstrom emission line
 * strength in (in J s^-1).
 * @param c4363 Variable to store the oxygen 4363 angstrom emission line
 * strength in (in J s^-1).
 * @param c5007 Variable to store the oxygen 5007 angstrom emission line
 * strength in (in J s^-1).
 * @param c52mu Variable to store the oxygen 52 micrometre emission line
 * strength in (in J s^-1).
 * @param c88mu Variable to store the oxygen 88 micrometre emission line
 * strength in (in J s^-1).
 * @param c5755 Variable to store the nitrogen 5755 angstrom emission line
 * strength in (in J s^-1).
 * @param c6584 Variable to store the nitrogen 6584 angstrom emission line
 * strength in (in J s^-1).
 * @param c4072 Variable to store the sulphur 4072 angstrom emission line
 * strength in (in J s^-1).
 * @param c6717 Variable to store the sulphur 6717 angstrom emission line
 * strength in (in J s^-1).
 * @param c6725 Variable to store the sulphur 6725 angstrom emission line
 * strength in (in J s^-1).
 * @param c3869 Variable to store the neon 3869 angstrom emission line strength
 * in (in J s^-1).
 * @param cniii57 Variable to store the nitrogen 57.3 micrometre emission line
 * strength in (in J s^-1).
 * @param cneii12 Variable to store the neon 12.8 micrometre emission line
 * strength in (in J s^-1).
 * @param cneiii15 Variable to store the neon 15.5 micrometre emission line
 * strength in (in J s^-1).
 * @param cnii122 Variable to store the nitrogen 122 micrometre emission line
 * strength in (in J s^-1).
 * @param cii2325 Variable to store the carbon 2325 angstrom emission line
 * strength in (in J s^-1).
 * @param ciii1908 Variable to store the carbon 1907 + 1909 angstrom emission
 * line strength in (in J s^-1).
 * @param coii7325 Variable to store the oxygen 7320 + 7330 angstrom emission
 * line strength in (in J s^-1).
 * @param csiv10 Variable to store the sulphur 10 micrometre (?) emission line
 * strength in (in J s^-1).
 */
void LineCoolingData::linestr(
    double temperature, double electron_density, const double *abundances,
    double &c6300_6363, double &c9405, double &c6312, double &c33mu,
    double &c19mu, double &c3729, double &c3727, double &c7330, double &c4363,
    double &c5007, double &c52mu, double &c88mu, double &c5755, double &c6584,
    double &c4072, double &c6717, double &c6725, double &c3869, double &cniii57,
    double &cneii12, double &cneiii15, double &cnii122, double &cii2325,
    double &ciii1908, double &coii7325, double &csiv10) const {

  /// initialize some variables

  // Boltzmann constant (in J s^-1)
  const double kb =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

  const double collision_strength_prefactor =
      _collision_strength_prefactor * electron_density / std::sqrt(temperature);
  const double T4 = temperature * 1.e-4;
  const double Tinv = 1. / temperature;

  /// 5 level elements

  // there are no lines for element 0 (NI), so we skip the iteration for that
  // element
  for (int j = 1; j < 10; ++j) {

    const LineCoolingDataFiveLevelElement element =
        static_cast< LineCoolingDataFiveLevelElement >(j);

    double level_populations[5];
    compute_level_populations(element, collision_strength_prefactor, Tinv, T4,
                              level_populations);

    const double prefactor = abundances[j] * kb;

    if (element == NII) {

      // Osterbrock & Ferland (2006), table 3.12
      // ground state: 3P0
      // excited states: 3P1, 3P2, 1D2, 1S0
      c5755 = prefactor * level_populations[4] *
              _transition_probability[NII][TRANSITION_3_to_4] *
              _energy_difference[NII][TRANSITION_3_to_4];
      c6584 = prefactor * level_populations[3] *
              _transition_probability[NII][TRANSITION_2_to_3] *
              _energy_difference[NII][TRANSITION_2_to_3];
      cnii122 = prefactor * level_populations[2] *
                _transition_probability[NII][TRANSITION_1_to_2] *
                _energy_difference[j][TRANSITION_1_to_2];

    } else if (element == OI) {

      // Osterbrock & Ferland (2006), table 3.14
      // ground state: 3P2
      // excited states: 3P1, 3P0, 1D2, 1S0
      // this is the sum of the 6300.3 and 6363.8 angstrom transitions
      c6300_6363 = prefactor * level_populations[3] *
                   (_transition_probability[OI][TRANSITION_0_to_3] *
                        _energy_difference[OI][TRANSITION_0_to_3] +
                    _transition_probability[OI][TRANSITION_1_to_3] *
                        _energy_difference[OI][TRANSITION_1_to_3]);

    } else if (element == OII) {

      // Osterbrock & Ferland (2006), table 3.13
      // ground state: 4S3/2
      // excited states: 2D5/2, 2D3/2, 2P3/2, 2P1/2

      c3729 = prefactor * level_populations[1] *
              _transition_probability[OII][TRANSITION_0_to_1] *
              _energy_difference[OII][TRANSITION_0_to_1];
      // this is the sum of the 3726.0 and 3728.8 angstrom transitions
      // note that Kenny's version wrongly included the 497.1 um transition as
      // well...
      c3727 = prefactor * (level_populations[1] *
                               _transition_probability[OII][TRANSITION_0_to_1] *
                               _energy_difference[OII][TRANSITION_0_to_1] +
                           level_populations[2] *
                               _transition_probability[OII][TRANSITION_0_to_2] *
                               _energy_difference[OII][TRANSITION_0_to_2]);
      // this is the sum of the transitions at 7319.9, 7330.7, 7318.8 and 7329.6
      // angstrom
      coii7325 =
          prefactor * (level_populations[4] *
                           (_transition_probability[OII][TRANSITION_1_to_4] *
                                _energy_difference[OII][TRANSITION_1_to_4] +
                            _transition_probability[OII][TRANSITION_2_to_4] *
                                _energy_difference[OII][TRANSITION_2_to_4]) +
                       level_populations[3] *
                           (_transition_probability[OII][TRANSITION_1_to_3] *
                                _energy_difference[OII][TRANSITION_1_to_3] +
                            _transition_probability[OII][TRANSITION_2_to_3] *
                                _energy_difference[OII][TRANSITION_2_to_3]));

    } else if (element == OIII) {

      // Osterbrock & Ferland (2006), table 3.12
      // ground state: 3P0
      // excited states: 3P1, 3P2, 1D2, 1S0
      c4363 = prefactor * level_populations[4] *
              _transition_probability[OIII][TRANSITION_3_to_4] *
              _energy_difference[OIII][TRANSITION_3_to_4];
      c5007 = prefactor * level_populations[3] *
              _transition_probability[OIII][TRANSITION_2_to_3] *
              _energy_difference[OIII][TRANSITION_2_to_3];
      c52mu = prefactor * level_populations[2] *
              _transition_probability[OIII][TRANSITION_1_to_2] *
              _energy_difference[OIII][TRANSITION_1_to_2];
      c88mu = prefactor * level_populations[1] *
              _transition_probability[OIII][TRANSITION_0_to_1] *
              _energy_difference[OIII][TRANSITION_0_to_1];

    } else if (element == NeIII) {

      // Osterbrock & Ferland (2006), table 3.14
      // ground state: 3P2
      // excited states: 3P1, 3P0, 1D2, 1S0
      c3869 = prefactor * level_populations[3] *
              _transition_probability[NeIII][TRANSITION_0_to_3] *
              _energy_difference[NeIII][TRANSITION_0_to_3];
      cneiii15 = prefactor * level_populations[1] *
                 _transition_probability[NeIII][TRANSITION_0_to_1] *
                 _energy_difference[NeIII][TRANSITION_0_to_1];

    } else if (element == SII) {

      // Osterbrock & Ferland (2006), table 3.13
      // ground state: 4S3/2
      // excited states: 2D3/2, 2D5/2, 2P1/2, 2P3/2
      // this is the sum of the 4068.6 and 4076.4 angstrom transitions
      c4072 = prefactor * (level_populations[3] *
                               _transition_probability[SII][TRANSITION_0_to_3] *
                               _energy_difference[SII][TRANSITION_0_to_3] +
                           level_populations[4] *
                               _transition_probability[SII][TRANSITION_0_to_4] *
                               _energy_difference[SII][TRANSITION_0_to_4]);
      // note that Kenny's version wrongly includes the 314.5 um transition...
      c6717 = prefactor * level_populations[2] *
              _transition_probability[SII][TRANSITION_0_to_2] *
              _energy_difference[SII][TRANSITION_0_to_2];
      // this is the sum of the 6716.5 and 6730.8 angstrom transitions
      // note that Kenny's version wrongly includes the 314.5 um transition...
      c6725 = prefactor * (level_populations[1] *
                               _transition_probability[SII][TRANSITION_0_to_1] *
                               _energy_difference[SII][TRANSITION_0_to_1] +
                           level_populations[2] *
                               _transition_probability[SII][TRANSITION_0_to_2] *
                               _energy_difference[SII][TRANSITION_0_to_2]);

    } else if (element == SIII) {

      // Osterbrock & Ferland (2006), table 3.12
      // ground state: 3P0
      // excited states: 3P1, 3P2, 1D2, 1S0
      // this is the sum of the 9531.0 and 9068.9 angstrom transitions
      c9405 = prefactor * level_populations[3] *
              (_transition_probability[SIII][TRANSITION_1_to_3] *
                   _energy_difference[SIII][TRANSITION_1_to_3] +
               _transition_probability[SIII][TRANSITION_2_to_3] *
                   _energy_difference[SIII][TRANSITION_2_to_3]);
      c6312 = prefactor * level_populations[4] *
              _transition_probability[SIII][TRANSITION_3_to_4] *
              _energy_difference[SIII][TRANSITION_3_to_4];
      c33mu = prefactor * level_populations[1] *
              _transition_probability[SIII][TRANSITION_0_to_1] *
              _energy_difference[SIII][TRANSITION_0_to_1];
      c19mu = prefactor * level_populations[2] *
              _transition_probability[SIII][TRANSITION_1_to_2] *
              _energy_difference[SIII][TRANSITION_1_to_2];

    } else if (element == CII) {

      // Osterbrock & Ferland (2006), table 3.9
      // ground state: 2P1/2
      // excited states: 2P3/2, 4P1/2, 4P3/2, 4P5/2
      // this should be the sum of all 4P to 2P transitions
      // note that Kenny's code wrongly includes some 4P to 4P transitions...
      cii2325 =
          prefactor * (level_populations[2] *
                           (_transition_probability[CII][TRANSITION_0_to_2] *
                                _energy_difference[CII][TRANSITION_0_to_2] +
                            _transition_probability[CII][TRANSITION_1_to_2] *
                                _energy_difference[CII][TRANSITION_1_to_2]) +
                       level_populations[3] *
                           (_transition_probability[CII][TRANSITION_0_to_3] *
                                _energy_difference[CII][TRANSITION_0_to_3] +
                            _transition_probability[CII][TRANSITION_1_to_3] *
                                _energy_difference[CII][TRANSITION_1_to_3]) +
                       level_populations[4] *
                           (_transition_probability[CII][TRANSITION_0_to_4] *
                                _energy_difference[CII][TRANSITION_0_to_4] +
                            _transition_probability[CII][TRANSITION_1_to_4] *
                                _energy_difference[CII][TRANSITION_1_to_4]));

    } else if (element == CIII) {

      // Osterbrock & Ferland (2006), table 3.8
      // ground state: 1S0
      // excited states: 3P0, 3P1, 3P2, 1P1
      // this is the sum of all 3P to 1S transitions
      // note that Kenny's code wrongly includes some 3P to 3P transitions...
      ciii1908 =
          prefactor * (level_populations[1] *
                           _transition_probability[CIII][TRANSITION_0_to_1] *
                           _energy_difference[CIII][TRANSITION_0_to_1] +
                       level_populations[2] *
                           _transition_probability[CIII][TRANSITION_0_to_2] *
                           _energy_difference[CIII][TRANSITION_0_to_2] +
                       level_populations[3] *
                           _transition_probability[CIII][TRANSITION_0_to_3] *
                           _energy_difference[CIII][TRANSITION_0_to_3]);
    }
  }

  /// 2 level elements

  // offset of two level elements in the abundances array
  const int offset = LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  for (int i = 0; i < LINECOOLINGDATA_NUMTWOLEVELELEMENTS; ++i) {

    const LineCoolingDataTwoLevelElement element =
        static_cast< LineCoolingDataTwoLevelElement >(i);

    const double level_population =
        compute_level_population(element, collision_strength_prefactor, Tinv);

    if (element == NIII) {

      // Osterbrock & Ferland (2006), table 3.9
      // ground state: 2P1/2
      // excited state: 2P3/2
      cniii57 =
          abundances[offset + NIII] * kb *
          _two_level_element_data[NIII][TWOLEVELFIELD_ENERGY_DIFFERENCE] *
          _two_level_element_data[NIII][TWOLEVELFIELD_TRANSITION_PROBABILITY] *
          level_population;

    } else if (element == NeII) {

      // Osterbrock & Ferland (2006), table 3.11
      // ground state: 2P3/2
      // excited state: 2P1/2
      cneii12 =
          abundances[offset + NeII] * kb *
          _two_level_element_data[NeII][TWOLEVELFIELD_ENERGY_DIFFERENCE] *
          _two_level_element_data[NeII][TWOLEVELFIELD_TRANSITION_PROBABILITY] *
          level_population;
    }
  }
}
