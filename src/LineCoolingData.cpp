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
    _five_level_inverse_statistical_weight[NI][0] = 0.25;
    _five_level_inverse_statistical_weight[NI][1] = 1. / 6.;
    _five_level_inverse_statistical_weight[NI][2] = 0.25;
    _five_level_inverse_statistical_weight[NI][3] = 0.5;
    _five_level_inverse_statistical_weight[NI][4] = 0.25;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[NI][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[NI][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we take the sum of all transition probabilities (in s^-1)
    _five_level_transition_probability[NI][TRANSITION_0_to_1] = 7.566e-6;
    _five_level_transition_probability[NI][TRANSITION_0_to_2] = 2.029e-5;
    _five_level_transition_probability[NI][TRANSITION_0_to_3] = 2.605e-3;
    _five_level_transition_probability[NI][TRANSITION_0_to_4] = 6.504e-3;
    _five_level_transition_probability[NI][TRANSITION_1_to_2] = 1.071e-8;
    _five_level_transition_probability[NI][TRANSITION_1_to_3] = 3.447e-2;
    _five_level_transition_probability[NI][TRANSITION_1_to_4] = 6.135e-2;
    _five_level_transition_probability[NI][TRANSITION_2_to_3] = 5.272e-2;
    _five_level_transition_probability[NI][TRANSITION_2_to_4] = 2.745e-2;
    _five_level_transition_probability[NI][TRANSITION_3_to_4] = 9.504e-17;
    // our own fits to the data of Tayal (2000)
    // these fits were made with the script data/linecooling/gamma_NI.py
    // (this script also outputs the code below)
    _five_level_collision_strength[NI][TRANSITION_0_to_1][0] = 4.350e-04;
    _five_level_collision_strength[NI][TRANSITION_0_to_1][1] = -4.500e-04;
    _five_level_collision_strength[NI][TRANSITION_0_to_1][2] = 2.206e-01;
    _five_level_collision_strength[NI][TRANSITION_0_to_1][3] = 4.958e-05;
    _five_level_collision_strength[NI][TRANSITION_0_to_1][4] = 7.786e-02;
    _five_level_collision_strength[NI][TRANSITION_0_to_1][5] = -1.099e-07;
    _five_level_collision_strength[NI][TRANSITION_0_to_1][6] = -8.271e-09;
    _five_level_collision_strength[NI][TRANSITION_0_to_2][0] = 2.751e-06;
    _five_level_collision_strength[NI][TRANSITION_0_to_2][1] = -3.007e-04;
    _five_level_collision_strength[NI][TRANSITION_0_to_2][2] = 1.474e-01;
    _five_level_collision_strength[NI][TRANSITION_0_to_2][3] = 3.312e-05;
    _five_level_collision_strength[NI][TRANSITION_0_to_2][4] = 3.882e-01;
    _five_level_collision_strength[NI][TRANSITION_0_to_2][5] = -1.470e-08;
    _five_level_collision_strength[NI][TRANSITION_0_to_2][6] = -1.106e-09;
    _five_level_collision_strength[NI][TRANSITION_0_to_3][0] = 2.254e-07;
    _five_level_collision_strength[NI][TRANSITION_0_to_3][1] = -3.336e-05;
    _five_level_collision_strength[NI][TRANSITION_0_to_3][2] = 2.277e-02;
    _five_level_collision_strength[NI][TRANSITION_0_to_3][3] = 3.636e-06;
    _five_level_collision_strength[NI][TRANSITION_0_to_3][4] = 8.211e-02;
    _five_level_collision_strength[NI][TRANSITION_0_to_3][5] = -5.171e-09;
    _five_level_collision_strength[NI][TRANSITION_0_to_3][6] = -3.695e-10;
    _five_level_collision_strength[NI][TRANSITION_0_to_4][0] = 6.567e-06;
    _five_level_collision_strength[NI][TRANSITION_0_to_4][1] = -6.814e-05;
    _five_level_collision_strength[NI][TRANSITION_0_to_4][2] = 4.632e-02;
    _five_level_collision_strength[NI][TRANSITION_0_to_4][3] = 7.429e-06;
    _five_level_collision_strength[NI][TRANSITION_0_to_4][4] = 5.470e-02;
    _five_level_collision_strength[NI][TRANSITION_0_to_4][5] = -1.612e-08;
    _five_level_collision_strength[NI][TRANSITION_0_to_4][6] = -1.155e-09;
    _five_level_collision_strength[NI][TRANSITION_1_to_2][0] = 7.907e-02;
    _five_level_collision_strength[NI][TRANSITION_1_to_2][1] = -8.603e-03;
    _five_level_collision_strength[NI][TRANSITION_1_to_2][2] = 6.355e+00;
    _five_level_collision_strength[NI][TRANSITION_1_to_2][3] = 9.366e-04;
    _five_level_collision_strength[NI][TRANSITION_1_to_2][4] = 1.259e-01;
    _five_level_collision_strength[NI][TRANSITION_1_to_2][5] = -1.315e-06;
    _five_level_collision_strength[NI][TRANSITION_1_to_2][6] = -9.950e-08;
    _five_level_collision_strength[NI][TRANSITION_1_to_3][0] = -1.420e-02;
    _five_level_collision_strength[NI][TRANSITION_1_to_3][1] = -8.743e-03;
    _five_level_collision_strength[NI][TRANSITION_1_to_3][2] = 3.581e+00;
    _five_level_collision_strength[NI][TRANSITION_1_to_3][3] = 9.893e-04;
    _five_level_collision_strength[NI][TRANSITION_1_to_3][4] = 8.734e-01;
    _five_level_collision_strength[NI][TRANSITION_1_to_3][5] = -2.580e-07;
    _five_level_collision_strength[NI][TRANSITION_1_to_3][6] = -1.983e-08;
    _five_level_collision_strength[NI][TRANSITION_1_to_4][0] = -3.258e-01;
    _five_level_collision_strength[NI][TRANSITION_1_to_4][1] = -1.290e+00;
    _five_level_collision_strength[NI][TRANSITION_1_to_4][2] = 5.455e+02;
    _five_level_collision_strength[NI][TRANSITION_1_to_4][3] = 1.451e-01;
    _five_level_collision_strength[NI][TRANSITION_1_to_4][4] = -1.122e+01;
    _five_level_collision_strength[NI][TRANSITION_1_to_4][5] = 2.792e-06;
    _five_level_collision_strength[NI][TRANSITION_1_to_4][6] = 2.138e-07;
    _five_level_collision_strength[NI][TRANSITION_2_to_3][0] = -3.592e-01;
    _five_level_collision_strength[NI][TRANSITION_2_to_3][1] = -9.058e-01;
    _five_level_collision_strength[NI][TRANSITION_2_to_3][2] = 3.817e+02;
    _five_level_collision_strength[NI][TRANSITION_2_to_3][3] = 1.019e-01;
    _five_level_collision_strength[NI][TRANSITION_2_to_3][4] = -1.055e+01;
    _five_level_collision_strength[NI][TRANSITION_2_to_3][5] = 2.096e-06;
    _five_level_collision_strength[NI][TRANSITION_2_to_3][6] = 1.606e-07;
    _five_level_collision_strength[NI][TRANSITION_2_to_4][0] = -1.910e-01;
    _five_level_collision_strength[NI][TRANSITION_2_to_4][1] = -9.570e-02;
    _five_level_collision_strength[NI][TRANSITION_2_to_4][2] = 4.070e+01;
    _five_level_collision_strength[NI][TRANSITION_2_to_4][3] = 1.077e-02;
    _five_level_collision_strength[NI][TRANSITION_2_to_4][4] = -5.010e+00;
    _five_level_collision_strength[NI][TRANSITION_2_to_4][5] = 4.665e-07;
    _five_level_collision_strength[NI][TRANSITION_2_to_4][6] = 3.574e-08;
    _five_level_collision_strength[NI][TRANSITION_3_to_4][0] = -2.313e-03;
    _five_level_collision_strength[NI][TRANSITION_3_to_4][1] = -2.099e-02;
    _five_level_collision_strength[NI][TRANSITION_3_to_4][2] = 1.122e+01;
    _five_level_collision_strength[NI][TRANSITION_3_to_4][3] = 2.337e-03;
    _five_level_collision_strength[NI][TRANSITION_3_to_4][4] = 2.240e+00;
    _five_level_collision_strength[NI][TRANSITION_3_to_4][5] = -2.197e-07;
    _five_level_collision_strength[NI][TRANSITION_3_to_4][6] = -1.682e-08;
  }

  /// NII
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 4 and 5
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {48.7, 130.8, 15316.3, 32688.9};
    _five_level_inverse_statistical_weight[NII][0] = 1.;
    _five_level_inverse_statistical_weight[NII][1] = 1. / 3.;
    _five_level_inverse_statistical_weight[NII][2] = 0.2;
    _five_level_inverse_statistical_weight[NII][3] = 0.2;
    _five_level_inverse_statistical_weight[NII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[NII][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[NII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _five_level_transition_probability[NII][TRANSITION_0_to_1] = 2.077e-6;
    _five_level_transition_probability[NII][TRANSITION_0_to_2] = 1.127e-12;
    _five_level_transition_probability[NII][TRANSITION_0_to_3] = 3.554e-7;
    _five_level_transition_probability[NII][TRANSITION_0_to_4] = 0.;
    _five_level_transition_probability[NII][TRANSITION_1_to_2] = 7.463e-6;
    _five_level_transition_probability[NII][TRANSITION_1_to_3] = 1.016e-3;
    _five_level_transition_probability[NII][TRANSITION_1_to_4] = 3.297e-2;
    _five_level_transition_probability[NII][TRANSITION_2_to_3] = 3.005e-3;
    _five_level_transition_probability[NII][TRANSITION_2_to_4] = 1.315e-4;
    _five_level_transition_probability[NII][TRANSITION_3_to_4] = 1.023;
    // our own fits to the data of Lennon & Burke (1994)
    // these fits were made with the script data/linecooling/gamma_NII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[NII][TRANSITION_0_to_1][0] = -3.483e-01;
    _five_level_collision_strength[NII][TRANSITION_0_to_1][1] = 4.384e-03;
    _five_level_collision_strength[NII][TRANSITION_0_to_1][2] = 2.059e+00;
    _five_level_collision_strength[NII][TRANSITION_0_to_1][3] = -4.032e-04;
    _five_level_collision_strength[NII][TRANSITION_0_to_1][4] = -3.391e-05;
    _five_level_collision_strength[NII][TRANSITION_0_to_1][5] = -1.192e-03;
    _five_level_collision_strength[NII][TRANSITION_0_to_1][6] = -8.966e-05;
    _five_level_collision_strength[NII][TRANSITION_0_to_2][0] = -2.779e-03;
    _five_level_collision_strength[NII][TRANSITION_0_to_2][1] = -2.640e-06;
    _five_level_collision_strength[NII][TRANSITION_0_to_2][2] = 2.127e-01;
    _five_level_collision_strength[NII][TRANSITION_0_to_2][3] = 1.368e-06;
    _five_level_collision_strength[NII][TRANSITION_0_to_2][4] = -7.504e-04;
    _five_level_collision_strength[NII][TRANSITION_0_to_2][5] = 1.570e-06;
    _five_level_collision_strength[NII][TRANSITION_0_to_2][6] = 1.229e-07;
    _five_level_collision_strength[NII][TRANSITION_0_to_3][0] = -1.472e-02;
    _five_level_collision_strength[NII][TRANSITION_0_to_3][1] = 2.665e-05;
    _five_level_collision_strength[NII][TRANSITION_0_to_3][2] = 2.965e-01;
    _five_level_collision_strength[NII][TRANSITION_0_to_3][3] = -2.560e-06;
    _five_level_collision_strength[NII][TRANSITION_0_to_3][4] = 1.796e-07;
    _five_level_collision_strength[NII][TRANSITION_0_to_3][5] = 1.600e-03;
    _five_level_collision_strength[NII][TRANSITION_0_to_3][6] = 1.209e-04;
    _five_level_collision_strength[NII][TRANSITION_0_to_4][0] = 2.426e-03;
    _five_level_collision_strength[NII][TRANSITION_0_to_4][1] = -1.170e-07;
    _five_level_collision_strength[NII][TRANSITION_0_to_4][2] = 3.066e-02;
    _five_level_collision_strength[NII][TRANSITION_0_to_4][3] = 2.803e-08;
    _five_level_collision_strength[NII][TRANSITION_0_to_4][4] = 2.259e-09;
    _five_level_collision_strength[NII][TRANSITION_0_to_4][5] = -2.924e-03;
    _five_level_collision_strength[NII][TRANSITION_0_to_4][6] = -2.031e-04;
    _five_level_collision_strength[NII][TRANSITION_1_to_2][0] = -1.015e-01;
    _five_level_collision_strength[NII][TRANSITION_1_to_2][1] = 6.264e-04;
    _five_level_collision_strength[NII][TRANSITION_1_to_2][2] = 1.559e+00;
    _five_level_collision_strength[NII][TRANSITION_1_to_2][3] = -5.445e-05;
    _five_level_collision_strength[NII][TRANSITION_1_to_2][4] = 1.428e-06;
    _five_level_collision_strength[NII][TRANSITION_1_to_2][5] = 5.929e-04;
    _five_level_collision_strength[NII][TRANSITION_1_to_2][6] = 3.183e-05;
    _five_level_collision_strength[NII][TRANSITION_1_to_3][0] = -1.472e-02;
    _five_level_collision_strength[NII][TRANSITION_1_to_3][1] = 7.995e-05;
    _five_level_collision_strength[NII][TRANSITION_1_to_3][2] = 8.894e-01;
    _five_level_collision_strength[NII][TRANSITION_1_to_3][3] = -7.680e-06;
    _five_level_collision_strength[NII][TRANSITION_1_to_3][4] = 5.191e-07;
    _five_level_collision_strength[NII][TRANSITION_1_to_3][5] = 1.661e-03;
    _five_level_collision_strength[NII][TRANSITION_1_to_3][6] = 1.255e-04;
    _five_level_collision_strength[NII][TRANSITION_1_to_4][0] = 2.426e-03;
    _five_level_collision_strength[NII][TRANSITION_1_to_4][1] = -3.510e-07;
    _five_level_collision_strength[NII][TRANSITION_1_to_4][2] = 9.197e-02;
    _five_level_collision_strength[NII][TRANSITION_1_to_4][3] = 8.408e-08;
    _five_level_collision_strength[NII][TRANSITION_1_to_4][4] = 7.197e-09;
    _five_level_collision_strength[NII][TRANSITION_1_to_4][5] = -2.753e-03;
    _five_level_collision_strength[NII][TRANSITION_1_to_4][6] = -1.912e-04;
    _five_level_collision_strength[NII][TRANSITION_2_to_3][0] = -1.473e-02;
    _five_level_collision_strength[NII][TRANSITION_2_to_3][1] = 1.333e-04;
    _five_level_collision_strength[NII][TRANSITION_2_to_3][2] = 1.482e+00;
    _five_level_collision_strength[NII][TRANSITION_2_to_3][3] = -1.280e-05;
    _five_level_collision_strength[NII][TRANSITION_2_to_3][4] = -1.036e-06;
    _five_level_collision_strength[NII][TRANSITION_2_to_3][5] = -1.390e-03;
    _five_level_collision_strength[NII][TRANSITION_2_to_3][6] = -1.048e-04;
    _five_level_collision_strength[NII][TRANSITION_2_to_4][0] = 2.425e-03;
    _five_level_collision_strength[NII][TRANSITION_2_to_4][1] = -5.846e-07;
    _five_level_collision_strength[NII][TRANSITION_2_to_4][2] = 1.533e-01;
    _five_level_collision_strength[NII][TRANSITION_2_to_4][3] = 1.401e-07;
    _five_level_collision_strength[NII][TRANSITION_2_to_4][4] = -1.952e-07;
    _five_level_collision_strength[NII][TRANSITION_2_to_4][5] = 1.690e-04;
    _five_level_collision_strength[NII][TRANSITION_2_to_4][6] = 1.175e-05;
    _five_level_collision_strength[NII][TRANSITION_3_to_4][0] = 2.226e-02;
    _five_level_collision_strength[NII][TRANSITION_3_to_4][1] = -3.918e-04;
    _five_level_collision_strength[NII][TRANSITION_3_to_4][2] = 1.121e+00;
    _five_level_collision_strength[NII][TRANSITION_3_to_4][3] = 3.934e-05;
    _five_level_collision_strength[NII][TRANSITION_3_to_4][4] = 6.245e-06;
    _five_level_collision_strength[NII][TRANSITION_3_to_4][5] = -7.473e-04;
    _five_level_collision_strength[NII][TRANSITION_3_to_4][6] = -5.553e-05;
  }

  /// OI
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 6 and 7
    // ground state: 3P2
    // excited states: 3P1, 3P0, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {158., 227., 15868., 33793.};
    _five_level_inverse_statistical_weight[OI][0] = 0.2;
    _five_level_inverse_statistical_weight[OI][1] = 1. / 3.;
    _five_level_inverse_statistical_weight[OI][2] = 1.;
    _five_level_inverse_statistical_weight[OI][3] = 0.2;
    _five_level_inverse_statistical_weight[OI][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[OI][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[OI][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _five_level_transition_probability[OI][TRANSITION_0_to_1] = 8.865e-5;
    _five_level_transition_probability[OI][TRANSITION_0_to_2] = 1.275e-10;
    _five_level_transition_probability[OI][TRANSITION_0_to_3] = 6.535e-3;
    _five_level_transition_probability[OI][TRANSITION_0_to_4] = 2.945e-4;
    _five_level_transition_probability[OI][TRANSITION_1_to_2] = 1.772e-5;
    _five_level_transition_probability[OI][TRANSITION_1_to_3] = 2.111e-3;
    _five_level_transition_probability[OI][TRANSITION_1_to_4] = 7.909e-2;
    _five_level_transition_probability[OI][TRANSITION_2_to_3] = 6.388e-7;
    _five_level_transition_probability[OI][TRANSITION_2_to_4] = 0.;
    _five_level_transition_probability[OI][TRANSITION_3_to_4] = 1.124;
    // our own fits to the data of Zatsarinny & Tayal (2003) and Berrington
    // (1988)
    // these fits were made with the script data/linecooling/gamma_OI.py
    // (this script also outputs the code below)
    // there is only data available for fine structure transitions below
    // 10,000 K, so we need to extrapolate for higher temperatures
    _five_level_collision_strength[OI][TRANSITION_0_to_1][0] = -8.839e-01;
    _five_level_collision_strength[OI][TRANSITION_0_to_1][1] = -6.783e-03;
    _five_level_collision_strength[OI][TRANSITION_0_to_1][2] = 7.448e-02;
    _five_level_collision_strength[OI][TRANSITION_0_to_1][3] = 1.578e-03;
    _five_level_collision_strength[OI][TRANSITION_0_to_1][4] = 3.239e-01;
    _five_level_collision_strength[OI][TRANSITION_0_to_1][5] = -2.549e-05;
    _five_level_collision_strength[OI][TRANSITION_0_to_1][6] = -3.727e-06;
    _five_level_collision_strength[OI][TRANSITION_0_to_2][0] = -7.697e-01;
    _five_level_collision_strength[OI][TRANSITION_0_to_2][1] = -9.788e-04;
    _five_level_collision_strength[OI][TRANSITION_0_to_2][2] = 1.200e-02;
    _five_level_collision_strength[OI][TRANSITION_0_to_2][3] = 2.477e-04;
    _five_level_collision_strength[OI][TRANSITION_0_to_2][4] = 9.490e-01;
    _five_level_collision_strength[OI][TRANSITION_0_to_2][5] = -4.462e-07;
    _five_level_collision_strength[OI][TRANSITION_0_to_2][6] = -7.755e-08;
    _five_level_collision_strength[OI][TRANSITION_0_to_3][0] = -1.093e+00;
    _five_level_collision_strength[OI][TRANSITION_0_to_3][1] = 1.413e+00;
    _five_level_collision_strength[OI][TRANSITION_0_to_3][2] = -1.373e+02;
    _five_level_collision_strength[OI][TRANSITION_0_to_3][3] = -1.998e-01;
    _five_level_collision_strength[OI][TRANSITION_0_to_3][4] = 9.365e+00;
    _five_level_collision_strength[OI][TRANSITION_0_to_3][5] = 3.197e-05;
    _five_level_collision_strength[OI][TRANSITION_0_to_3][6] = 2.514e-06;
    _five_level_collision_strength[OI][TRANSITION_0_to_4][0] = -1.219e+00;
    _five_level_collision_strength[OI][TRANSITION_0_to_4][1] = 6.529e-01;
    _five_level_collision_strength[OI][TRANSITION_0_to_4][2] = -5.311e+01;
    _five_level_collision_strength[OI][TRANSITION_0_to_4][3] = -9.393e-02;
    _five_level_collision_strength[OI][TRANSITION_0_to_4][4] = 4.506e+00;
    _five_level_collision_strength[OI][TRANSITION_0_to_4][5] = 2.706e-05;
    _five_level_collision_strength[OI][TRANSITION_0_to_4][6] = 2.088e-06;
    _five_level_collision_strength[OI][TRANSITION_1_to_2][0] = -1.016e+00;
    _five_level_collision_strength[OI][TRANSITION_1_to_2][1] = -6.509e-03;
    _five_level_collision_strength[OI][TRANSITION_1_to_2][2] = 6.884e-02;
    _five_level_collision_strength[OI][TRANSITION_1_to_2][3] = 1.457e-03;
    _five_level_collision_strength[OI][TRANSITION_1_to_2][4] = 1.015e+00;
    _five_level_collision_strength[OI][TRANSITION_1_to_2][5] = -1.947e-05;
    _five_level_collision_strength[OI][TRANSITION_1_to_2][6] = -2.390e-06;
    _five_level_collision_strength[OI][TRANSITION_1_to_3][0] = -1.093e+00;
    _five_level_collision_strength[OI][TRANSITION_1_to_3][1] = 8.480e-01;
    _five_level_collision_strength[OI][TRANSITION_1_to_3][2] = -8.239e+01;
    _five_level_collision_strength[OI][TRANSITION_1_to_3][3] = -1.199e-01;
    _five_level_collision_strength[OI][TRANSITION_1_to_3][4] = 3.496e+00;
    _five_level_collision_strength[OI][TRANSITION_1_to_3][5] = 5.140e-05;
    _five_level_collision_strength[OI][TRANSITION_1_to_3][6] = 4.042e-06;
    _five_level_collision_strength[OI][TRANSITION_1_to_4][0] = -1.219e+00;
    _five_level_collision_strength[OI][TRANSITION_1_to_4][1] = 3.917e-01;
    _five_level_collision_strength[OI][TRANSITION_1_to_4][2] = -3.187e+01;
    _five_level_collision_strength[OI][TRANSITION_1_to_4][3] = -5.636e-02;
    _five_level_collision_strength[OI][TRANSITION_1_to_4][4] = 3.302e+00;
    _five_level_collision_strength[OI][TRANSITION_1_to_4][5] = 2.216e-05;
    _five_level_collision_strength[OI][TRANSITION_1_to_4][6] = 1.709e-06;
    _five_level_collision_strength[OI][TRANSITION_2_to_3][0] = -1.093e+00;
    _five_level_collision_strength[OI][TRANSITION_2_to_3][1] = 2.827e-01;
    _five_level_collision_strength[OI][TRANSITION_2_to_3][2] = -2.746e+01;
    _five_level_collision_strength[OI][TRANSITION_2_to_3][3] = -3.996e-02;
    _five_level_collision_strength[OI][TRANSITION_2_to_3][4] = 2.681e+00;
    _five_level_collision_strength[OI][TRANSITION_2_to_3][5] = 2.234e-05;
    _five_level_collision_strength[OI][TRANSITION_2_to_3][6] = 1.757e-06;
    _five_level_collision_strength[OI][TRANSITION_2_to_4][0] = -1.219e+00;
    _five_level_collision_strength[OI][TRANSITION_2_to_4][1] = 1.306e-01;
    _five_level_collision_strength[OI][TRANSITION_2_to_4][2] = -1.062e+01;
    _five_level_collision_strength[OI][TRANSITION_2_to_4][3] = -1.879e-02;
    _five_level_collision_strength[OI][TRANSITION_2_to_4][4] = 1.964e+00;
    _five_level_collision_strength[OI][TRANSITION_2_to_4][5] = 1.242e-05;
    _five_level_collision_strength[OI][TRANSITION_2_to_4][6] = 9.578e-07;
    _five_level_collision_strength[OI][TRANSITION_3_to_4][0] = -1.056e+00;
    _five_level_collision_strength[OI][TRANSITION_3_to_4][1] = -3.587e-01;
    _five_level_collision_strength[OI][TRANSITION_3_to_4][2] = 3.534e+01;
    _five_level_collision_strength[OI][TRANSITION_3_to_4][3] = 4.778e-02;
    _five_level_collision_strength[OI][TRANSITION_3_to_4][4] = 1.295e+00;
    _five_level_collision_strength[OI][TRANSITION_3_to_4][5] = 1.089e-05;
    _five_level_collision_strength[OI][TRANSITION_3_to_4][6] = 6.545e-07;
  }

  /// OII
  {
    // data from Froese Fischer & Tachiev (2004), table 4
    // ground state: 4S3/2
    // excited states: 2D5/2, 2D3/2, 2P3/2, 2P1/2
    // in cm^-1
    const double energy_levels[4] = {26810.73, 26830.45, 40468.36, 40470.96};
    _five_level_inverse_statistical_weight[OII][0] = 0.25;
    _five_level_inverse_statistical_weight[OII][1] = 1. / 6.;
    _five_level_inverse_statistical_weight[OII][2] = 0.25;
    _five_level_inverse_statistical_weight[OII][3] = 0.25;
    _five_level_inverse_statistical_weight[OII][4] = 0.5;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[OII][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[OII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we sum contributions of all types (in s^-1)
    _five_level_transition_probability[OII][TRANSITION_0_to_1] = 4.124e-5;
    _five_level_transition_probability[OII][TRANSITION_0_to_2] = 1.635e-4;
    _five_level_transition_probability[OII][TRANSITION_0_to_3] = 5.646e-2;
    _five_level_transition_probability[OII][TRANSITION_0_to_4] = 2.265e-2;
    _five_level_transition_probability[OII][TRANSITION_1_to_2] = 1.241e-7;
    _five_level_transition_probability[OII][TRANSITION_1_to_3] = 1.106e-1;
    _five_level_transition_probability[OII][TRANSITION_1_to_4] = 5.824e-2;
    _five_level_transition_probability[OII][TRANSITION_2_to_3] = 5.871e-2;
    _five_level_transition_probability[OII][TRANSITION_2_to_4] = 9.668e-2;
    _five_level_transition_probability[OII][TRANSITION_3_to_4] = 3.158e-10;
    // our own fits to the data of Kisielius et al. (2009)
    // these fits were made with the script data/linecooling/gamma_OII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[OII][TRANSITION_0_to_1][0] = -2.523e-03;
    _five_level_collision_strength[OII][TRANSITION_0_to_1][1] = -6.241e-06;
    _five_level_collision_strength[OII][TRANSITION_0_to_1][2] = 8.502e-01;
    _five_level_collision_strength[OII][TRANSITION_0_to_1][3] = 7.233e-07;
    _five_level_collision_strength[OII][TRANSITION_0_to_1][4] = 1.444e-03;
    _five_level_collision_strength[OII][TRANSITION_0_to_1][5] = 9.709e-09;
    _five_level_collision_strength[OII][TRANSITION_0_to_1][6] = 1.522e-09;
    _five_level_collision_strength[OII][TRANSITION_0_to_2][0] = -1.125e-02;
    _five_level_collision_strength[OII][TRANSITION_0_to_2][1] = -8.446e-07;
    _five_level_collision_strength[OII][TRANSITION_0_to_2][2] = 6.036e-01;
    _five_level_collision_strength[OII][TRANSITION_0_to_2][3] = 2.181e-07;
    _five_level_collision_strength[OII][TRANSITION_0_to_2][4] = -1.635e-04;
    _five_level_collision_strength[OII][TRANSITION_0_to_2][5] = 4.311e-08;
    _five_level_collision_strength[OII][TRANSITION_0_to_2][6] = -2.238e-10;
    _five_level_collision_strength[OII][TRANSITION_0_to_3][0] = 1.161e-02;
    _five_level_collision_strength[OII][TRANSITION_0_to_3][1] = -4.213e-06;
    _five_level_collision_strength[OII][TRANSITION_0_to_3][2] = 2.279e-01;
    _five_level_collision_strength[OII][TRANSITION_0_to_3][3] = 4.878e-07;
    _five_level_collision_strength[OII][TRANSITION_0_to_3][4] = -1.033e-11;
    _five_level_collision_strength[OII][TRANSITION_0_to_3][5] = 1.000e+00;
    _five_level_collision_strength[OII][TRANSITION_0_to_3][6] = 2.265e+00;
    _five_level_collision_strength[OII][TRANSITION_0_to_4][0] = 7.789e-02;
    _five_level_collision_strength[OII][TRANSITION_0_to_4][1] = -1.026e-05;
    _five_level_collision_strength[OII][TRANSITION_0_to_4][2] = 7.464e-02;
    _five_level_collision_strength[OII][TRANSITION_0_to_4][3] = 1.046e-06;
    _five_level_collision_strength[OII][TRANSITION_0_to_4][4] = 1.431e-07;
    _five_level_collision_strength[OII][TRANSITION_0_to_4][5] = -8.913e-04;
    _five_level_collision_strength[OII][TRANSITION_0_to_4][6] = -6.598e-05;
    _five_level_collision_strength[OII][TRANSITION_1_to_2][0] = -8.588e-01;
    _five_level_collision_strength[OII][TRANSITION_1_to_2][1] = -4.648e+00;
    _five_level_collision_strength[OII][TRANSITION_1_to_2][2] = -3.736e+02;
    _five_level_collision_strength[OII][TRANSITION_1_to_2][3] = -1.190e+00;
    _five_level_collision_strength[OII][TRANSITION_1_to_2][4] = -1.466e-06;
    _five_level_collision_strength[OII][TRANSITION_1_to_2][5] = -5.766e+06;
    _five_level_collision_strength[OII][TRANSITION_1_to_2][6] = -9.308e-01;
    _five_level_collision_strength[OII][TRANSITION_1_to_3][0] = -1.028e-02;
    _five_level_collision_strength[OII][TRANSITION_1_to_3][1] = -7.342e-05;
    _five_level_collision_strength[OII][TRANSITION_1_to_3][2] = 8.844e-01;
    _five_level_collision_strength[OII][TRANSITION_1_to_3][3] = 9.292e-06;
    _five_level_collision_strength[OII][TRANSITION_1_to_3][4] = -3.079e-04;
    _five_level_collision_strength[OII][TRANSITION_1_to_3][5] = 7.387e-06;
    _five_level_collision_strength[OII][TRANSITION_1_to_3][6] = 5.594e-07;
    _five_level_collision_strength[OII][TRANSITION_1_to_4][0] = -5.239e-02;
    _five_level_collision_strength[OII][TRANSITION_1_to_4][1] = -6.736e-06;
    _five_level_collision_strength[OII][TRANSITION_1_to_4][2] = 4.666e-01;
    _five_level_collision_strength[OII][TRANSITION_1_to_4][3] = 1.976e-06;
    _five_level_collision_strength[OII][TRANSITION_1_to_4][4] = -8.040e-04;
    _five_level_collision_strength[OII][TRANSITION_1_to_4][5] = 1.289e-06;
    _five_level_collision_strength[OII][TRANSITION_1_to_4][6] = 9.867e-08;
    _five_level_collision_strength[OII][TRANSITION_2_to_3][0] = -3.494e-02;
    _five_level_collision_strength[OII][TRANSITION_2_to_3][1] = -2.313e-05;
    _five_level_collision_strength[OII][TRANSITION_2_to_3][2] = 5.770e-01;
    _five_level_collision_strength[OII][TRANSITION_2_to_3][3] = 3.760e-06;
    _five_level_collision_strength[OII][TRANSITION_2_to_3][4] = -1.210e-03;
    _five_level_collision_strength[OII][TRANSITION_2_to_3][5] = 1.098e-06;
    _five_level_collision_strength[OII][TRANSITION_2_to_3][6] = 8.364e-08;
    _five_level_collision_strength[OII][TRANSITION_2_to_4][0] = 1.038e-02;
    _five_level_collision_strength[OII][TRANSITION_2_to_4][1] = -3.425e-05;
    _five_level_collision_strength[OII][TRANSITION_2_to_4][2] = 2.996e-01;
    _five_level_collision_strength[OII][TRANSITION_2_to_4][3] = 4.003e-06;
    _five_level_collision_strength[OII][TRANSITION_2_to_4][4] = -2.143e-04;
    _five_level_collision_strength[OII][TRANSITION_2_to_4][5] = 3.819e-06;
    _five_level_collision_strength[OII][TRANSITION_2_to_4][6] = 2.881e-07;
    _five_level_collision_strength[OII][TRANSITION_3_to_4][0] = 6.465e-02;
    _five_level_collision_strength[OII][TRANSITION_3_to_4][1] = -3.643e-05;
    _five_level_collision_strength[OII][TRANSITION_3_to_4][2] = 1.838e-01;
    _five_level_collision_strength[OII][TRANSITION_3_to_4][3] = 3.856e-06;
    _five_level_collision_strength[OII][TRANSITION_3_to_4][4] = -1.932e-07;
    _five_level_collision_strength[OII][TRANSITION_3_to_4][5] = 2.875e-03;
    _five_level_collision_strength[OII][TRANSITION_3_to_4][6] = 2.147e-04;
  }

  /// OIII
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 4 and 5
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {114., 307., 20274., 43186.};
    _five_level_inverse_statistical_weight[OIII][0] = 1.;
    _five_level_inverse_statistical_weight[OIII][1] = 1. / 3.;
    _five_level_inverse_statistical_weight[OIII][2] = 0.2;
    _five_level_inverse_statistical_weight[OIII][3] = 0.2;
    _five_level_inverse_statistical_weight[OIII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[OIII][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[OIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _five_level_transition_probability[OIII][TRANSITION_0_to_1] = 2.664e-5;
    _five_level_transition_probability[OIII][TRANSITION_0_to_2] = 3.094e-11;
    _five_level_transition_probability[OIII][TRANSITION_0_to_3] = 1.69e-6;
    _five_level_transition_probability[OIII][TRANSITION_0_to_4] = 0.;
    _five_level_transition_probability[OIII][TRANSITION_1_to_2] = 9.695e-5;
    _five_level_transition_probability[OIII][TRANSITION_1_to_3] = 6.995e-3;
    _five_level_transition_probability[OIII][TRANSITION_1_to_4] = 2.268e-1;
    _five_level_transition_probability[OIII][TRANSITION_2_to_3] = 2.041e-2;
    _five_level_transition_probability[OIII][TRANSITION_2_to_4] = 6.091e-4;
    _five_level_transition_probability[OIII][TRANSITION_3_to_4] = 1.561;
    // our own fits to the data of Lennon & Burke (1994)
    // these fits were made with the script data/linecooling/gamma_OIII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[OIII][TRANSITION_0_to_1][0] = -2.517e-01;
    _five_level_collision_strength[OIII][TRANSITION_0_to_1][1] = 2.326e-03;
    _five_level_collision_strength[OIII][TRANSITION_0_to_1][2] = 2.174e+00;
    _five_level_collision_strength[OIII][TRANSITION_0_to_1][3] = -2.246e-04;
    _five_level_collision_strength[OIII][TRANSITION_0_to_1][4] = -4.349e-03;
    _five_level_collision_strength[OIII][TRANSITION_0_to_1][5] = -5.803e-06;
    _five_level_collision_strength[OIII][TRANSITION_0_to_1][6] = -4.351e-07;
    _five_level_collision_strength[OIII][TRANSITION_0_to_2][0] = -2.886e-01;
    _five_level_collision_strength[OIII][TRANSITION_0_to_2][1] = 1.462e-03;
    _five_level_collision_strength[OIII][TRANSITION_0_to_2][2] = 1.373e+00;
    _five_level_collision_strength[OIII][TRANSITION_0_to_2][3] = -1.364e-04;
    _five_level_collision_strength[OIII][TRANSITION_0_to_2][4] = 6.062e-03;
    _five_level_collision_strength[OIII][TRANSITION_0_to_2][5] = 2.338e-06;
    _five_level_collision_strength[OIII][TRANSITION_0_to_2][6] = 1.764e-07;
    _five_level_collision_strength[OIII][TRANSITION_0_to_3][0] = -7.878e-01;
    _five_level_collision_strength[OIII][TRANSITION_0_to_3][1] = -3.075e-02;
    _five_level_collision_strength[OIII][TRANSITION_0_to_3][2] = 5.015e+01;
    _five_level_collision_strength[OIII][TRANSITION_0_to_3][3] = 7.448e-03;
    _five_level_collision_strength[OIII][TRANSITION_0_to_3][4] = -3.189e-03;
    _five_level_collision_strength[OIII][TRANSITION_0_to_3][5] = 6.917e-04;
    _five_level_collision_strength[OIII][TRANSITION_0_to_3][6] = 5.198e-05;
    _five_level_collision_strength[OIII][TRANSITION_0_to_4][0] = -5.259e-01;
    _five_level_collision_strength[OIII][TRANSITION_0_to_4][1] = 1.618e-04;
    _five_level_collision_strength[OIII][TRANSITION_0_to_4][2] = 1.045e+00;
    _five_level_collision_strength[OIII][TRANSITION_0_to_4][3] = 2.345e-05;
    _five_level_collision_strength[OIII][TRANSITION_0_to_4][4] = -7.038e-04;
    _five_level_collision_strength[OIII][TRANSITION_0_to_4][5] = 3.326e-05;
    _five_level_collision_strength[OIII][TRANSITION_0_to_4][6] = 2.552e-06;
    _five_level_collision_strength[OIII][TRANSITION_1_to_2][0] = -2.714e-01;
    _five_level_collision_strength[OIII][TRANSITION_1_to_2][1] = 6.306e-03;
    _five_level_collision_strength[OIII][TRANSITION_1_to_2][2] = 5.834e+00;
    _five_level_collision_strength[OIII][TRANSITION_1_to_2][3] = -5.990e-04;
    _five_level_collision_strength[OIII][TRANSITION_1_to_2][4] = -5.646e-03;
    _five_level_collision_strength[OIII][TRANSITION_1_to_2][5] = -1.140e-05;
    _five_level_collision_strength[OIII][TRANSITION_1_to_2][6] = -8.562e-07;
    _five_level_collision_strength[OIII][TRANSITION_1_to_3][0] = -7.879e-01;
    _five_level_collision_strength[OIII][TRANSITION_1_to_3][1] = -9.239e-02;
    _five_level_collision_strength[OIII][TRANSITION_1_to_3][2] = 1.505e+02;
    _five_level_collision_strength[OIII][TRANSITION_1_to_3][3] = 2.237e-02;
    _five_level_collision_strength[OIII][TRANSITION_1_to_3][4] = 1.919e-02;
    _five_level_collision_strength[OIII][TRANSITION_1_to_3][5] = -3.453e-04;
    _five_level_collision_strength[OIII][TRANSITION_1_to_3][6] = -2.593e-05;
    _five_level_collision_strength[OIII][TRANSITION_1_to_4][0] = -5.259e-01;
    _five_level_collision_strength[OIII][TRANSITION_1_to_4][1] = 4.853e-04;
    _five_level_collision_strength[OIII][TRANSITION_1_to_4][2] = 3.136e+00;
    _five_level_collision_strength[OIII][TRANSITION_1_to_4][3] = 7.039e-05;
    _five_level_collision_strength[OIII][TRANSITION_1_to_4][4] = -1.210e-03;
    _five_level_collision_strength[OIII][TRANSITION_1_to_4][5] = 5.806e-05;
    _five_level_collision_strength[OIII][TRANSITION_1_to_4][6] = 4.455e-06;
    _five_level_collision_strength[OIII][TRANSITION_2_to_3][0] = -7.879e-01;
    _five_level_collision_strength[OIII][TRANSITION_2_to_3][1] = -1.541e-01;
    _five_level_collision_strength[OIII][TRANSITION_2_to_3][2] = 2.510e+02;
    _five_level_collision_strength[OIII][TRANSITION_2_to_3][3] = 3.730e-02;
    _five_level_collision_strength[OIII][TRANSITION_2_to_3][4] = 6.860e-02;
    _five_level_collision_strength[OIII][TRANSITION_2_to_3][5] = -1.610e-04;
    _five_level_collision_strength[OIII][TRANSITION_2_to_3][6] = -1.209e-05;
    _five_level_collision_strength[OIII][TRANSITION_2_to_4][0] = -5.260e-01;
    _five_level_collision_strength[OIII][TRANSITION_2_to_4][1] = 8.087e-04;
    _five_level_collision_strength[OIII][TRANSITION_2_to_4][2] = 5.229e+00;
    _five_level_collision_strength[OIII][TRANSITION_2_to_4][3] = 1.175e-04;
    _five_level_collision_strength[OIII][TRANSITION_2_to_4][4] = -1.729e-03;
    _five_level_collision_strength[OIII][TRANSITION_2_to_4][5] = 6.777e-05;
    _five_level_collision_strength[OIII][TRANSITION_2_to_4][6] = 5.200e-06;
    _five_level_collision_strength[OIII][TRANSITION_3_to_4][0] = 6.242e-02;
    _five_level_collision_strength[OIII][TRANSITION_3_to_4][1] = 1.702e-04;
    _five_level_collision_strength[OIII][TRANSITION_3_to_4][2] = 1.766e-01;
    _five_level_collision_strength[OIII][TRANSITION_3_to_4][3] = -1.776e-05;
    _five_level_collision_strength[OIII][TRANSITION_3_to_4][4] = -2.314e-03;
    _five_level_collision_strength[OIII][TRANSITION_3_to_4][5] = -1.121e-06;
    _five_level_collision_strength[OIII][TRANSITION_3_to_4][6] = -8.415e-08;
  }

  /// NeIII
  {
    // data from Galavis, Mendoza & Zeippen (1997), tables 6 and 7
    // ground state: 3P2
    // excited states: 3P1, 3P0, 1D2, 1S0
    // in cm^-1
    const double energy_levels[4] = {643., 921., 25841., 55751.};
    _five_level_inverse_statistical_weight[NeIII][0] = 0.2;
    _five_level_inverse_statistical_weight[NeIII][1] = 1. / 3.;
    _five_level_inverse_statistical_weight[NeIII][2] = 1.;
    _five_level_inverse_statistical_weight[NeIII][3] = 0.2;
    _five_level_inverse_statistical_weight[NeIII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[NeIII][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[NeIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // in s^-1
    _five_level_transition_probability[NeIII][TRANSITION_0_to_1] = 5.974e-3;
    _five_level_transition_probability[NeIII][TRANSITION_0_to_2] = 2.081e-8;
    _five_level_transition_probability[NeIII][TRANSITION_0_to_3] = 0.173;
    _five_level_transition_probability[NeIII][TRANSITION_0_to_4] = 3.985e-3;
    _five_level_transition_probability[NeIII][TRANSITION_1_to_2] = 1.159e-3;
    _five_level_transition_probability[NeIII][TRANSITION_1_to_3] = 5.344e-2;
    _five_level_transition_probability[NeIII][TRANSITION_1_to_4] = 2.028;
    _five_level_transition_probability[NeIII][TRANSITION_2_to_3] = 8.269e-6;
    _five_level_transition_probability[NeIII][TRANSITION_2_to_4] = 0.;
    _five_level_transition_probability[NeIII][TRANSITION_3_to_4] = 2.563;
    // our own fits to the data of Butler & Zeippen (1994)
    // these fits were made with the script data/linecooling/gamma_NeIII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[NeIII][TRANSITION_0_to_1][0] = -6.515e-01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_1][1] = 1.958e-01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_1][2] = -3.645e+01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_1][3] = -1.812e-02;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_1][4] = -1.815e-03;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_1][5] = -1.043e-03;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_1][6] = -7.735e-05;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_2][0] = -7.577e-01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_2][1] = 1.184e-01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_2][2] = -2.862e+01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_2][3] = -1.047e-02;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_2][4] = -9.100e-04;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_2][5] = -1.094e-03;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_2][6] = -8.029e-05;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_3][0] = -5.448e-03;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_3][1] = -9.752e-06;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_3][2] = 8.134e-01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_3][3] = 7.937e-07;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_3][4] = 4.823e-06;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_3][5] = 3.627e-05;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_3][6] = 3.005e-06;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_4][0] = 1.019e-01;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_4][1] = -8.068e-06;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_4][2] = 4.190e-02;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_4][3] = 8.074e-07;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_4][4] = 1.758e-08;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_4][5] = -5.000e-03;
    _five_level_collision_strength[NeIII][TRANSITION_0_to_4][6] = -3.672e-04;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_2][0] = -6.173e-01;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_2][1] = 4.546e-02;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_2][2] = -8.088e+00;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_2][3] = -4.215e-03;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_2][4] = -2.781e-04;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_2][5] = -1.537e-03;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_2][6] = -1.139e-04;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_3][0] = -1.932e-03;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_3][1] = -8.154e-06;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_3][2] = 4.758e-01;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_3][3] = 7.006e-07;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_3][4] = -9.413e-05;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_3][5] = -8.679e-07;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_3][6] = -7.422e-08;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_4][0] = 1.161e-01;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_4][1] = -5.770e-06;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_4][2] = 2.333e-02;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_4][3] = 5.826e-07;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_4][4] = -6.898e-04;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_4][5] = 1.004e-07;
    _five_level_collision_strength[NeIII][TRANSITION_1_to_4][6] = 7.445e-09;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_3][0] = -2.118e-02;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_3][1] = 7.071e-06;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_3][2] = 1.788e-01;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_3][3] = -7.768e-07;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_3][4] = 1.263e-07;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_3][5] = 1.534e-03;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_3][6] = 1.191e-04;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_4][0] = 9.843e-02;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_4][1] = -1.648e-06;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_4][2] = 8.827e-03;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_4][3] = 1.628e-07;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_4][4] = 1.792e-08;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_4][5] = -8.454e-04;
    _five_level_collision_strength[NeIII][TRANSITION_2_to_4][6] = -6.126e-05;
    _five_level_collision_strength[NeIII][TRANSITION_3_to_4][0] = 5.812e-02;
    _five_level_collision_strength[NeIII][TRANSITION_3_to_4][1] = -4.021e-05;
    _five_level_collision_strength[NeIII][TRANSITION_3_to_4][2] = 1.876e-01;
    _five_level_collision_strength[NeIII][TRANSITION_3_to_4][3] = 4.268e-06;
    _five_level_collision_strength[NeIII][TRANSITION_3_to_4][4] = 3.941e-03;
    _five_level_collision_strength[NeIII][TRANSITION_3_to_4][5] = -1.725e-07;
    _five_level_collision_strength[NeIII][TRANSITION_3_to_4][6] = -1.298e-08;
  }

  /// SII
  {
    // data from Tayal & Zatsarinny (2010), tables 1 and 3
    // ground state: 4S3/2
    // excited states: 2D3/2, 2D5/2, 2P1/2, 2P3/2
    // in eV
    const double energy_levels[4] = {1.842, 1.845, 3.041, 3.046};
    _five_level_inverse_statistical_weight[SII][0] = 0.25;
    _five_level_inverse_statistical_weight[SII][1] = 0.25;
    _five_level_inverse_statistical_weight[SII][2] = 1. / 6.;
    _five_level_inverse_statistical_weight[SII][3] = 0.5;
    _five_level_inverse_statistical_weight[SII][4] = 0.25;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[SII][TRANSITION_0_to_1] =
        energy_levels[0] * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_0_to_2] =
        energy_levels[1] * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_0_to_3] =
        energy_levels[2] * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_0_to_4] =
        energy_levels[3] * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * eV_over_k;
    _five_level_energy_difference[SII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * eV_over_k;
    // we take the sum of all contributions (in s^-1)
    _five_level_transition_probability[SII][TRANSITION_0_to_1] = 6.32e-4;
    _five_level_transition_probability[SII][TRANSITION_0_to_2] = 2.20e-4;
    _five_level_transition_probability[SII][TRANSITION_0_to_3] = 7.64e-2;
    _five_level_transition_probability[SII][TRANSITION_0_to_4] = 1.90e-1;
    _five_level_transition_probability[SII][TRANSITION_1_to_2] = 1.71e-7;
    _five_level_transition_probability[SII][TRANSITION_1_to_3] = 1.47e-1;
    _five_level_transition_probability[SII][TRANSITION_1_to_4] = 1.165e-1;
    _five_level_transition_probability[SII][TRANSITION_2_to_3] = 7.16e-2;
    _five_level_transition_probability[SII][TRANSITION_2_to_4] = 1.61e-1;
    _five_level_transition_probability[SII][TRANSITION_3_to_4] = 2.43e-7;
    // our own fits to the data of Tayal & Zatsarinny (2010)
    // these fits were made with the script data/linecooling/gamma_SII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[SII][TRANSITION_0_to_1][0] = -2.378e-02;
    _five_level_collision_strength[SII][TRANSITION_0_to_1][1] = -3.780e-05;
    _five_level_collision_strength[SII][TRANSITION_0_to_1][2] = 3.350e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_1][3] = 2.326e-06;
    _five_level_collision_strength[SII][TRANSITION_0_to_1][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_1][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_1][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_2][0] = -2.591e-02;
    _five_level_collision_strength[SII][TRANSITION_0_to_2][1] = -5.072e-05;
    _five_level_collision_strength[SII][TRANSITION_0_to_2][2] = 5.088e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_2][3] = 3.003e-06;
    _five_level_collision_strength[SII][TRANSITION_0_to_2][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_2][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_2][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_3][0] = -4.744e-01;
    _five_level_collision_strength[SII][TRANSITION_0_to_3][1] = 1.594e-02;
    _five_level_collision_strength[SII][TRANSITION_0_to_3][2] = 1.361e+01;
    _five_level_collision_strength[SII][TRANSITION_0_to_3][3] = -1.274e-03;
    _five_level_collision_strength[SII][TRANSITION_0_to_3][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_3][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_3][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_4][0] = -4.730e-01;
    _five_level_collision_strength[SII][TRANSITION_0_to_4][1] = 3.177e-02;
    _five_level_collision_strength[SII][TRANSITION_0_to_4][2] = 2.699e+01;
    _five_level_collision_strength[SII][TRANSITION_0_to_4][3] = -2.544e-03;
    _five_level_collision_strength[SII][TRANSITION_0_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_0_to_4][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_2][0] = -8.272e-02;
    _five_level_collision_strength[SII][TRANSITION_1_to_2][1] = 3.140e-05;
    _five_level_collision_strength[SII][TRANSITION_1_to_2][2] = 1.492e+01;
    _five_level_collision_strength[SII][TRANSITION_1_to_2][3] = -5.602e-06;
    _five_level_collision_strength[SII][TRANSITION_1_to_2][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_2][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_2][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_3][0] = -5.368e-02;
    _five_level_collision_strength[SII][TRANSITION_1_to_3][1] = 8.155e-05;
    _five_level_collision_strength[SII][TRANSITION_1_to_3][2] = 2.216e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_3][3] = -6.673e-06;
    _five_level_collision_strength[SII][TRANSITION_1_to_3][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_3][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_3][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_4][0] = -2.183e-02;
    _five_level_collision_strength[SII][TRANSITION_1_to_4][1] = 5.571e-05;
    _five_level_collision_strength[SII][TRANSITION_1_to_4][2] = 2.822e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_4][3] = -4.915e-06;
    _five_level_collision_strength[SII][TRANSITION_1_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_1_to_4][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_3][0] = -1.569e-02;
    _five_level_collision_strength[SII][TRANSITION_2_to_3][1] = 2.911e-05;
    _five_level_collision_strength[SII][TRANSITION_2_to_3][2] = 2.012e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_3][3] = -2.687e-06;
    _five_level_collision_strength[SII][TRANSITION_2_to_3][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_3][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_3][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_4][0] = -3.268e-02;
    _five_level_collision_strength[SII][TRANSITION_2_to_4][1] = 1.339e-04;
    _five_level_collision_strength[SII][TRANSITION_2_to_4][2] = 5.190e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_4][3] = -1.130e-05;
    _five_level_collision_strength[SII][TRANSITION_2_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_2_to_4][6] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_3_to_4][0] = -5.249e-01;
    _five_level_collision_strength[SII][TRANSITION_3_to_4][1] = 6.595e-02;
    _five_level_collision_strength[SII][TRANSITION_3_to_4][2] = 4.606e+01;
    _five_level_collision_strength[SII][TRANSITION_3_to_4][3] = -5.206e-03;
    _five_level_collision_strength[SII][TRANSITION_3_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_3_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SII][TRANSITION_3_to_4][6] = 0.000e+00;
  }

  /// SIII
  {
    // data from Mendoza & Zeippen (1982), table 2 and 3
    // ground state: 3P0
    // excited states: 3P1, 3P2, 1D2, 1S0
    // in Ry
    const double energy_levels[4] = {0.002708, 0.007586, 0.103157, 0.247532};
    _five_level_inverse_statistical_weight[SIII][0] = 1.;
    _five_level_inverse_statistical_weight[SIII][1] = 1. / 3.;
    _five_level_inverse_statistical_weight[SIII][2] = 0.2;
    _five_level_inverse_statistical_weight[SIII][3] = 0.2;
    _five_level_inverse_statistical_weight[SIII][4] = 1.;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[SIII][TRANSITION_0_to_1] =
        energy_levels[0] * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_0_to_2] =
        energy_levels[1] * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_0_to_3] =
        energy_levels[2] * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_0_to_4] =
        energy_levels[3] * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * Ry_over_k;
    _five_level_energy_difference[SIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * Ry_over_k;
    // we take the sum of all contributions (in s^-1)
    _five_level_transition_probability[SIII][TRANSITION_0_to_1] = 4.72e-4;
    _five_level_transition_probability[SIII][TRANSITION_0_to_2] = 4.61e-8;
    _five_level_transition_probability[SIII][TRANSITION_0_to_3] = 5.82e-6;
    _five_level_transition_probability[SIII][TRANSITION_0_to_4] = 0.;
    _five_level_transition_probability[SIII][TRANSITION_1_to_2] = 2.07e-3;
    _five_level_transition_probability[SIII][TRANSITION_1_to_3] = 2.20e-2;
    _five_level_transition_probability[SIII][TRANSITION_1_to_4] = 7.96e-1;
    _five_level_transition_probability[SIII][TRANSITION_2_to_3] = 5.76e-2;
    _five_level_transition_probability[SIII][TRANSITION_2_to_4] = 1.05e-2;
    _five_level_transition_probability[SIII][TRANSITION_3_to_4] = 2.22;
    // our own fits to the data of Hudson, Ramsbottom & Scott (2012)
    // these fits were made with the script data/linecooling/gamma_SIII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[SIII][TRANSITION_0_to_1][0] = 8.018e-02;
    _five_level_collision_strength[SIII][TRANSITION_0_to_1][1] = -4.930e-05;
    _five_level_collision_strength[SIII][TRANSITION_0_to_1][2] = 1.206e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_1][3] = 3.800e-06;
    _five_level_collision_strength[SIII][TRANSITION_0_to_1][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_1][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_1][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_2][0] = -2.250e-01;
    _five_level_collision_strength[SIII][TRANSITION_0_to_2][1] = 1.746e-03;
    _five_level_collision_strength[SIII][TRANSITION_0_to_2][2] = 3.814e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_2][3] = -1.408e-04;
    _five_level_collision_strength[SIII][TRANSITION_0_to_2][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_2][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_2][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_3][0] = -3.470e-02;
    _five_level_collision_strength[SIII][TRANSITION_0_to_3][1] = 6.561e-05;
    _five_level_collision_strength[SIII][TRANSITION_0_to_3][2] = 8.851e-01;
    _five_level_collision_strength[SIII][TRANSITION_0_to_3][3] = -5.617e-06;
    _five_level_collision_strength[SIII][TRANSITION_0_to_3][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_3][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_3][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_4][0] = 1.134e-01;
    _five_level_collision_strength[SIII][TRANSITION_0_to_4][1] = 3.383e-06;
    _five_level_collision_strength[SIII][TRANSITION_0_to_4][2] = 3.649e-02;
    _five_level_collision_strength[SIII][TRANSITION_0_to_4][3] = -2.851e-07;
    _five_level_collision_strength[SIII][TRANSITION_0_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_0_to_4][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_2][0] = -4.587e-03;
    _five_level_collision_strength[SIII][TRANSITION_1_to_2][1] = 2.257e-04;
    _five_level_collision_strength[SIII][TRANSITION_1_to_2][2] = 4.846e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_2][3] = -1.949e-05;
    _five_level_collision_strength[SIII][TRANSITION_1_to_2][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_2][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_2][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_3][0] = -4.454e-02;
    _five_level_collision_strength[SIII][TRANSITION_1_to_3][1] = 2.365e-04;
    _five_level_collision_strength[SIII][TRANSITION_1_to_3][2] = 2.823e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_3][3] = -2.014e-05;
    _five_level_collision_strength[SIII][TRANSITION_1_to_3][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_3][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_3][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_4][0] = 1.665e-01;
    _five_level_collision_strength[SIII][TRANSITION_1_to_4][1] = 2.041e-06;
    _five_level_collision_strength[SIII][TRANSITION_1_to_4][2] = 6.825e-02;
    _five_level_collision_strength[SIII][TRANSITION_1_to_4][3] = -1.900e-07;
    _five_level_collision_strength[SIII][TRANSITION_1_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_1_to_4][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_2_to_3][0] = -3.303e-02;
    _five_level_collision_strength[SIII][TRANSITION_2_to_3][1] = 3.185e-04;
    _five_level_collision_strength[SIII][TRANSITION_2_to_3][2] = 4.890e+00;
    _five_level_collision_strength[SIII][TRANSITION_2_to_3][3] = -2.742e-05;
    _five_level_collision_strength[SIII][TRANSITION_2_to_3][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_2_to_3][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_2_to_3][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_2_to_4][0] = 1.765e-01;
    _five_level_collision_strength[SIII][TRANSITION_2_to_4][1] = 2.565e-06;
    _five_level_collision_strength[SIII][TRANSITION_2_to_4][2] = 1.041e-01;
    _five_level_collision_strength[SIII][TRANSITION_2_to_4][3] = -2.451e-07;
    _five_level_collision_strength[SIII][TRANSITION_2_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_2_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_2_to_4][6] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_3_to_4][0] = -2.273e-01;
    _five_level_collision_strength[SIII][TRANSITION_3_to_4][1] = 3.338e-03;
    _five_level_collision_strength[SIII][TRANSITION_3_to_4][2] = 2.390e+00;
    _five_level_collision_strength[SIII][TRANSITION_3_to_4][3] = -2.696e-04;
    _five_level_collision_strength[SIII][TRANSITION_3_to_4][4] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_3_to_4][5] = 0.000e+00;
    _five_level_collision_strength[SIII][TRANSITION_3_to_4][6] = 0.000e+00;
  }

  /// CII
  {
    // data from Froese Fischer & Tachiev (2004), table 2
    // ground state: 2P1/2
    // excited states: 2P3/2, 4P1/2, 4P3/2, 4P5/2
    // in cm^-1
    const double energy_levels[4] = {63.67, 43057.99, 43080.21, 43108.66};
    _five_level_inverse_statistical_weight[CII][0] = 0.5;
    _five_level_inverse_statistical_weight[CII][1] = 0.25;
    _five_level_inverse_statistical_weight[CII][2] = 0.5;
    _five_level_inverse_statistical_weight[CII][3] = 0.25;
    _five_level_inverse_statistical_weight[CII][4] = 1. / 6.;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[CII][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[CII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we sum contributions of all types (in s^-1)
    _five_level_transition_probability[CII][TRANSITION_0_to_1] = 2.321e-6;
    _five_level_transition_probability[CII][TRANSITION_0_to_2] = 6.136e1;
    _five_level_transition_probability[CII][TRANSITION_0_to_3] = 1.463;
    _five_level_transition_probability[CII][TRANSITION_0_to_4] = 8.177e-4;
    _five_level_transition_probability[CII][TRANSITION_1_to_2] = 6.929e1;
    _five_level_transition_probability[CII][TRANSITION_1_to_3] = 8.853;
    _five_level_transition_probability[CII][TRANSITION_1_to_4] = 4.477e1;
    _five_level_transition_probability[CII][TRANSITION_2_to_3] = 2.467e-7;
    _five_level_transition_probability[CII][TRANSITION_2_to_4] = 3.571e-14;
    _five_level_transition_probability[CII][TRANSITION_3_to_4] = 3.725e-7;
    // our own fits to the data of Tayal (2008)
    // these fits were made with the script data/linecooling/gamma_CII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[CII][TRANSITION_0_to_1][0] = -9.519e-01;
    _five_level_collision_strength[CII][TRANSITION_0_to_1][1] = -5.097e+00;
    _five_level_collision_strength[CII][TRANSITION_0_to_1][2] = 1.436e+03;
    _five_level_collision_strength[CII][TRANSITION_0_to_1][3] = 7.424e-01;
    _five_level_collision_strength[CII][TRANSITION_0_to_1][4] = -5.391e+01;
    _five_level_collision_strength[CII][TRANSITION_0_to_1][5] = 3.673e-06;
    _five_level_collision_strength[CII][TRANSITION_0_to_1][6] = 2.806e-07;
    _five_level_collision_strength[CII][TRANSITION_0_to_2][0] = 1.247e-01;
    _five_level_collision_strength[CII][TRANSITION_0_to_2][1] = -2.710e-05;
    _five_level_collision_strength[CII][TRANSITION_0_to_2][2] = 1.110e-01;
    _five_level_collision_strength[CII][TRANSITION_0_to_2][3] = 2.784e-06;
    _five_level_collision_strength[CII][TRANSITION_0_to_2][4] = 1.410e-02;
    _five_level_collision_strength[CII][TRANSITION_0_to_2][5] = -3.198e-08;
    _five_level_collision_strength[CII][TRANSITION_0_to_2][6] = -2.437e-09;
    _five_level_collision_strength[CII][TRANSITION_0_to_3][0] = 3.496e-07;
    _five_level_collision_strength[CII][TRANSITION_0_to_3][1] = -2.076e-05;
    _five_level_collision_strength[CII][TRANSITION_0_to_3][2] = 3.990e-01;
    _five_level_collision_strength[CII][TRANSITION_0_to_3][3] = 2.426e-06;
    _five_level_collision_strength[CII][TRANSITION_0_to_3][4] = -2.832e+00;
    _five_level_collision_strength[CII][TRANSITION_0_to_3][5] = 2.765e-10;
    _five_level_collision_strength[CII][TRANSITION_0_to_3][6] = 2.159e-11;
    _five_level_collision_strength[CII][TRANSITION_0_to_4][0] = -2.814e-06;
    _five_level_collision_strength[CII][TRANSITION_0_to_4][1] = -2.407e-05;
    _five_level_collision_strength[CII][TRANSITION_0_to_4][2] = 2.474e-01;
    _five_level_collision_strength[CII][TRANSITION_0_to_4][3] = 2.906e-06;
    _five_level_collision_strength[CII][TRANSITION_0_to_4][4] = 1.393e+00;
    _five_level_collision_strength[CII][TRANSITION_0_to_4][5] = -6.962e-10;
    _five_level_collision_strength[CII][TRANSITION_0_to_4][6] = -5.438e-11;
    _five_level_collision_strength[CII][TRANSITION_1_to_2][0] = 1.949e-01;
    _five_level_collision_strength[CII][TRANSITION_1_to_2][1] = -9.811e-06;
    _five_level_collision_strength[CII][TRANSITION_1_to_2][2] = 4.074e-02;
    _five_level_collision_strength[CII][TRANSITION_1_to_2][3] = 1.003e-06;
    _five_level_collision_strength[CII][TRANSITION_1_to_2][4] = 7.165e-03;
    _five_level_collision_strength[CII][TRANSITION_1_to_2][5] = -2.276e-08;
    _five_level_collision_strength[CII][TRANSITION_1_to_2][6] = -1.738e-09;
    _five_level_collision_strength[CII][TRANSITION_1_to_3][0] = 1.449e-01;
    _five_level_collision_strength[CII][TRANSITION_1_to_3][1] = -3.673e-05;
    _five_level_collision_strength[CII][TRANSITION_1_to_3][2] = 1.713e-01;
    _five_level_collision_strength[CII][TRANSITION_1_to_3][3] = 3.763e-06;
    _five_level_collision_strength[CII][TRANSITION_1_to_3][4] = -9.812e-03;
    _five_level_collision_strength[CII][TRANSITION_1_to_3][5] = 6.304e-08;
    _five_level_collision_strength[CII][TRANSITION_1_to_3][6] = 4.812e-09;
    _five_level_collision_strength[CII][TRANSITION_1_to_4][0] = -2.743e-06;
    _five_level_collision_strength[CII][TRANSITION_1_to_4][1] = -1.291e-04;
    _five_level_collision_strength[CII][TRANSITION_1_to_4][2] = 1.183e+00;
    _five_level_collision_strength[CII][TRANSITION_1_to_4][3] = 1.437e-05;
    _five_level_collision_strength[CII][TRANSITION_1_to_4][4] = 1.598e+00;
    _five_level_collision_strength[CII][TRANSITION_1_to_4][5] = -2.263e-09;
    _five_level_collision_strength[CII][TRANSITION_1_to_4][6] = -1.753e-10;
    _five_level_collision_strength[CII][TRANSITION_2_to_3][0] = -1.195e-03;
    _five_level_collision_strength[CII][TRANSITION_2_to_3][1] = -2.027e-04;
    _five_level_collision_strength[CII][TRANSITION_2_to_3][2] = 6.479e-01;
    _five_level_collision_strength[CII][TRANSITION_2_to_3][3] = 2.660e-05;
    _five_level_collision_strength[CII][TRANSITION_2_to_3][4] = -7.636e-01;
    _five_level_collision_strength[CII][TRANSITION_2_to_3][5] = 1.257e-08;
    _five_level_collision_strength[CII][TRANSITION_2_to_3][6] = 9.795e-10;
    _five_level_collision_strength[CII][TRANSITION_2_to_4][0] = -4.787e-01;
    _five_level_collision_strength[CII][TRANSITION_2_to_4][1] = -1.490e-03;
    _five_level_collision_strength[CII][TRANSITION_2_to_4][2] = 1.501e+01;
    _five_level_collision_strength[CII][TRANSITION_2_to_4][3] = 9.750e-04;
    _five_level_collision_strength[CII][TRANSITION_2_to_4][4] = 3.103e+00;
    _five_level_collision_strength[CII][TRANSITION_2_to_4][5] = -2.437e-07;
    _five_level_collision_strength[CII][TRANSITION_2_to_4][6] = -1.908e-08;
    _five_level_collision_strength[CII][TRANSITION_3_to_4][0] = -2.430e-03;
    _five_level_collision_strength[CII][TRANSITION_3_to_4][1] = -3.839e-04;
    _five_level_collision_strength[CII][TRANSITION_3_to_4][2] = 1.605e+00;
    _five_level_collision_strength[CII][TRANSITION_3_to_4][3] = 5.173e-05;
    _five_level_collision_strength[CII][TRANSITION_3_to_4][4] = -2.446e-01;
    _five_level_collision_strength[CII][TRANSITION_3_to_4][5] = 8.372e-08;
    _five_level_collision_strength[CII][TRANSITION_3_to_4][6] = 6.553e-09;
  }

  /// CIII
  {
    // data from Froese Fischer & Tachiev (2004), table 1
    // ground state: 1S0
    // excited states: 3P0, 3P1, 3P2, 1P1
    // in cm^-1
    const double energy_levels[4] = {52391.87, 52415.53, 52472.16, 102446.94};
    _five_level_inverse_statistical_weight[CIII][0] = 1.;
    _five_level_inverse_statistical_weight[CIII][1] = 1.;
    _five_level_inverse_statistical_weight[CIII][2] = 1. / 3.;
    _five_level_inverse_statistical_weight[CIII][3] = 0.2;
    _five_level_inverse_statistical_weight[CIII][4] = 1. / 3.;
    // convert energy levels to energy differences (and convert units)
    _five_level_energy_difference[CIII][TRANSITION_0_to_1] =
        energy_levels[0] * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_0_to_2] =
        energy_levels[1] * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_0_to_3] =
        energy_levels[2] * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_0_to_4] =
        energy_levels[3] * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_1_to_2] =
        (energy_levels[1] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_1_to_3] =
        (energy_levels[2] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_1_to_4] =
        (energy_levels[3] - energy_levels[0]) * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_2_to_3] =
        (energy_levels[2] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_2_to_4] =
        (energy_levels[3] - energy_levels[1]) * hc_over_k;
    _five_level_energy_difference[CIII][TRANSITION_3_to_4] =
        (energy_levels[3] - energy_levels[2]) * hc_over_k;
    // we sum contributions of all types (in s^-1)
    _five_level_transition_probability[CIII][TRANSITION_0_to_1] = 0.;
    _five_level_transition_probability[CIII][TRANSITION_0_to_2] = 1.040e2;
    _five_level_transition_probability[CIII][TRANSITION_0_to_3] = 5.216e-3;
    _five_level_transition_probability[CIII][TRANSITION_0_to_4] = 1.769e9;
    _five_level_transition_probability[CIII][TRANSITION_1_to_2] = 2.381e-7;
    _five_level_transition_probability[CIII][TRANSITION_1_to_3] = 1.486e-13;
    _five_level_transition_probability[CIII][TRANSITION_1_to_4] = 1.595e-3;
    _five_level_transition_probability[CIII][TRANSITION_2_to_3] = 2.450e-6;
    _five_level_transition_probability[CIII][TRANSITION_2_to_4] = 1.269e-3;
    _five_level_transition_probability[CIII][TRANSITION_3_to_4] = 2.036e-3;
    // our own fits to the data of Berrington et al. (1985)
    // these fits were made with the script data/linecooling/gamma_CIII.py
    // (this script also outputs the code below)
    _five_level_collision_strength[CIII][TRANSITION_0_to_1][0] = 1.265e-07;
    _five_level_collision_strength[CIII][TRANSITION_0_to_1][1] = 5.762e-06;
    _five_level_collision_strength[CIII][TRANSITION_0_to_1][2] = 1.080e-01;
    _five_level_collision_strength[CIII][TRANSITION_0_to_1][3] = -5.626e-07;
    _five_level_collision_strength[CIII][TRANSITION_0_to_1][4] = 1.412e+00;
    _five_level_collision_strength[CIII][TRANSITION_0_to_1][5] = 1.950e-11;
    _five_level_collision_strength[CIII][TRANSITION_0_to_1][6] = 1.323e-12;
    _five_level_collision_strength[CIII][TRANSITION_0_to_2][0] = 3.793e-07;
    _five_level_collision_strength[CIII][TRANSITION_0_to_2][1] = 1.727e-05;
    _five_level_collision_strength[CIII][TRANSITION_0_to_2][2] = 3.239e-01;
    _five_level_collision_strength[CIII][TRANSITION_0_to_2][3] = -1.686e-06;
    _five_level_collision_strength[CIII][TRANSITION_0_to_2][4] = 1.095e+00;
    _five_level_collision_strength[CIII][TRANSITION_0_to_2][5] = 7.533e-11;
    _five_level_collision_strength[CIII][TRANSITION_0_to_2][6] = 5.111e-12;
    _five_level_collision_strength[CIII][TRANSITION_0_to_3][0] = 3.942e-07;
    _five_level_collision_strength[CIII][TRANSITION_0_to_3][1] = 2.877e-05;
    _five_level_collision_strength[CIII][TRANSITION_0_to_3][2] = 5.399e-01;
    _five_level_collision_strength[CIII][TRANSITION_0_to_3][3] = -2.809e-06;
    _five_level_collision_strength[CIII][TRANSITION_0_to_3][4] = 6.095e-01;
    _five_level_collision_strength[CIII][TRANSITION_0_to_3][5] = 2.255e-10;
    _five_level_collision_strength[CIII][TRANSITION_0_to_3][6] = 1.530e-11;
    _five_level_collision_strength[CIII][TRANSITION_0_to_4][0] = 8.424e-06;
    _five_level_collision_strength[CIII][TRANSITION_0_to_4][1] = 7.004e-05;
    _five_level_collision_strength[CIII][TRANSITION_0_to_4][2] = 3.767e+00;
    _five_level_collision_strength[CIII][TRANSITION_0_to_4][3] = -5.480e-06;
    _five_level_collision_strength[CIII][TRANSITION_0_to_4][4] = -4.018e-01;
    _five_level_collision_strength[CIII][TRANSITION_0_to_4][5] = -4.228e-10;
    _five_level_collision_strength[CIII][TRANSITION_0_to_4][6] = -2.906e-11;
    _five_level_collision_strength[CIII][TRANSITION_1_to_2][0] = 5.033e-06;
    _five_level_collision_strength[CIII][TRANSITION_1_to_2][1] = 1.611e-04;
    _five_level_collision_strength[CIII][TRANSITION_1_to_2][2] = 6.600e-01;
    _five_level_collision_strength[CIII][TRANSITION_1_to_2][3] = -1.477e-05;
    _five_level_collision_strength[CIII][TRANSITION_1_to_2][4] = 3.546e+00;
    _five_level_collision_strength[CIII][TRANSITION_1_to_2][5] = 1.577e-10;
    _five_level_collision_strength[CIII][TRANSITION_1_to_2][6] = 1.063e-11;
    _five_level_collision_strength[CIII][TRANSITION_1_to_3][0] = -6.251e-05;
    _five_level_collision_strength[CIII][TRANSITION_1_to_3][1] = 3.034e-04;
    _five_level_collision_strength[CIII][TRANSITION_1_to_3][2] = 1.039e-01;
    _five_level_collision_strength[CIII][TRANSITION_1_to_3][3] = -2.744e-05;
    _five_level_collision_strength[CIII][TRANSITION_1_to_3][4] = 2.856e-01;
    _five_level_collision_strength[CIII][TRANSITION_1_to_3][5] = 3.444e-09;
    _five_level_collision_strength[CIII][TRANSITION_1_to_3][6] = 2.317e-10;
    _five_level_collision_strength[CIII][TRANSITION_1_to_4][0] = 1.741e-06;
    _five_level_collision_strength[CIII][TRANSITION_1_to_4][1] = -1.148e-05;
    _five_level_collision_strength[CIII][TRANSITION_1_to_4][2] = 4.633e-01;
    _five_level_collision_strength[CIII][TRANSITION_1_to_4][3] = 7.853e-07;
    _five_level_collision_strength[CIII][TRANSITION_1_to_4][4] = -4.244e+00;
    _five_level_collision_strength[CIII][TRANSITION_1_to_4][5] = -5.205e-12;
    _five_level_collision_strength[CIII][TRANSITION_1_to_4][6] = -3.781e-13;
    _five_level_collision_strength[CIII][TRANSITION_2_to_3][0] = 1.932e-05;
    _five_level_collision_strength[CIII][TRANSITION_2_to_3][1] = 8.772e-04;
    _five_level_collision_strength[CIII][TRANSITION_2_to_3][2] = 1.032e+00;
    _five_level_collision_strength[CIII][TRANSITION_2_to_3][3] = -7.960e-05;
    _five_level_collision_strength[CIII][TRANSITION_2_to_3][4] = 3.735e+00;
    _five_level_collision_strength[CIII][TRANSITION_2_to_3][5] = 7.737e-10;
    _five_level_collision_strength[CIII][TRANSITION_2_to_3][6] = 5.206e-11;
    _five_level_collision_strength[CIII][TRANSITION_2_to_4][0] = 3.777e-04;
    _five_level_collision_strength[CIII][TRANSITION_2_to_4][1] = -3.470e-05;
    _five_level_collision_strength[CIII][TRANSITION_2_to_4][2] = 1.386e+00;
    _five_level_collision_strength[CIII][TRANSITION_2_to_4][3] = 2.382e-06;
    _five_level_collision_strength[CIII][TRANSITION_2_to_4][4] = 4.838e-02;
    _five_level_collision_strength[CIII][TRANSITION_2_to_4][5] = 1.336e-09;
    _five_level_collision_strength[CIII][TRANSITION_2_to_4][6] = 9.718e-11;
    _five_level_collision_strength[CIII][TRANSITION_3_to_4][0] = 7.721e-06;
    _five_level_collision_strength[CIII][TRANSITION_3_to_4][1] = -5.752e-05;
    _five_level_collision_strength[CIII][TRANSITION_3_to_4][2] = 2.317e+00;
    _five_level_collision_strength[CIII][TRANSITION_3_to_4][3] = 3.938e-06;
    _five_level_collision_strength[CIII][TRANSITION_3_to_4][4] = 5.449e-01;
    _five_level_collision_strength[CIII][TRANSITION_3_to_4][5] = 2.013e-10;
    _five_level_collision_strength[CIII][TRANSITION_3_to_4][6] = 1.463e-11;
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
    _two_level_energy_difference[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS] =
        0.00159 * Ry_over_k;
    // Galavis, Mendoza & Zeippen (1998), table 4, 1 to 2 transition (in s^-1)
    _two_level_transition_probability[NIII -
                                      LINECOOLINGDATA_NUMFIVELEVELELEMENTS] =
        4.736e-5;
    // our own fits to the data of Blum & Pradhan (1992), table 3
    // these fits were made with the script data/linecooling/gamma_NIII.py
    // (this script also outputs the code below)
    _two_level_collision_strength[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [0] = -3.257e-02;
    _two_level_collision_strength[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [1] = -3.573e-04;
    _two_level_collision_strength[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [2] = 1.653e+00;
    _two_level_collision_strength[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [3] = 4.822e-05;
    _two_level_collision_strength[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [4] = -7.529e-06;
    _two_level_collision_strength[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [5] = 3.018e-03;
    _two_level_collision_strength[NIII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [6] = 2.453e-04;
    // statistical weights: level 0 is a P_{1/2} level, while level 1 is a
    // P_{3/2}
    _two_level_inverse_statistical_weight[NIII -
                                          LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                         [0] = 0.5;
    _two_level_inverse_statistical_weight[NIII -
                                          LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                         [1] = 0.25;
  }

  /// NeII
  {
    // Saraph & Tully (1994), table 2, fine structure splitting energy for
    // Z = 10 (in Ry)
    // ground state: 2P3/2
    // excited state: 2P1/2
    _two_level_energy_difference[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS] =
        0.0071 * Ry_over_k;
    // Kaufman & Sugar (1986), table 7 (in s^-1)
    _two_level_transition_probability[NeII -
                                      LINECOOLINGDATA_NUMFIVELEVELELEMENTS] =
        8.55e-3;
    // our own fits to the data of Griffin, Mitnik & Badnell (2001), table 4
    // these fits were made with the script data/linecooling/gamma_NeII.py
    // (this script also outputs the code below)
    _two_level_collision_strength[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [0] = -8.548e-01;
    _two_level_collision_strength[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [1] = 1.938e-01;
    _two_level_collision_strength[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [2] = -1.038e+01;
    _two_level_collision_strength[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [3] = -1.254e-02;
    _two_level_collision_strength[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [4] = 2.715e-02;
    _two_level_collision_strength[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [5] = 5.755e-05;
    _two_level_collision_strength[NeII - LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                 [6] = 4.171e-06;
    // statistical weights: level 0 is a P_{3/2} level, while level 1 is a
    // P_{1/2}
    _two_level_inverse_statistical_weight[NeII -
                                          LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                         [0] = 0.25;
    _two_level_inverse_statistical_weight[NeII -
                                          LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                         [1] = 0.5;
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
  return _five_level_transition_probability[element][transition];
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
  return _five_level_energy_difference[element][transition];
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
  return 1. / _five_level_inverse_statistical_weight[element][level];
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
int LineCoolingData::solve_system_of_linear_equations(double A[5][5],
                                                      double B[5]) {

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
 * @param logT Natural logarithm of the temperature in K.
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
        std::pow(T, 1. + _five_level_collision_strength[element][i][0]) *
        (_five_level_collision_strength[element][i][1] +
         _five_level_collision_strength[element][i][2] * Tinv +
         _five_level_collision_strength[element][i][3] * logT +
         _five_level_collision_strength[element][i][4] * T *
             (1. +
              (_five_level_collision_strength[element][i][5] - 1.) *
                  std::pow(T, _five_level_collision_strength[element][i][6])));
    collision_rate_down[i] = collision_strength;
    collision_rate_up[i] =
        collision_strength *
        std::exp(-_five_level_energy_difference[element][i] * Tinv);
  }

  level_matrix[1][0] = collision_rate_up[TRANSITION_0_to_1] *
                       _five_level_inverse_statistical_weight[element][0];
  level_matrix[1][1] =
      -(_five_level_transition_probability[element][TRANSITION_0_to_1] +
        _five_level_inverse_statistical_weight[element][1] *
            (collision_rate_down[TRANSITION_0_to_1] +
             collision_rate_up[TRANSITION_1_to_2] +
             collision_rate_up[TRANSITION_1_to_3] +
             collision_rate_up[TRANSITION_1_to_4]));
  level_matrix[1][2] =
      _five_level_transition_probability[element][TRANSITION_1_to_2] +
      _five_level_inverse_statistical_weight[element][2] *
          collision_rate_down[TRANSITION_1_to_2];
  level_matrix[1][3] =
      _five_level_transition_probability[element][TRANSITION_1_to_3] +
      _five_level_inverse_statistical_weight[element][3] *
          collision_rate_down[TRANSITION_1_to_3];
  level_matrix[1][4] =
      _five_level_transition_probability[element][TRANSITION_1_to_4] +
      _five_level_inverse_statistical_weight[element][4] *
          collision_rate_down[TRANSITION_1_to_4];

  level_matrix[2][0] = collision_rate_up[TRANSITION_0_to_2] *
                       _five_level_inverse_statistical_weight[element][0];
  level_matrix[2][1] = collision_rate_up[TRANSITION_1_to_2] *
                       _five_level_inverse_statistical_weight[element][1];
  level_matrix[2][2] =
      -(_five_level_transition_probability[element][TRANSITION_0_to_2] +
        _five_level_transition_probability[element][TRANSITION_1_to_2] +
        _five_level_inverse_statistical_weight[element][2] *
            (collision_rate_down[TRANSITION_0_to_2] +
             collision_rate_down[TRANSITION_1_to_2] +
             collision_rate_up[TRANSITION_2_to_3] +
             collision_rate_up[TRANSITION_2_to_4]));
  level_matrix[2][3] =
      _five_level_transition_probability[element][TRANSITION_2_to_3] +
      collision_rate_down[TRANSITION_2_to_3] *
          _five_level_inverse_statistical_weight[element][3];
  level_matrix[2][4] =
      _five_level_transition_probability[element][TRANSITION_2_to_4] +
      collision_rate_down[TRANSITION_2_to_4] *
          _five_level_inverse_statistical_weight[element][4];

  level_matrix[3][0] = collision_rate_up[TRANSITION_0_to_3] *
                       _five_level_inverse_statistical_weight[element][0];
  level_matrix[3][1] = collision_rate_up[TRANSITION_1_to_3] *
                       _five_level_inverse_statistical_weight[element][1];
  level_matrix[3][2] = collision_rate_up[TRANSITION_2_to_3] *
                       _five_level_inverse_statistical_weight[element][2];
  level_matrix[3][3] =
      -(_five_level_transition_probability[element][TRANSITION_0_to_3] +
        _five_level_transition_probability[element][TRANSITION_1_to_3] +
        _five_level_transition_probability[element][TRANSITION_2_to_3] +
        _five_level_inverse_statistical_weight[element][3] *
            (collision_rate_down[TRANSITION_0_to_3] +
             collision_rate_down[TRANSITION_1_to_3] +
             collision_rate_down[TRANSITION_2_to_3] +
             collision_rate_up[TRANSITION_3_to_4]));
  level_matrix[3][4] =
      _five_level_transition_probability[element][TRANSITION_3_to_4] +
      collision_rate_down[TRANSITION_3_to_4] *
          _five_level_inverse_statistical_weight[element][4];

  level_matrix[4][0] = collision_rate_up[TRANSITION_0_to_4] *
                       _five_level_inverse_statistical_weight[element][0];
  level_matrix[4][1] = collision_rate_up[TRANSITION_1_to_4] *
                       _five_level_inverse_statistical_weight[element][1];
  level_matrix[4][2] = collision_rate_up[TRANSITION_2_to_4] *
                       _five_level_inverse_statistical_weight[element][2];
  level_matrix[4][3] = collision_rate_up[TRANSITION_3_to_4] *
                       _five_level_inverse_statistical_weight[element][3];
  level_matrix[4][4] =
      -(_five_level_transition_probability[element][TRANSITION_0_to_4] +
        _five_level_transition_probability[element][TRANSITION_1_to_4] +
        _five_level_transition_probability[element][TRANSITION_2_to_4] +
        _five_level_transition_probability[element][TRANSITION_3_to_4] +
        _five_level_inverse_statistical_weight[element][4] *
            (collision_rate_down[TRANSITION_0_to_4] +
             collision_rate_down[TRANSITION_1_to_4] +
             collision_rate_down[TRANSITION_2_to_4] +
             collision_rate_down[TRANSITION_3_to_4]));

  // find level populations
  const int status =
      solve_system_of_linear_equations(level_matrix, level_populations);
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
 * @param T Temperature (in K).
 * @param Tinv Inverse temperature (in K^-1).
 * @param logT Natural logarithm of the temperature in K.
 * @return Level population of the second level.
 */
double LineCoolingData::compute_level_population(
    LineCoolingDataTwoLevelElement element, double collision_strength_prefactor,
    double T, double Tinv, double logT) const {

  // note that we need to remap the element index
  const int i = element - LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  const double ksi = _two_level_energy_difference[i];
  const double A = _two_level_transition_probability[i];
  const double collision_strength =
      collision_strength_prefactor *
      std::pow(T, 1. + _two_level_collision_strength[i][0]) *
      (_two_level_collision_strength[i][1] +
       _two_level_collision_strength[i][2] * Tinv +
       _two_level_collision_strength[i][3] * logT +
       _two_level_collision_strength[i][4] * T *
           (1. +
            (_two_level_collision_strength[i][5] - 1.) *
                std::pow(T, _two_level_collision_strength[i][6])));
  const double inv_omega_1 = _two_level_inverse_statistical_weight[i][0];
  const double inv_omega_2 = _two_level_inverse_statistical_weight[i][1];
  const double Texp = std::exp(-ksi * Tinv);
  return collision_strength * Texp * inv_omega_1 /
         (A + collision_strength * (inv_omega_2 + Texp * inv_omega_1));
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
double LineCoolingData::get_cooling(
    double temperature, double electron_density,
    const double abundances[LINECOOLINGDATA_NUMELEMENTS]) const {

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
    const double cl2 =
        level_populations[1] *
        _five_level_transition_probability[j][TRANSITION_0_to_1] *
        _five_level_energy_difference[j][TRANSITION_0_to_1];
    const double cl3 =
        level_populations[2] *
        (_five_level_transition_probability[j][TRANSITION_0_to_2] *
             _five_level_energy_difference[j][TRANSITION_0_to_2] +
         _five_level_transition_probability[j][TRANSITION_1_to_2] *
             _five_level_energy_difference[j][TRANSITION_1_to_2]);
    const double cl4 =
        level_populations[3] *
        (_five_level_transition_probability[j][TRANSITION_0_to_3] *
             _five_level_energy_difference[j][TRANSITION_0_to_3] +
         _five_level_transition_probability[j][TRANSITION_1_to_3] *
             _five_level_energy_difference[j][TRANSITION_1_to_3] +
         _five_level_transition_probability[j][TRANSITION_2_to_3] *
             _five_level_energy_difference[j][TRANSITION_2_to_3]);
    const double cl5 =
        level_populations[4] *
        (_five_level_transition_probability[j][TRANSITION_0_to_4] *
             _five_level_energy_difference[j][TRANSITION_0_to_4] +
         _five_level_transition_probability[j][TRANSITION_1_to_4] *
             _five_level_energy_difference[j][TRANSITION_1_to_4] +
         _five_level_transition_probability[j][TRANSITION_2_to_4] *
             _five_level_energy_difference[j][TRANSITION_2_to_4] +
         _five_level_transition_probability[j][TRANSITION_3_to_4] *
             _five_level_energy_difference[j][TRANSITION_3_to_4]);

    cooling += abundances[j] * kb * (cl2 + cl3 + cl4 + cl5);
  }

  /// 2 level atoms

  // offset of two level elements in the abundances array
  const int offset = LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  for (int i = 0; i < LINECOOLINGDATA_NUMTWOLEVELELEMENTS; ++i) {

    const int index = i + offset;
    const LineCoolingDataTwoLevelElement element =
        static_cast< LineCoolingDataTwoLevelElement >(index);
    const double level_population = compute_level_population(
        element, collision_strength_prefactor, temperature, Tinv, logT);
    cooling += abundances[index] * kb * _two_level_energy_difference[i] *
               _two_level_transition_probability[i] * level_population;
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
std::vector< std::vector< double > > LineCoolingData::get_line_strengths(
    double temperature, double electron_density,
    const double abundances[LINECOOLINGDATA_NUMELEMENTS]) const {

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
      LINECOOLINGDATA_NUMELEMENTS);

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
        _five_level_transition_probability[j][TRANSITION_0_to_1] *
        _five_level_energy_difference[j][TRANSITION_0_to_1];
    line_strengths[j][TRANSITION_0_to_2] =
        prefactor * level_populations[2] *
        _five_level_transition_probability[j][TRANSITION_0_to_2] *
        _five_level_energy_difference[j][TRANSITION_0_to_2];
    line_strengths[j][TRANSITION_1_to_2] =
        prefactor * level_populations[2] *
        _five_level_transition_probability[j][TRANSITION_1_to_2] *
        _five_level_energy_difference[j][TRANSITION_1_to_2];
    line_strengths[j][TRANSITION_0_to_3] =
        prefactor * level_populations[3] *
        _five_level_transition_probability[j][TRANSITION_0_to_3] *
        _five_level_energy_difference[j][TRANSITION_0_to_3];
    line_strengths[j][TRANSITION_1_to_3] =
        prefactor * level_populations[3] *
        _five_level_transition_probability[j][TRANSITION_1_to_3] *
        _five_level_energy_difference[j][TRANSITION_1_to_3];
    line_strengths[j][TRANSITION_2_to_3] =
        prefactor * level_populations[3] *
        _five_level_transition_probability[j][TRANSITION_2_to_3] *
        _five_level_energy_difference[j][TRANSITION_2_to_3];
    line_strengths[j][TRANSITION_0_to_4] =
        prefactor * level_populations[4] *
        _five_level_transition_probability[j][TRANSITION_0_to_4] *
        _five_level_energy_difference[j][TRANSITION_0_to_4];
    line_strengths[j][TRANSITION_1_to_4] =
        prefactor * level_populations[4] *
        _five_level_transition_probability[j][TRANSITION_1_to_4] *
        _five_level_energy_difference[j][TRANSITION_1_to_4];
    line_strengths[j][TRANSITION_2_to_4] =
        prefactor * level_populations[4] *
        _five_level_transition_probability[j][TRANSITION_2_to_4] *
        _five_level_energy_difference[j][TRANSITION_2_to_4];
    line_strengths[j][TRANSITION_3_to_4] =
        prefactor * level_populations[4] *
        _five_level_transition_probability[j][TRANSITION_3_to_4] *
        _five_level_energy_difference[j][TRANSITION_3_to_4];
  }

  /// 2 level elements

  // offset of two level elements in the abundances array
  const int offset = LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  for (int i = 0; i < LINECOOLINGDATA_NUMTWOLEVELELEMENTS; ++i) {

    const int index = i + offset;

    line_strengths[index].resize(1);

    const LineCoolingDataTwoLevelElement element =
        static_cast< LineCoolingDataTwoLevelElement >(index);

    const double level_population = compute_level_population(
        element, collision_strength_prefactor, temperature, Tinv, logT);

    line_strengths[index][0] = abundances[index] * kb * level_population *
                               _two_level_energy_difference[i] *
                               _two_level_transition_probability[i];
  }

  return line_strengths;
}
