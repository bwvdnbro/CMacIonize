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
 * @brief Read a given number of values from the given string into the given
 * array.
 *
 * @param line string containing comma separated values.
 * @param array Pointer to a double array of at least the given size which will
 * be filled with values.
 * @param size Number of values to read.
 * @return True if the read was successful, false otherwise.
 */
bool LineCoolingData::read_values(std::string line, double *array,
                                  unsigned int size) {
  int start, end;
  unsigned int index;
  start = 0;
  end = line.find(',', start);
  index = 0;
  while (end > 0 && index < size) {
    array[index] = atof(line.substr(start, end - start).c_str());
    start = end + 1;
    end = line.find(',', start);
    ++index;
  }
  if (end < 0) {
    // try to extract a value from the remainder of the string
    array[index] = atof(line.substr(start, end - start).c_str());
    ++index;
  }
  // if index == size, we know we were able to read size values
  // if not, we failed to fill the array and we return false
  return (index == size);
}

/**
 * @brief Constructor.
 *
 * Reads the data file and stores the necessary quantities in internal arrays.
 */
LineCoolingData::LineCoolingData() {

  // five level elements

  std::ifstream file(LINECOOLINGDATALOCATION);
  std::string line;

  // there are 10 groups of 5 rows
  for (unsigned int i = 0; i < LINECOOLINGDATA_NUMFIVELEVELELEMENTS; ++i) {

    // the first row holds energy level information (in cm^-1)
    std::getline(file, line);
    double enlev[5];
    read_values(line, enlev, 5);
    // convert from cm^-1 to K and compute energy differences for the
    // transitions
    _energy_difference[i][TRANSITION_0_to_1] = enlev[1] * 1.439;
    _energy_difference[i][TRANSITION_0_to_2] = enlev[2] * 1.439;
    _energy_difference[i][TRANSITION_0_to_3] = enlev[3] * 1.439;
    _energy_difference[i][TRANSITION_0_to_4] = enlev[4] * 1.439;
    _energy_difference[i][TRANSITION_1_to_2] = (enlev[2] - enlev[1]) * 1.439;
    _energy_difference[i][TRANSITION_1_to_3] = (enlev[3] - enlev[1]) * 1.439;
    _energy_difference[i][TRANSITION_1_to_4] = (enlev[4] - enlev[1]) * 1.439;
    _energy_difference[i][TRANSITION_2_to_3] = (enlev[3] - enlev[2]) * 1.439;
    _energy_difference[i][TRANSITION_2_to_4] = (enlev[4] - enlev[2]) * 1.439;
    _energy_difference[i][TRANSITION_3_to_4] = (enlev[4] - enlev[3]) * 1.439;

    // the second row contains the transition probabilities for deexcitation
    // between different levels (in s^-1)
    getline(file, line);
    read_values(line, _transition_probability[i], NUMBER_OF_TRANSITIONS);

    // the third row contains the velocity-averaged collision strengths (at
    // 10,000 K)
    getline(file, line);
    read_values(line, _collision_strength[i], NUMBER_OF_TRANSITIONS);

    // the fourth row contains exponents for the temperature variation of the
    // collision strengths
    getline(file, line);
    read_values(line, _collision_strength_exponent[i], NUMBER_OF_TRANSITIONS);

    // the fifth row contains the statistical weights of the different levels
    getline(file, line);
    read_values(line, _inverse_statistical_weight[i], 5);
    // invert the sw values, since they are only used in divisions
    for (unsigned char j = 0; j < 5; ++j) {
      _inverse_statistical_weight[i][j] =
          1. / _inverse_statistical_weight[i][j];
    }
  }

  // we can only change to the correct conversion factor once all energy values
  // have been updated in the old code, so that we can change it in the old code
  // as well...
  // conversion factor from cm^-1 to K (for energy differences)
  //  const double hc_over_k = 100. *
  //  PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK)
  //      *PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED)/
  //      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
  const double hc_over_k = 1.439;
  // data from Galavis, Mendoza & Zeippen (1997), tables 4 and 5
  // ground state: 3P0
  // excited states: 3P1, 3P2, 1D2, 1S0
  const double energy_levels[4] = {48.7, 130.8, 15316.3, 32688.9};
  _inverse_statistical_weight[NII][0] = 1.;
  _inverse_statistical_weight[NII][1] = 1. / 3.;
  _inverse_statistical_weight[NII][2] = 0.2;
  _inverse_statistical_weight[NII][3] = 0.2;
  _inverse_statistical_weight[NII][4] = 1.;
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
  // need to use data from Lennon & Burke (1994)
  //  _collision_strength[NII][TRANSITION_0_to_1] = 0.;
  //  _collision_strength[NII][TRANSITION_0_to_2] = 0.;
  //  _collision_strength[NII][TRANSITION_0_to_3] = 0.;
  //  _collision_strength[NII][TRANSITION_0_to_4] = 0.;
  //  _collision_strength[NII][TRANSITION_1_to_2] = 0.;
  //  _collision_strength[NII][TRANSITION_1_to_3] = 0.;
  //  _collision_strength[NII][TRANSITION_1_to_4] = 0.;
  //  _collision_strength[NII][TRANSITION_2_to_3] = 0.;
  //  _collision_strength[NII][TRANSITION_2_to_4] = 0.;
  //  _collision_strength[NII][TRANSITION_3_to_4] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_0_to_1] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_0_to_2] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_0_to_3] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_0_to_4] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_1_to_2] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_1_to_3] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_1_to_4] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_2_to_3] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_2_to_4] = 0.;
  //  _collision_strength_exponent[NII][TRANSITION_3_to_4] = 0.;

  // two level elements

  const double econst =
      PhysicalConstants::get_physical_constant(
          PHYSICALCONSTANT_RYDBERG_ENERGY) /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
  // Blum & Pradhan (1992), table 5, first energy level
  _two_level_element_data[NIII][TWOLEVELFIELD_ENERGY_DIFFERENCE] =
      0.00159 * econst;
  // Galavis, Mendoza & Zeippen (1998), table 4, 1 to 2 transition
  _two_level_element_data[NIII][TWOLEVELFIELD_TRANSITION_PROBABILITY] =
      4.736e-5;
  // Blum & Pradhan (1992), table 3, value for 10,000 K, 1 to 2 transition
  _two_level_element_data[NIII][TWOLEVELFIELD_COLLISION_STRENGTH] = 1.4454;
  // statistical weights: level 0 is a P_{1/2} level, while level 1 is a P_{3/2}
  _two_level_element_data[NIII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0] =
      0.5;
  _two_level_element_data[NIII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1] =
      0.25;

  // Saraph & Tully (1994), table 2, fine structure splitting energy for Z = 10
  _two_level_element_data[NeII][TWOLEVELFIELD_ENERGY_DIFFERENCE] =
      0.0071 * econst;
  // Kaufman & Sugar (1986), table 7
  _two_level_element_data[NeII][TWOLEVELFIELD_TRANSITION_PROBABILITY] = 8.55e-3;
  // Griffin, Mitnik & Badnell (2001), table 4, value for 10,000 K
  _two_level_element_data[NeII][TWOLEVELFIELD_COLLISION_STRENGTH] = 0.314;
  // statistical weights: level 0 is a P_{3/2} level, while level 1 is a P_{1/2}
  _two_level_element_data[NeII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0] =
      0.25;
  _two_level_element_data[NeII][TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1] =
      0.5;
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

  for (unsigned int j = 0; j < 5; ++j) {
    // find the next row with the largest coefficient
    unsigned int imax = 0;
    double Amax = 0.;
    for (unsigned int i = j; i < 5; ++i) {
      if (std::abs(A[i][j]) > std::abs(Amax)) {
        Amax = A[i][j];
        imax = i;
      }
    }
    // check that the matrix is non-singular
    if (Amax == 0.) {
      return 1;
    }
    // imax now contains the index of the row with the largest coefficient
    // interchange rows if necessary to make sure that the row with the largest
    // coefficient is the next row to eliminate with
    for (unsigned int k = 0; k < 5; ++k) {
      if (imax != j) {
        double save = A[j][k];
        A[j][k] = A[imax][k];
        A[imax][k] = save;
      }
      A[j][k] /= Amax;
    }
    if (imax != j) {
      double save = B[j];
      B[j] = B[imax];
      B[imax] = save;
    }
    B[j] /= Amax;
    // row j now contains a 1 at position j
    // all elements in columns with indices smaller than j are supposed to be
    // zero due to previous eliminations (we do not set them to zero however to
    // save computations)
    if (j < 4) {
      // we are not finished yet: use row j to eliminate all rows below row j
      for (unsigned int i = j + 1; i < 5; ++i) {
        for (unsigned int k = j + 1; k < 5; ++k) {
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
  for (unsigned int i = 0; i < 4; ++i) {
    for (unsigned int j = 0; j < i + 1; ++j) {
      B[3 - i] -= B[4 - j] * A[3 - i][4 - j];
    }
  }
  return 0;
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

  // Boltzmann constant (in J s^-1)
  const double kb =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
  // Planck constant (in J s)
  const double h =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  // electron mass (in kg)
  const double m_e =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_ELECTRON_MASS);

  const double numfac =
      h * h / (std::sqrt(kb) * std::pow(2. * M_PI * m_e, 1.5));
  const double cfac = numfac * electron_density / std::sqrt(temperature);
  const double T4 = temperature * 1.e-4;
  const double Tinv = 1. / temperature;

  double Om[10][10];
  double cs[10][10], cse[10][10];
  for (unsigned int j = 0; j < 10; ++j) {
    for (unsigned int mm = 0; mm < 10; ++mm) {
      cs[j][mm] = _collision_strength[j][mm];
      cse[j][mm] = _collision_strength_exponent[j][mm];
    }
  }
  double T1;
  const double A1 = std::pow(T4, 0.91);
  const double A2 = std::pow(T4, 1.11);
  const double A3 = std::pow(T4, 0.8);
  cs[0][0] = 0.29 * A1;
  cs[0][1] = 0.194 * A1;
  cs[0][2] = 0.113 * A1;
  cs[0][3] = 0.0567 * A1;
  cs[0][4] = 0.269 * A2;
  cs[0][5] = 0.266 * std::pow(T4, 0.715);
  cs[0][6] = 0.109 * A3;
  cs[0][7] = 0.147 * A3;
  cs[0][8] = 0.097 * std::pow(T4, 0.69);
  cs[0][9] = 0.071 * A2;
  T1 = 0.266 * A2;
  const double T2 = 0.0324 * A2;
  cs[2][0] = 0.0987 * A2;
  cs[2][1] = 0.0292 * T4;
  cs[2][2] = 0.55556 * T1;
  cs[2][3] = 0.55556 * T2;
  cs[2][4] = 0.0265 * std::pow(T4, 1.24);
  cs[2][5] = 0.333 * T1;
  cs[2][6] = 0.333 * T2;
  cs[2][7] = T1 / 9.;
  cs[2][8] = T2 / 9.;
  cs[2][9] = 0.105 * std::pow(T4, 0.52);

  double cooling = 0.;
  double alev[5][5], lev[5];
  for (unsigned int j = 0; j < 10; ++j) {
    for (unsigned int mm = 0; mm < 5; ++mm) {
      alev[0][mm] = 1.;
      lev[mm] = 0.;
    }
    lev[0] = 1.;
    for (unsigned int mm = 0; mm < 10; ++mm) {
      Om[j][mm] = cs[j][mm] * std::pow(T4, cse[j][mm]);
    }

    alev[1][0] = cfac * Om[j][TRANSITION_0_to_1] *
                 _inverse_statistical_weight[j][0] *
                 std::exp(-_energy_difference[j][TRANSITION_0_to_1] * Tinv);
    alev[1][1] =
        -(_transition_probability[j][TRANSITION_0_to_1] +
          cfac * _inverse_statistical_weight[j][1] *
              (Om[j][TRANSITION_0_to_1] +
               Om[j][TRANSITION_1_to_2] *
                   std::exp(-_energy_difference[j][TRANSITION_1_to_2] * Tinv) +
               Om[j][TRANSITION_1_to_3] *
                   std::exp(-_energy_difference[j][TRANSITION_1_to_3] * Tinv) +
               Om[j][TRANSITION_1_to_4] *
                   std::exp(-_energy_difference[j][TRANSITION_1_to_4] * Tinv)));
    alev[1][2] =
        _transition_probability[j][TRANSITION_1_to_2] +
        cfac * _inverse_statistical_weight[j][2] * Om[j][TRANSITION_1_to_2];
    alev[1][3] =
        _transition_probability[j][TRANSITION_1_to_3] +
        cfac * _inverse_statistical_weight[j][3] * Om[j][TRANSITION_1_to_3];
    alev[1][4] =
        _transition_probability[j][TRANSITION_1_to_4] +
        cfac * _inverse_statistical_weight[j][4] * Om[j][TRANSITION_1_to_4];

    alev[2][0] = cfac * Om[j][TRANSITION_0_to_2] *
                 std::exp(-_energy_difference[j][TRANSITION_0_to_2] * Tinv) *
                 _inverse_statistical_weight[j][0];
    alev[2][1] = cfac * Om[j][TRANSITION_1_to_2] *
                 std::exp(-_energy_difference[j][TRANSITION_1_to_2] * Tinv) *
                 _inverse_statistical_weight[j][1];
    alev[2][2] =
        -(_transition_probability[j][TRANSITION_0_to_2] +
          _transition_probability[j][TRANSITION_1_to_2] +
          cfac * _inverse_statistical_weight[j][2] *
              (Om[j][TRANSITION_0_to_2] + Om[j][TRANSITION_1_to_2] +
               Om[j][TRANSITION_2_to_3] *
                   std::exp(-_energy_difference[j][TRANSITION_2_to_3] * Tinv) +
               Om[j][TRANSITION_2_to_4] *
                   std::exp(-_energy_difference[j][TRANSITION_2_to_4] * Tinv)));
    alev[2][3] =
        _transition_probability[j][TRANSITION_2_to_3] +
        cfac * Om[j][TRANSITION_2_to_3] * _inverse_statistical_weight[j][3];
    alev[2][4] =
        _transition_probability[j][TRANSITION_2_to_4] +
        cfac * Om[j][TRANSITION_2_to_4] * _inverse_statistical_weight[j][4];

    alev[3][0] = cfac * Om[j][TRANSITION_0_to_3] *
                 std::exp(-_energy_difference[j][TRANSITION_0_to_3] * Tinv) *
                 _inverse_statistical_weight[j][0];
    alev[3][1] = cfac * Om[j][TRANSITION_1_to_3] *
                 std::exp(-_energy_difference[j][TRANSITION_1_to_3] * Tinv) *
                 _inverse_statistical_weight[j][1];
    alev[3][2] = cfac * Om[j][TRANSITION_2_to_3] *
                 std::exp(-_energy_difference[j][TRANSITION_2_to_3] * Tinv) *
                 _inverse_statistical_weight[j][2];
    alev[3][3] =
        -(_transition_probability[j][TRANSITION_0_to_3] +
          _transition_probability[j][TRANSITION_1_to_3] +
          _transition_probability[j][TRANSITION_2_to_3] +
          cfac * _inverse_statistical_weight[j][3] *
              (Om[j][TRANSITION_0_to_3] + Om[j][TRANSITION_1_to_3] +
               Om[j][TRANSITION_2_to_3] +
               Om[j][TRANSITION_3_to_4] *
                   std::exp(-_energy_difference[j][TRANSITION_3_to_4] * Tinv)));
    alev[3][4] =
        _transition_probability[j][TRANSITION_3_to_4] +
        cfac * Om[j][TRANSITION_3_to_4] * _inverse_statistical_weight[j][4];

    alev[4][0] = cfac * Om[j][TRANSITION_0_to_4] *
                 std::exp(-_energy_difference[j][TRANSITION_0_to_4] * Tinv) *
                 _inverse_statistical_weight[j][0];
    alev[4][1] = cfac * Om[j][TRANSITION_1_to_4] *
                 std::exp(-_energy_difference[j][TRANSITION_1_to_4] * Tinv) *
                 _inverse_statistical_weight[j][1];
    alev[4][2] = cfac * Om[j][TRANSITION_2_to_4] *
                 std::exp(-_energy_difference[j][TRANSITION_2_to_4] * Tinv) *
                 _inverse_statistical_weight[j][2];
    alev[4][3] = cfac * Om[j][TRANSITION_3_to_4] *
                 std::exp(-_energy_difference[j][TRANSITION_3_to_4] * Tinv) *
                 _inverse_statistical_weight[j][3];
    alev[4][4] = -(_transition_probability[j][TRANSITION_0_to_4] +
                   _transition_probability[j][TRANSITION_1_to_4] +
                   _transition_probability[j][TRANSITION_2_to_4] +
                   _transition_probability[j][TRANSITION_3_to_4] +
                   cfac * _inverse_statistical_weight[j][4] *
                       (Om[j][TRANSITION_0_to_4] + Om[j][TRANSITION_1_to_4] +
                        Om[j][TRANSITION_2_to_4] + Om[j][TRANSITION_3_to_4]));

    // find level populations
    int status = simq(alev, lev);
    if (status != 0) {
      // something went wrong
      cmac_warning("Singular matrix given to simq!");
      cmac_warning("Temperature: %g", temperature);
      cmac_warning("Electron density: %g", electron_density);
      cmac_error("We better stop!");
    }

    const double cl2 = abundances[j] * kb * lev[1] *
                       _transition_probability[j][TRANSITION_0_to_1] *
                       _energy_difference[j][TRANSITION_0_to_1];
    const double cl3 = abundances[j] * kb * lev[2] *
                       (_transition_probability[j][TRANSITION_0_to_2] *
                            _energy_difference[j][TRANSITION_0_to_2] +
                        _transition_probability[j][TRANSITION_1_to_2] *
                            _energy_difference[j][TRANSITION_1_to_2]);
    const double cl4 = abundances[j] * kb * lev[3] *
                       (_transition_probability[j][TRANSITION_0_to_3] *
                            _energy_difference[j][TRANSITION_0_to_3] +
                        _transition_probability[j][TRANSITION_1_to_3] *
                            _energy_difference[j][TRANSITION_1_to_3] +
                        _transition_probability[j][TRANSITION_2_to_3] *
                            _energy_difference[j][TRANSITION_2_to_3]);
    const double cl5 = abundances[j] * kb * lev[4] *
                       (_transition_probability[j][TRANSITION_0_to_4] *
                            _energy_difference[j][TRANSITION_0_to_4] +
                        _transition_probability[j][TRANSITION_1_to_4] *
                            _energy_difference[j][TRANSITION_1_to_4] +
                        _transition_probability[j][TRANSITION_2_to_4] *
                            _energy_difference[j][TRANSITION_2_to_4] +
                        _transition_probability[j][TRANSITION_3_to_4] *
                            _energy_difference[j][TRANSITION_3_to_4]);

    const double coolj = cl2 + cl3 + cl4 + cl5;
    cooling += coolj;
  }

  // 2 level atoms

  // offset of two level elements in the abundances array
  const int offset = LINECOOLINGDATA_NUMFIVELEVELELEMENTS;
  for (int i = 0; i < LINECOOLINGDATA_NUMTWOLEVELELEMENTS; ++i) {
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
    const double n = abundances[offset + i];
    const double Texp = std::exp(-ksi * Tinv);
    const double Celement =
        n * kb * cfac * ksi * Gamma * Texp * A * inv_omega_1 /
        (A + cfac * Gamma * (inv_omega_2 + Texp * inv_omega_1));
    cooling += Celement;
  }

  return cooling;
}

/**
 * @brief Calculate the strength of a number of emission lines for the given
 * temperature, electron density and ion abundances.
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

  // Boltzmann constant (in J s^-1)
  const double kb =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
  // Planck constant (in J s)
  const double h =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  // electron mass (in kg)
  const double m_e =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_ELECTRON_MASS);

  const double numfac =
      h * h / (std::sqrt(kb) * std::pow(2. * M_PI * m_e, 1.5));
  const double cfac = numfac * electron_density / std::sqrt(temperature);
  const double T4 = temperature * 1.e-4;
  const double Tinv = 1. / temperature;

  double Om[10][10];
  double cs[10][10], cse[10][10];
  for (unsigned int j = 0; j < 10; ++j) {
    for (unsigned int mm = 0; mm < 10; ++mm) {
      cs[j][mm] = _collision_strength[j][mm];
      cse[j][mm] = _collision_strength_exponent[j][mm];
    }
  }
  const double A1 = std::pow(T4, 0.91);
  const double A2 = std::pow(T4, 1.11);
  const double A3 = std::pow(T4, 0.8);
  cs[0][0] = 0.29 * A1;
  cs[0][1] = 0.194 * A1;
  cs[0][2] = 0.113 * A1;
  cs[0][3] = 0.0567 * A1;
  cs[0][4] = 0.269 * A2;
  cs[0][5] = 0.266 * std::pow(T4, 0.715);
  cs[0][6] = 0.109 * A3;
  cs[0][7] = 0.147 * A3;
  cs[0][8] = 0.097 * std::pow(T4, 0.69);
  cs[0][9] = 0.071 * A2;
  double T1 = 0.266 * A2;
  const double T2 = 0.0324 * A2;
  cs[2][0] = 0.0987 * A2;
  cs[2][1] = 0.0292 * T4;
  cs[2][2] = 0.55556 * T1;
  cs[2][3] = 0.55556 * T2;
  cs[2][4] = 0.0265 * std::pow(T4, 1.24);
  cs[2][5] = 0.333 * T1;
  cs[2][6] = 0.333 * T2;
  cs[2][7] = T1 / 9.;
  cs[2][8] = T2 / 9.;
  cs[2][9] = 0.105 * std::pow(T4, 0.52);

  double alev[5][5], lev[5];
  for (unsigned int j = 0; j < 10; ++j) {
    for (unsigned int mm = 0; mm < 5; ++mm) {
      alev[0][mm] = 1.;
      lev[mm] = 0.;
    }
    lev[0] = 1.;
    for (unsigned int mm = 0; mm < 10; ++mm) {
      Om[j][mm] = cs[j][mm] * std::pow(T4, cse[j][mm]);
    }

    alev[1][0] = cfac * Om[j][0] * _inverse_statistical_weight[j][0] *
                 std::exp(-_energy_difference[j][0] * Tinv);
    alev[1][1] = -(_transition_probability[j][0] +
                   cfac * _inverse_statistical_weight[j][1] *
                       (Om[j][0] +
                        Om[j][4] * std::exp(-_energy_difference[j][4] * Tinv) +
                        Om[j][5] * std::exp(-_energy_difference[j][5] * Tinv) +
                        Om[j][6] * std::exp(-_energy_difference[j][6] * Tinv)));
    alev[1][2] = _transition_probability[j][4] +
                 cfac * _inverse_statistical_weight[j][2] * Om[j][4];
    alev[1][3] = _transition_probability[j][5] +
                 cfac * _inverse_statistical_weight[j][3] * Om[j][5];
    alev[1][4] = _transition_probability[j][6] +
                 cfac * _inverse_statistical_weight[j][4] * Om[j][6];

    alev[2][0] = cfac * Om[j][1] * std::exp(-_energy_difference[j][1] * Tinv) *
                 _inverse_statistical_weight[j][0];
    alev[2][1] = cfac * Om[j][4] * std::exp(-_energy_difference[j][4] * Tinv) *
                 _inverse_statistical_weight[j][1];
    alev[2][2] =
        -(_transition_probability[j][1] + _transition_probability[j][4] +
          cfac * _inverse_statistical_weight[j][2] *
              (Om[j][1] + Om[j][4] +
               Om[j][7] * std::exp(-_energy_difference[j][7] * Tinv) +
               Om[j][8] * std::exp(-_energy_difference[j][8] * Tinv)));
    alev[2][3] = _transition_probability[j][7] +
                 cfac * Om[j][7] * _inverse_statistical_weight[j][3];
    alev[2][4] = _transition_probability[j][8] +
                 cfac * Om[j][8] * _inverse_statistical_weight[j][4];

    alev[3][0] = cfac * Om[j][2] * std::exp(-_energy_difference[j][2] * Tinv) *
                 _inverse_statistical_weight[j][0];
    alev[3][1] = cfac * Om[j][5] * std::exp(-_energy_difference[j][5] * Tinv) *
                 _inverse_statistical_weight[j][1];
    alev[3][2] = cfac * Om[j][7] * std::exp(-_energy_difference[j][7] * Tinv) *
                 _inverse_statistical_weight[j][2];
    alev[3][3] =
        -(_transition_probability[j][2] + _transition_probability[j][5] +
          _transition_probability[j][7] +
          cfac * _inverse_statistical_weight[j][3] *
              (Om[j][2] + Om[j][5] + Om[j][7] +
               Om[j][9] * std::exp(-_energy_difference[j][9] * Tinv)));
    alev[3][4] = _transition_probability[j][9] +
                 cfac * Om[j][9] * _inverse_statistical_weight[j][4];

    alev[4][0] = cfac * Om[j][3] * std::exp(-_energy_difference[j][3] * Tinv) *
                 _inverse_statistical_weight[j][0];
    alev[4][1] = cfac * Om[j][6] * std::exp(-_energy_difference[j][6] * Tinv) *
                 _inverse_statistical_weight[j][1];
    alev[4][2] = cfac * Om[j][8] * std::exp(-_energy_difference[j][8] * Tinv) *
                 _inverse_statistical_weight[j][2];
    alev[4][3] = cfac * Om[j][9] * std::exp(-_energy_difference[j][9] * Tinv) *
                 _inverse_statistical_weight[j][3];
    alev[4][4] =
        -(_transition_probability[j][3] + _transition_probability[j][6] +
          _transition_probability[j][8] + _transition_probability[j][9] +
          cfac * _inverse_statistical_weight[j][4] *
              (Om[j][3] + Om[j][6] + Om[j][8] + Om[j][9]));

    // find level populations
    const int status = simq(alev, lev);
    if (status != 0) {
      // something went wrong
      cmac_warning("Singular matrix given to simq!");
      cmac_warning("Temperature: %g", temperature);
      cmac_warning("Electron density: %g", electron_density);
      cmac_error("We better stop!");
    }

    const double cl2 = abundances[j] * kb * lev[1] *
                       _transition_probability[j][0] * _energy_difference[j][0];
    const double cl3 =
        abundances[j] * kb * lev[2] *
        (_transition_probability[j][1] * _energy_difference[j][1] +
         _transition_probability[j][4] * _energy_difference[j][4]);
    const double cl4 =
        abundances[j] * kb * lev[3] *
        (_transition_probability[j][2] * _energy_difference[j][2] +
         _transition_probability[j][5] * _energy_difference[j][5] +
         _transition_probability[j][7] * _energy_difference[j][7]);
    const double cl5 =
        abundances[j] * kb * lev[4] *
        (_transition_probability[j][3] * _energy_difference[j][3] +
         _transition_probability[j][6] * _energy_difference[j][6] +
         _transition_probability[j][8] * _energy_difference[j][8] +
         _transition_probability[j][9] * _energy_difference[j][9]);

    if (j == 1) {
      c5755 = abundances[j] * kb * lev[4] * _transition_probability[j][9] *
              _energy_difference[j][9];
      c6584 = abundances[j] * kb * lev[3] * _transition_probability[j][7] *
              _energy_difference[j][7];
      cnii122 = abundances[j] * kb * lev[2] * _transition_probability[j][4] *
                _energy_difference[j][4];
    }
    if (j == 2) {
      c6300_6363 = abundances[j] * kb * lev[3] *
                   (_transition_probability[j][2] * _energy_difference[j][2] +
                    _transition_probability[j][5] * _energy_difference[j][5]);
    }
    if (j == 3) {
      c3729 = cl2;
      c3727 = cl2 + cl3;
      coii7325 =
          abundances[j] * kb *
          (lev[4] * (_transition_probability[j][6] * _energy_difference[j][6] +
                     _transition_probability[j][8] * _energy_difference[j][8]) +
           lev[3] * (_transition_probability[j][5] * _energy_difference[j][5] +
                     _transition_probability[j][7] * _energy_difference[j][7]));
    }
    if (j == 4) {
      c4363 = abundances[j] * kb * lev[4] * _transition_probability[j][9] *
              _energy_difference[j][9];
      c5007 = abundances[j] * kb * lev[3] * _transition_probability[j][7] *
              _energy_difference[j][7];
      c52mu = abundances[j] * kb * lev[2] * _transition_probability[j][4] *
              _energy_difference[j][4];
      c88mu = abundances[j] * kb * lev[1] * _transition_probability[j][0] *
              _energy_difference[j][0];
    }
    if (j == 5) {
      c3869 = abundances[j] * kb * lev[3] * _transition_probability[j][2] *
              _energy_difference[j][2];
      cneiii15 = cl2;
    }
    if (j == 6) {
      c4072 =
          abundances[j] * kb *
          (lev[3] * _transition_probability[j][2] * _energy_difference[j][2] +
           lev[4] * _transition_probability[j][3] * _energy_difference[j][3]);
      c6717 = cl3;
      c6725 = cl2 + cl3;
    }
    if (j == 7) {
      c9405 = abundances[j] * kb * lev[3] *
              (_transition_probability[j][5] * _energy_difference[j][5] +
               _transition_probability[j][7] * _energy_difference[j][7]);
      c6312 = abundances[j] * kb * lev[4] * _transition_probability[j][9] *
              _energy_difference[j][9];
      c33mu = abundances[j] * kb * lev[1] * _transition_probability[j][0] *
              _energy_difference[j][0];
      c19mu = abundances[j] * kb * lev[2] * _transition_probability[j][4] *
              _energy_difference[j][4];
    }
    if (j == 8) {
      cii2325 = cl3 + cl4 + cl5;
    }
    if (j == 9) {
      ciii1908 = cl2 + cl3 + cl4;
    }
  }

  // 2 level atoms

  const double econst = PhysicalConstants::get_physical_constant(
                            PHYSICALCONSTANT_RYDBERG_ENERGY) /
                        kb;
  // Blum & Pradhan (1992), table 5, first energy level
  const double EnNIII = 0.00159 * econst;
  // Galavis, Mendoza & Zeippen (1998), table 4, 1 to 2 transition
  const double EaNIII = 4.736e-5;
  // Blum & Pradhan (1992), table 3, value for 10,000 K, 1 to 2 transition
  const double OmNIII = 1.4454;

  // Saraph & Tully (1994), table 2, fine structure splitting energy for Z = 10
  const double EnNeII = 0.0071 * econst;
  // Kaufman & Sugar (1986), table 7
  const double EaNeII = 8.55e-3;
  // Griffin, Mitnik & Badnell (2001), table 4, value for 10,000 K
  const double OmNeII = 0.314;

  double sw1 = 2.;
  double sw2 = 4.;
  T1 = std::exp(-EnNIII * Tinv);
  cniii57 = abundances[10] * kb * cfac * EnNIII * OmNIII * T1 * EaNIII /
            (sw1 * (EaNIII + cfac * OmNIII * (1. / sw2 + T1 / sw1)));
  sw1 = 4.;
  sw2 = 2.;
  T1 = std::exp(-EnNeII * Tinv);
  cneii12 = abundances[11] * kb * cfac * OmNeII * EnNeII * T1 * EaNeII /
            (sw1 * (EaNeII + cfac * OmNeII * (1. / sw2 + T1 / sw1)));
}
