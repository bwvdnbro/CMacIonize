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
#include "UnitConverter.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

/**
 * @brief Read a given number of values from the given string into the given
 * array
 *
 * @param line string containing comma separated values.
 * @param array Pointer to a double array of at least the given size which will
 * be filled with values.
 * @param size Number of values to read.
 * @return True if the read was successful, false otherwise.
 */
bool LineCoolingData::read_values(string line, double *array,
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
 * @brief Constructor
 *
 * Reads the data file and stores the necessary quantities in internal arrays.
 */
LineCoolingData::LineCoolingData() {
  double enlev[5];

  ifstream file(LINECOOLINGDATALOCATION);
  string line;

  // there are 10 groups of 5 rows
  for (unsigned int i = 0; i < LINECOOLINGDATA_NUMELEMENTS; ++i) {
    // the first row holds energy level information
    getline(file, line);
    read_values(line, enlev, 5);
    // these values need to be converted before we can store them
    unsigned int l = 0;
    for (unsigned int j = 0; j < 5; ++j) {
      for (unsigned int k = j + 1; k < 5; ++k) {
        // the factor converts from cm^-1 to K
        _en[i][l] = (enlev[k] - enlev[j]) * 1.439;
        ++l;
      }
    }
    // the second row contains the Einstein A values
    getline(file, line);
    read_values(line, _ea[i], 10);
    // the third row contains the omega values
    getline(file, line);
    read_values(line, _cs[i], 10);
    // the fourth row contains the exponential omega values
    getline(file, line);
    read_values(line, _cse[i], 10);
    // the fifth row contains the sw values
    getline(file, line);
    read_values(line, _sw[i], 5);
  }
}

/**
 * @brief Get the cs value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return cs value
 */
double LineCoolingData::get_cs(unsigned int element, unsigned int level) {
  return _cs[element][level];
}

/**
 * @brief Get the cse value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return cse value
 */
double LineCoolingData::get_cse(unsigned int element, unsigned int level) {
  return _cse[element][level];
}

/**
 * @brief Get the ea value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return ea value
 */
double LineCoolingData::get_ea(unsigned int element, unsigned int level) {
  return _ea[element][level];
}

/**
 * @brief Get the en value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return en value
 */
double LineCoolingData::get_en(unsigned int element, unsigned int level) {
  return _en[element][level];
}

/**
 * @brief Get the sw value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return sw value
 */
double LineCoolingData::get_sw(unsigned int element, unsigned int level) {
  return _sw[element][level];
}

/**
 * @brief Solve a system of 5 coupled linear equations.
 *
 * We assume a system of equations of the form
 * \f[
 *   A.X = B,
 * \f]
 * with
 * \f[
 * A = \begin{pmatrix}
 * A00 & A01 & A02 & A03 & A04 \\
 * A10 & A11 & A12 & A13 & A14 \\
 * A20 & A21 & A22 & A23 & A24 \\
 * A30 & A31 & A32 & A33 & A34 \\
 * A40 & A41 & A42 & A43 & A44
 * \end{pmatrix},
 * \f]
 * and
 * \f[
 * B = \begin{pmatrix}
 * B0 \\
 * B1 \\
 * B2 \\
 * B3 \\
 * B4
 * \end{pmatrix}.
 * \f]
 * We want to solve for the elements of the column matrix \f$X\f$.
 *
 * To do this, we first reduce the combined matrix \f$AB\f$ (the matrix \f$A\f$
 * with the matrix \f$B\f$ added as an extra column) to an upper triangular
 * matrix using Gaussian elimination. This automatically gives us the solution
 * for the last element of \f$X\f$. We then use this element to recursively
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
 * Note that both matrix \f$A\f$ and matrix \f$B\f$ are modified in place. When
 * the method returns, \f$B\f$ contains the elements of the matrix \f$X\f$.
 *
 * @param A Elements of the matrix \f$A\f$.
 * @param B Elements of the matrix \f$B\f$, and elements of the solution on
 * exit.
 */
void LineCoolingData::simq(double A[5][5], double B[5]) {

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
      error("Singular matrix given to simq!");
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
  return;
}

/**
 * @brief Get the radiative energy losses due to line cooling at the given
 * temperature, electron density and coolant abundances.
 *
 * @param temperature Temperature (in K).
 * @param electron_density Electron density (in m^-3).
 * @param abundances Abdunances of coolants.
 * @return Radiative cooling per hydrogen atom (in kg m^2s^-3).
 */
double LineCoolingData::get_cooling(double temperature, double electron_density,
                                    double *abundances) {
  double EnNIII = 251.;
  double EaNIII = 4.77e-5;
  double OmNIII = 0.701;
  double EnNeII = 1125.;
  double EaNeII = 8.55e-3;
  double OmNeII = 0.368;
  double kb = 1.38e-16;

  electron_density = UnitConverter< QUANTITY_NUMBER_DENSITY >::to_unit(
      electron_density, "cm^-3");

  double cfac = 8.63e-6 * electron_density / std::sqrt(temperature);
  double T4 = temperature * 1.e-4;

  double Om[10][10];

  double cooling = 0.;
  double alev[5][5], lev[5];
  for (unsigned int j = 0; j < 10; ++j) {
    for (unsigned int mm = 0; mm < 5; ++mm) {
      alev[0][mm] = 1.;
      lev[mm] = 0.;
    }
    lev[0] = 1.;
    for (unsigned int mm = 0; mm < 10; ++mm) {
      Om[j][mm] = _cs[j][mm] * std::pow(T4, _cse[j][mm]);
    }

    alev[1][0] =
        cfac * Om[j][0] / (_sw[j][0] * std::exp(_en[j][0] / temperature));
    alev[1][1] =
        -(_ea[j][0] +
          (cfac / _sw[j][1]) *
              (Om[j][0] + Om[j][4] * std::exp(-_en[j][4] / temperature) +
               Om[j][5] * std::exp(-_en[j][5] / temperature) +
               Om[j][6] * std::exp(-_en[j][6] / temperature)));
    alev[1][2] = _ea[j][4] + (cfac / _sw[j][2]) * Om[j][4];
    alev[1][3] = _ea[j][5] + (cfac / _sw[j][3]) * Om[j][5];
    alev[1][4] = _ea[j][6] + (cfac / _sw[j][4]) * Om[j][6];

    alev[2][0] =
        cfac * Om[j][1] * std::exp(-_en[j][1] / temperature) / _sw[j][0];
    alev[2][1] =
        cfac * Om[j][4] * std::exp(-_en[j][4] / temperature) / _sw[j][1];
    alev[2][2] =
        -(_ea[j][1] + _ea[j][4] +
          (cfac / _sw[j][2]) * (Om[j][1] + Om[j][4] +
                                Om[j][7] * std::exp(-_en[j][7] / temperature) +
                                Om[j][8] * std::exp(-_en[j][8] / temperature)));
    alev[2][3] = _ea[j][7] + cfac * Om[j][7] / _sw[j][3];
    alev[2][4] = _ea[j][8] + cfac * Om[j][8] / _sw[j][4];

    alev[3][0] =
        cfac * Om[j][2] * std::exp(-_en[j][2] / temperature) / _sw[j][0];
    alev[3][1] =
        cfac * Om[j][5] * std::exp(-_en[j][5] / temperature) / _sw[j][1];
    alev[3][2] =
        cfac * Om[j][7] * std::exp(-_en[j][7] / temperature) / _sw[j][2];
    alev[3][3] =
        -(_ea[j][2] + _ea[j][5] + _ea[j][7] +
          (cfac / _sw[j][3]) * (Om[j][2] + Om[j][5] + Om[j][7] +
                                Om[j][9] * std::exp(-_en[j][9] / temperature)));
    alev[3][4] = _ea[j][9] + cfac * Om[j][9] / _sw[j][4];

    alev[4][0] =
        cfac * Om[j][3] * std::exp(-_en[j][3] / temperature) / _sw[j][0];
    alev[4][1] =
        cfac * Om[j][6] * std::exp(-_en[j][6] / temperature) / _sw[j][1];
    alev[4][2] =
        cfac * Om[j][8] * std::exp(-_en[j][8] / temperature) / _sw[j][2];
    alev[4][3] =
        cfac * Om[j][9] * std::exp(-_en[j][9] / temperature) / _sw[j][3];
    alev[4][4] =
        -(_ea[j][3] + _ea[j][6] + _ea[j][8] + _ea[j][9] +
          (cfac / _sw[j][4]) * (Om[j][3] + Om[j][6] + Om[j][8] + Om[j][9]));

    // find level populations
    simq(alev, lev);

    double cl2 = abundances[j] * kb * lev[1] * _ea[j][0] * _en[j][0];
    double cl3 = abundances[j] * kb * lev[2] *
                 (_ea[j][1] * _en[j][1] + _ea[j][4] * _en[j][4]);
    double cl4 =
        abundances[j] * kb * lev[3] *
        (_ea[j][2] * _en[j][2] + _ea[j][5] * _en[j][5] + _ea[j][7] * _en[j][7]);
    double cl5 = abundances[j] * kb * lev[4] *
                 (_ea[j][3] * _en[j][3] + _ea[j][6] * _en[j][6] +
                  _ea[j][8] * _en[j][8] + _ea[j][9] * _en[j][9]);

    double coolj = cl2 + cl3 + cl4 + cl5;
    cooling += coolj;
  }

  // 2 level atoms
  double sw1 = 2.;
  double sw2 = 4.;
  double T1 = std::exp(-EnNIII / temperature);
  double CNIII =
      abundances[10] * kb * cfac * EnNIII * OmNIII * std::pow(T4, 0.136) * T1 *
      EaNIII /
      (sw1 *
       (EaNIII + cfac * OmNIII * std::pow(T4, 0.136) * (1. / sw2 + T1 / sw1)));
  sw1 = 4.;
  sw2 = 2.;
  T1 = std::exp(-EnNeII / temperature);
  double CNeII = abundances[11] * kb * cfac * OmNeII * EnNeII * T1 * EaNeII /
                 (sw1 * (EaNeII + cfac * OmNeII * (1. / sw2 + T1 / sw1)));
  cooling += CNIII + CNeII;

  // cooling now is in erg/s/hydrogen atom
  // convert to J/s (kg m^2s^-3)
  cooling = UnitConverter< QUANTITY_ENERGY_RATE >::to_SI(cooling, "erg s^-1");

  return cooling;
}

/**
 * @brief Calculate the strength of a number of emission lines for the given
 * temperature, electron density and ion abundances.
 *
 * @param temperature Temperature (in K).
 * @param electron_density Electron density (in m^-3).
 * @param abundances Ion abundances.
 * @param c6300 Variable to store the oxygen 6300 angstrom emission line
 * strength in.
 * @param c9405 Variable to store the sulfur 9405 angstrom emission line
 * strength in.
 * @param c6312 Variable to store the sulfur 6312 angstrom emission line
 * strength in.
 * @param c33mu Variable that is not used in the output.
 * @param c19mu Variable to store the sulfur 18.7 micrometre emission line
 * strength in.
 * @param c3729 Variable to store the oxygen 3729 angstrom emission line
 * strength in.
 * @param c3727 Variable to store the oxygen 3727 angstrom emission line
 * strength in.
 * @param c7330 Variable to store the oxygen 7330 angstrom emission line
 * strength in.
 * @param c4363 Variable to store the oxygen 4363 angstrom emission line
 * strength in.
 * @param c5007 Variable to store the oxygen 5007 angstrom emission line
 * strength in.
 * @param c52mu Variable to store the oxygen 52 micrometre emission line
 * strength in.
 * @param c88mu Variable to store the oxygen 88 micrometre emission line
 * strength in.
 * @param c5755 Variable to store the nitrogen 5755 angstrom emission line
 * strength in.
 * @param c6584 Variable to store the nitrogen 6584 angstrom emission line
 * strength in.
 * @param c4072 Variable to store the sulfur 4072 angstrom emission line
 * strength in.
 * @param c6717 Variable to store the sulfur 6717 angstrom emission line
 * strength in.
 * @param c6725 Variable to store the sulfur 6725 angstrom emission line
 * strength in.
 * @param c3869 Variable to store the neon 3869 angstrom emission line strength
 * in.
 * @param cniii57 Variable to store the nitrogen 57.3 micrometre emission line
 * strength in.
 * @param cneii12 Variable to store the neon 12.8 micrometre emission line
 * strength in.
 * @param cneiii15 Variable to store the neon 15.5 micrometre emission line
 * strength in.
 * @param cnii122 Variable to store the nitrogen 122 micrometre emission line
 * strength in.
 * @param cii2325 Variable to store the carbon 2325 angstrom emission line
 * strength in.
 * @param ciii1908 Variable to store the carbon 1907 + 1909 angstrom emission
 * line strength in.
 * @param coii7325 Variable to store the oxygen 7320 + 7330 angstrom emission
 * line strength in.
 * @param csiv10 Variable to store the sulfur 10 micrometre (?) emission line
 * strength in.
 */
void LineCoolingData::linestr(
    double temperature, double electron_density, double *abundances,
    double &c6300, double &c9405, double &c6312, double &c33mu, double &c19mu,
    double &c3729, double &c3727, double &c7330, double &c4363, double &c5007,
    double &c52mu, double &c88mu, double &c5755, double &c6584, double &c4072,
    double &c6717, double &c6725, double &c3869, double &cniii57,
    double &cneii12, double &cneiii15, double &cnii122, double &cii2325,
    double &ciii1908, double &coii7325, double &csiv10) {
  double EnNIII = 251.;
  double EaNIII = 4.77e-5;
  double OmNIII = 1.45;
  double EnNeII = 1125.;
  double EaNeII = 8.55e-3;
  double OmNeII = 0.303;
  double kb = 1.38e-16;

  electron_density = UnitConverter< QUANTITY_NUMBER_DENSITY >::to_unit(
      electron_density, "cm^-3");

  double cfac = 8.63e-6 * electron_density / std::sqrt(temperature);
  double T4 = temperature * 1.e-4;

  double Om[10][10];
  double cs[10][10], cse[10][10];
  for (unsigned int j = 0; j < 10; ++j) {
    for (unsigned int mm = 0; mm < 10; ++mm) {
      cs[j][mm] = _cs[j][mm];
      cse[j][mm] = _cse[j][mm];
    }
  }
  double A1 = std::pow(T4, 0.91);
  double A2 = std::pow(T4, 1.11);
  double A3 = std::pow(T4, 0.8);
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
  double T2 = 0.0324 * A2;
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

    alev[1][0] =
        cfac * Om[j][0] / (_sw[j][0] * std::exp(_en[j][0] / temperature));
    alev[1][1] =
        -(_ea[j][0] +
          (cfac / _sw[j][1]) *
              (Om[j][0] + Om[j][4] * std::exp(-_en[j][4] / temperature) +
               Om[j][5] * std::exp(-_en[j][5] / temperature) +
               Om[j][6] * std::exp(-_en[j][6] / temperature)));
    alev[1][2] = _ea[j][4] + (cfac / _sw[j][2]) * Om[j][4];
    alev[1][3] = _ea[j][5] + (cfac / _sw[j][3]) * Om[j][5];
    alev[1][4] = _ea[j][6] + (cfac / _sw[j][4]) * Om[j][6];

    alev[2][0] =
        cfac * Om[j][1] * std::exp(-_en[j][1] / temperature) / _sw[j][0];
    alev[2][1] =
        cfac * Om[j][4] * std::exp(-_en[j][4] / temperature) / _sw[j][1];
    alev[2][2] =
        -(_ea[j][1] + _ea[j][4] +
          (cfac / _sw[j][2]) * (Om[j][1] + Om[j][4] +
                                Om[j][7] * std::exp(-_en[j][7] / temperature) +
                                Om[j][8] * std::exp(-_en[j][8] / temperature)));
    alev[2][3] = _ea[j][7] + cfac * Om[j][7] / _sw[j][3];
    alev[2][4] = _ea[j][8] + cfac * Om[j][8] / _sw[j][4];

    alev[3][0] =
        cfac * Om[j][2] * std::exp(-_en[j][2] / temperature) / _sw[j][0];
    alev[3][1] =
        cfac * Om[j][5] * std::exp(-_en[j][5] / temperature) / _sw[j][1];
    alev[3][2] =
        cfac * Om[j][7] * std::exp(-_en[j][7] / temperature) / _sw[j][2];
    alev[3][3] =
        -(_ea[j][2] + _ea[j][5] + _ea[j][7] +
          (cfac / _sw[j][3]) * (Om[j][2] + Om[j][5] + Om[j][7] +
                                Om[j][9] * std::exp(-_en[j][9] / temperature)));
    alev[3][4] = _ea[j][9] + cfac * Om[j][9] / _sw[j][4];

    alev[4][0] =
        cfac * Om[j][3] * std::exp(-_en[j][3] / temperature) / _sw[j][0];
    alev[4][1] =
        cfac * Om[j][6] * std::exp(-_en[j][6] / temperature) / _sw[j][1];
    alev[4][2] =
        cfac * Om[j][8] * std::exp(-_en[j][8] / temperature) / _sw[j][2];
    alev[4][3] =
        cfac * Om[j][9] * std::exp(-_en[j][9] / temperature) / _sw[j][3];
    alev[4][4] =
        -(_ea[j][3] + _ea[j][6] + _ea[j][8] + _ea[j][9] +
          (cfac / _sw[j][4]) * (Om[j][3] + Om[j][6] + Om[j][8] + Om[j][9]));

    // find level populations
    simq(alev, lev);

    double cl2 = abundances[j] * kb * lev[1] * _ea[j][0] * _en[j][0];
    double cl3 = abundances[j] * kb * lev[2] *
                 (_ea[j][1] * _en[j][1] + _ea[j][4] * _en[j][4]);
    double cl4 =
        abundances[j] * kb * lev[3] *
        (_ea[j][2] * _en[j][2] + _ea[j][5] * _en[j][5] + _ea[j][7] * _en[j][7]);
    double cl5 = abundances[j] * kb * lev[4] *
                 (_ea[j][3] * _en[j][3] + _ea[j][6] * _en[j][6] +
                  _ea[j][8] * _en[j][8] + _ea[j][9] * _en[j][9]);

    status("%u: %g %g %g %g %g", j, lev[0], lev[1], lev[2], lev[3], lev[4]);

    if (j == 1) {
      c5755 = abundances[j] * kb * lev[4] * _ea[j][9] * _en[j][9];
      c6584 = abundances[j] * kb * lev[3] * _ea[j][7] * _en[j][7];
      cnii122 = abundances[j] * kb * lev[2] * _ea[j][4] * _en[j][4];
    }
    if (j == 2) {
      c6300 = abundances[j] * kb * lev[3] *
              (_ea[j][2] * _en[j][2] + _ea[j][5] * _en[j][5]);
      status("%g %g %g %g %g %g %g %g", abundances[j], kb, lev[3], _ea[j][2],
             _en[j][2], _ea[j][5], _en[j][5], c6300);
    }
    if (j == 3) {
      c3729 = cl2;
      c3727 = cl2 + cl3;
      coii7325 = abundances[j] * kb *
                 (lev[4] * (_ea[j][6] * _en[j][6] + _ea[j][8] * _en[j][8]) +
                  lev[3] * (_ea[j][5] * _en[j][5] + _ea[j][7] * _en[j][7]));
    }
    if (j == 4) {
      c4363 = abundances[j] * kb * lev[4] * _ea[j][9] * _en[j][9];
      c5007 = abundances[j] * kb * lev[3] * _ea[j][7] * _en[j][7];
      c52mu = abundances[j] * kb * lev[2] * _ea[j][4] * _en[j][4];
      c88mu = abundances[j] * kb * lev[1] * _ea[j][0] * _en[j][0];
    }
    if (j == 5) {
      c3869 = abundances[j] * kb * lev[3] * _ea[j][2] * _en[j][2];
      cneiii15 = cl2;
    }
    if (j == 6) {
      c4072 = abundances[j] * kb *
              (lev[3] * _ea[j][2] * _en[j][2] + lev[4] * _ea[j][3] * _en[j][3]);
      c6717 = cl3;
      c6725 = cl2 + cl3;
    }
    if (j == 7) {
      c9405 = abundances[j] * kb * lev[3] *
              (_ea[j][5] * _en[j][5] + _ea[j][7] * _en[j][7]);
      c6312 = abundances[j] * kb * lev[4] * _ea[j][9] * _en[j][9];
      c33mu = abundances[j] * kb * lev[1] * _ea[j][0] * _en[j][0];
      c19mu = abundances[j] * kb * lev[2] * _ea[j][4] * _en[j][4];
    }
    if (j == 8) {
      cii2325 = cl3 + cl4 + cl5;
    }
    if (j == 9) {
      ciii1908 = cl2 + cl3 + cl4;
    }
  }

  // 2 level atoms
  double sw1 = 2.;
  double sw2 = 4.;
  T1 = std::exp(-EnNIII / temperature);
  cniii57 = abundances[10] * kb * cfac * EnNIII * OmNIII * T1 * EaNIII /
            (sw1 * (EaNIII + cfac * OmNIII * (1. / sw2 + T1 / sw1)));
  sw1 = 4.;
  sw2 = 2.;
  T1 = std::exp(-EnNeII / temperature);
  cneii12 = abundances[11] * kb * cfac * OmNeII * EnNeII * T1 * EaNeII /
            (sw1 * (EaNeII + cfac * OmNeII * (1. / sw2 + T1 / sw1)));
}
