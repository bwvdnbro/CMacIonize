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
 * @param alev Coefficients of the equations.
 * @param lev Resulting solutions.
 */
void LineCoolingData::simq(double alev[5][5], double lev[5]) {
  // the Fortran code uses strange logic to access the 5x5 array as a 1D array
  // we cannot apply the same logic in C++, as the ordering of arrays is
  // different...
  // for simplicity, we therefore map the 5x5 array onto the equivalent 1D
  // Fortran array
  double A[25];
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < 5; ++j) {
      A[5 * j + i] = alev[i][j];
    }
  }

  /// attempt at translating the fortran code into readable C++
  //  for(unsigned int j = 0; j < 5; ++j){
  //    // find the column with the maximum coefficient
  //    unsigned int imax = 0;
  //    double Amax = 0.;
  //    for(unsigned int i = 0; i < 5; ++i){
  //      if(std::abs(A[j][i]) > Amax){
  //        Amax = A[j][i];
  //        imax = i;
  //      }
  //    }
  //    // check that the matrix is non-singular
  //    if(Amax == 0.){
  //      error("Singular matrix given to simq!");
  //    }
  //    // imax now contains the index of the column with the largest
  //    coefficient
  //    // interchange rows if necessary
  //    for(unsigned int k = 0; k < 5; ++k){
  //      if(imax != i){
  //        double save = A[i][k];
  //        A[i][k] = A[imax][k];
  //        A[imax][k] = save;
  //      }
  //      A[i][k] /= Amax;
  //    }
  //    if(imax != i){
  //      double save = lev[i];
  //      lev[i] = lev[imax];
  //      lev[imax] = save;
  //    }
  //    lev[i] /= Amax;
  //    if(j < 5){

  //      int iqs = 5*(j-1);
  //      for(unsigned int ix = jy; ix < 6; ++ix){
  //        int ixj = iqs + ix;
  //        it = j-ix;
  //        for(unsigned int jx = jy; jx < 6; ++jx){
  //          int ixjx = 5*(jx-1)+ix;
  //          int jjx = ixjx + it;
  //          A[ixjx-1] = A[ixjx-1] - (A[ixj-1]*A[jjx-1]);
  //        }
  //        lev[ix-1] = lev[ix-1] - (lev[j-1]*A[ixj-1]);
  //      }
  //    }
  //  }

  double tol = 0.;
  int jj = -5;
  for (unsigned int j = 1; j < 6; ++j) {
    unsigned int jy = j + 1;
    jj += 6;
    double bigA = 0.;
    unsigned int imax = 0;
    int it = jj - j;
    for (unsigned int i = j; i < 6; ++i) {
      int ij = it + i;
      if (std::abs(bigA) < std::abs(A[ij - 1])) {
        bigA = A[ij - 1];
        imax = i;
      }
    }
    if (std::abs(bigA) <= tol) {
      error("Just an error, no idea what it means...");
    }
    unsigned int i1 = j + 5 * (j - 2);
    it = imax - j;
    for (unsigned int k = j; k < 6; ++k) {
      i1 += 5;
      int i2 = i1 + it;
      double save = A[i1 - 1];
      A[i1 - 1] = A[i2 - 1];
      A[i2 - 1] = save;
      A[i1 - 1] /= bigA;
    }
    double save = lev[imax - 1];
    lev[imax - 1] = lev[j - 1];
    lev[j - 1] = save / bigA;
    if (j == 5) {
      it = 25;
      for (unsigned int l = 1; l < 5; ++l) {
        int ia = it - l;
        int ib = 5 - l;
        int ic = 5;
        for (unsigned int k = 1; k < l + 1; ++k) {
          lev[ib - 1] = lev[ib - 1] - A[ia - 1] * lev[ic - 1];
          ia -= 5;
          --ic;
        }
      }
      return;
    }
    int iqs = 5 * (j - 1);
    for (unsigned int ix = jy; ix < 6; ++ix) {
      int ixj = iqs + ix;
      it = j - ix;
      for (unsigned int jx = jy; jx < 6; ++jx) {
        int ixjx = 5 * (jx - 1) + ix;
        int jjx = ixjx + it;
        A[ixjx - 1] = A[ixjx - 1] - (A[ixj - 1] * A[jjx - 1]);
      }
      lev[ix - 1] = lev[ix - 1] - (lev[j - 1] * A[ixj - 1]);
    }
  }
}

/**
 * @brief Get the radiative energy losses due to line cooling at the given
 * temperature, electron density and coolant abundances.
 *
 * @param temperature Temperature (in K).
 * @param electron_density Electron density (in m^-3).
 * @param abundances Abdunances of coolants.
 * @return Radiative cooling in (kg m^-1s^-3).
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

  double cfac = 8.63e-6 * electron_density / std::sqrt(temperature);
  double T4 = temperature * 1.e-4;
  double A1 = std::pow(T4, 0.91);
  double A2 = std::pow(T4, 1.11);
  double A3 = std::pow(T4, 0.8);

  double Om[10][10];
  Om[0][0] = 0.29 * A1;
  Om[0][1] = 0.194 * A1;
  Om[0][2] = 0.0567 * A1;
  Om[0][3] = 0.113 * A1;
  Om[0][4] = 0.269 * A2;
  Om[0][5] = 0.109 * A3;
  Om[0][6] = 0.266 * std::pow(T4, 0.715);
  Om[0][7] = 0.097 * std::pow(T4, 0.69);
  Om[0][8] = 0.147 * A3;
  Om[0][9] = 0.071 * A2;

  double T1 = 0.266 * A2;
  double T2 = 0.0324 * A2;
  Om[2][0] = 0.0987 * A2;
  Om[2][1] = 0.0292 * T4;
  Om[2][2] = 0.55556 * T1;
  Om[2][3] = 0.55556 * T2;
  Om[2][4] = 0.0264 * std::pow(T4, 1.24);
  Om[2][5] = 0.333 * T1;
  Om[2][6] = 0.333 * T2;
  Om[2][7] = T1 / 9.;
  Om[2][8] = T2 / 9.;
  Om[2][9] = 0.105 * std::pow(T4, 0.52);

  double cooling = 0.;
  double alev[5][5], lev[5];
  for (unsigned int j = 0; j < 10; ++j) {
    for (unsigned int mm = 0; mm < 5; ++mm) {
      alev[0][mm] = 1.;
      lev[mm] = 0.;
    }
    lev[0] = 1.;
    // for the attentive reader: this effectively overwrites the values assigned
    // above. No idea what we're doing...
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
        cfac * Om[j][3] * std::exp(-_en[j][2] / temperature) / _sw[j][0];
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
  T1 = std::exp(-EnNIII / temperature);
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

  return cooling;
}
