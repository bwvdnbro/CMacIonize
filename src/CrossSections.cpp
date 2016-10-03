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
 * @file CrossSections.cpp
 *
 * @brief Scattering cross sections: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CrossSections.hpp"
#include "CrossSectionsDataLocation.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

/**
 * @brief Constructor.
 *
 * Reads in the data from a file.
 */
CrossSections::CrossSections() {
  ifstream file(CROSSSECTIONSDATALOCATION);
  string line;
  // skip initial comment line
  getline(file, line);
  // read L values
  {
    // skip comment line
    getline(file, line);
    getline(file, line);
    stringstream linestream(line);
    for (unsigned int i = 0; i < 7; ++i) {
      // we cannot simply use linestream >> _L[i], since then the integer value
      // is converted to a char (so we add the numeric value of '0' to all
      // L values...)
      unsigned int integer;
      linestream >> integer;
      _L[i] = integer;
    }
  }

  // read NINN values
  {
    // skip comment line
    getline(file, line);
    getline(file, line);
    stringstream linestream(line);
    for (unsigned int i = 0; i < 30; ++i) {
      unsigned int integer;
      linestream >> integer;
      _NINN[i] = integer;
    }
  }

  // read NTOT values
  {
    // skip comment line
    getline(file, line);
    getline(file, line);
    stringstream linestream(line);
    for (unsigned int i = 0; i < 30; ++i) {
      unsigned int integer;
      linestream >> integer;
      _NTOT[i] = integer;
    }
  }

  // read PH1 values
  {
    // skip comment line
    getline(file, line);
    // we have a triple loop of lines
    for (unsigned int i1 = 0; i1 < 6; i1++) {
      for (unsigned int i2 = 0; i2 < 30; i2++) {
        for (unsigned int i3 = 0; i3 < 30; i3++) {
          getline(file, line);
          stringstream linestream(line);
          for (unsigned int i4 = 0; i4 < 7; ++i4) {
            linestream >> _PH1[i1][i2][i3][i4];
          }
        }
      }
    }
  }

  // read PH2 values
  {
    // skip comment line
    getline(file, line);
    // we have a double loop
    for (unsigned int i1 = 0; i1 < 7; ++i1) {
      for (unsigned int i2 = 0; i2 < 30; ++i2) {
        getline(file, line);
        stringstream linestream(line);
        for (unsigned int i3 = 0; i3 < 30; ++i3) {
          linestream >> _PH2[i1][i2][i3];
        }
      }
    }
  }
}

/**
 * @brief C++ version of Verner's phfit2 routine.
 *
 * @param nz Atomic number.
 * @param ne Number of electrons.
 * @param is Shell number.
 * @param e Photon energy (in eV).
 * @return Photoionization cross section.
 */
double CrossSections::get_cross_section_verner(unsigned char nz,
                                               unsigned char ne,
                                               unsigned char is, double e) {
  double s = 0.;
  if (nz < 1 || nz > 30) {
    return 0.;
  }
  if (ne < 1 || ne > nz) {
    return 0.;
  }
  // Fortran counts from 1, we count from 0
  unsigned char nout = _NTOT[ne - 1];
  if (nz == ne && nz > 18) {
    nout = 7;
  }
  if (nz == (ne + 1) &&
      (nz == 20 || nz == 21 || nz == 22 || nz == 25 || nz == 26)) {
    nout = 7;
  }
  if (is > nout) {
    return 0.;
  }
  if (e < _PH1[0][nz - 1][ne - 1][is - 1]) {
    return 0.;
  }
  unsigned char nint = _NINN[ne - 1];
  double einn;
  if (nz == 15 || nz == 17 || nz == 19 || (nz > 20 && nz != 26)) {
    einn = 0.;
  } else {
    if (ne < 3) {
      einn = 1.e30;
    } else {
      einn = _PH1[0][nz - 1][ne - 1][nint - 1];
    }
  }
  if (is < nout && is > nint && e < einn) {
    return 0.;
  }
  if (is <= nint || e >= einn) {
    double p1 = -_PH1[4][nz - 1][ne - 1][is - 1];
    double y = e / _PH1[1][nz - 1][ne - 1][is - 1];
    double q = -0.5 * p1 - _L[is - 1] - 5.5;
    double a =
        _PH1[2][nz - 1][ne - 1][is - 1] *
        ((y - 1.) * (y - 1.) +
         _PH1[5][nz - 1][ne - 1][is - 1] * _PH1[5][nz - 1][ne - 1][is - 1]);
    double b = sqrt(y / _PH1[3][nz - 1][ne - 1][is - 1]) + 1.;
    s = a * pow(y, q) * pow(b, p1);
  } else {
    double p1 = -_PH2[3][nz - 1][ne - 1];
    double q = -0.5 * p1 - 5.5;
    double x = e / _PH2[0][nz - 1][ne - 1] - _PH2[5][nz - 1][ne - 1];
    double z = sqrt(x * x + _PH2[6][nz - 1][ne - 1] * _PH2[6][nz - 1][ne - 1]);
    double a = _PH2[1][nz - 1][ne - 1] *
               ((x - 1.) * (x - 1.) +
                _PH2[4][nz - 1][ne - 1] * _PH2[4][nz - 1][ne - 1]);
    double b = 1. + sqrt(z / _PH2[2][nz - 1][ne - 1]);
    s = a * pow(z, q) * pow(b, p1);
  }
  return s;
}

/**
 * @brief Get the photoionization cross section of the given element.
 *
 * @param element CrossSectionElements index of an element.
 * @param energy Photon energy.
 * @return Photoionization cross section for the given element and for the given
 * photon energy.
 */
double CrossSections::get_cross_section(int element, double energy) {
  switch (element) {
  case ELEMENT_H:
    return get_cross_section_verner(1, 1, 1, energy);
  case ELEMENT_He:
    return get_cross_section_verner(2, 2, 1, energy);
  }
  return 0.;
}
