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
 * @file VernerCrossSections.cpp
 *
 * @brief Verner photoionization cross sections: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "VernerCrossSections.hpp"
#include "Error.hpp"
#include "UnitConverter.hpp"
#include "VernerCrossSectionsDataLocation.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

/**
 * @brief Constructor.
 *
 * Reads in the data from the data file.
 */
VernerCrossSections::VernerCrossSections() {
  ifstream file(VERNERCROSSSECTIONSDATALOCATION);
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
 * @param e Photon energy (in Hz).
 * @return Photoionization cross section (in m^2).
 */
double VernerCrossSections::get_cross_section_verner(unsigned char nz,
                                                     unsigned char ne,
                                                     unsigned char is,
                                                     double e) {
  e = UnitConverter< QUANTITY_FREQUENCY >::to_unit(e, "eV");

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
  return UnitConverter< QUANTITY_SURFACE_AREA >::to_SI(1.e-18 * s, "cm^2");
}

/**
 * @brief Get the photoionization cross section of the given ion.
 *
 * @param ion IonName of a valid ion.
 * @param energy Photon energy (in Hz).
 * @return Photoionization cross section for the given ion and for the given
 * photon energy (in m^2).
 */
double VernerCrossSections::get_cross_section(IonName ion, double energy) {
  switch (ion) {

  case ION_H_n:
    return get_cross_section_verner(1, 1, 1, energy);

  case ION_He_n:
    return get_cross_section_verner(2, 2, 1, energy);

  case ION_C_p1:
    return get_cross_section_verner(6, 5, 3, energy) +
           get_cross_section_verner(6, 5, 2, energy);
  case ION_C_p2:
    return get_cross_section_verner(6, 4, 2, energy);

  case ION_N_n:
    return get_cross_section_verner(7, 7, 3, energy) +
           get_cross_section_verner(7, 7, 2, energy);
  case ION_N_p1:
    return get_cross_section_verner(7, 6, 3, energy) +
           get_cross_section_verner(7, 6, 2, energy);
  case ION_N_p2:
    return get_cross_section_verner(7, 5, 3, energy);

  case ION_O_n:
    return get_cross_section_verner(8, 8, 3, energy) +
           get_cross_section_verner(8, 8, 2, energy);
  case ION_O_p1:
    return get_cross_section_verner(8, 7, 3, energy) +
           get_cross_section_verner(8, 7, 2, energy);

  case ION_Ne_n:
    return get_cross_section_verner(10, 10, 3, energy) +
           get_cross_section_verner(10, 10, 2, energy);
  case ION_Ne_p1:
    return get_cross_section_verner(10, 9, 3, energy);

  case ION_S_p1:
    return get_cross_section_verner(16, 15, 5, energy) +
           get_cross_section_verner(16, 15, 4, energy);
  case ION_S_p2:
    return get_cross_section_verner(16, 14, 5, energy) +
           get_cross_section_verner(16, 14, 4, energy);
  case ION_S_p3:
    return get_cross_section_verner(16, 13, 5, energy);

  default:
    error("Unknown ion: %i", ion);
  }
  return 0.;
}
