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
#include "VernerCrossSectionsDataLocation.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

/**
 * @brief Constructor.
 *
 * Reads in the data from the data files.
 */
VernerCrossSections::VernerCrossSections() {
  std::ifstream fileA(VERNERCROSSSECTIONSDATALOCATION_A);
  std::string lineA;
  while (getline(fileA, lineA)) {
    if (lineA[0] != '#') {
      std::istringstream lstream(lineA);

      unsigned int Z, N, n, l;
      double E_th, E_0, sigma_0, y_a, P, y_w;

      lstream >> Z >> N >> n >> l >> E_th >> E_0 >> sigma_0 >> y_a >> P >> y_w;

      unsigned int iZ = Z - 1;
      unsigned int iN = N - 1;

      // n and l need to be combined to get the shell number
      if (n < 3) {
        n += l;
      } else {
        ++n;
        if (n < 5) {
          n += l;
        } else {
          n += 2;
        }
      }

      unsigned int in = n - 1;
      if (Z > _data_A.size()) {
        _data_A.resize(Z);
      }
      if (N > _data_A[iZ].size()) {
        _data_A[iZ].resize(N);
      }
      if (n > _data_A[iZ][iN].size()) {
        _data_A[iZ][iN].resize(n);
      }
      _data_A[iZ][iN][in].resize(VERNERDATA_A_NUMELEMENTS);
      _data_A[iZ][iN][in][VERNERDATA_A_l] = l;
      _data_A[iZ][iN][in][VERNERDATA_A_E_th] = E_th;
      _data_A[iZ][iN][in][VERNERDATA_A_E_0] = E_0;
      _data_A[iZ][iN][in][VERNERDATA_A_sigma_0] = sigma_0;
      _data_A[iZ][iN][in][VERNERDATA_A_y_a] = y_a;
      _data_A[iZ][iN][in][VERNERDATA_A_P] = P;
      _data_A[iZ][iN][in][VERNERDATA_A_y_w_squared] = y_w * y_w;
    }
  }

  std::ifstream fileB(VERNERCROSSSECTIONSDATALOCATION_B);
  std::string lineB;
  while (getline(fileB, lineB)) {
    if (lineB[0] != '#') {
      std::istringstream lstream(lineB);

      unsigned int Z, N;
      double E_th, E_max, E_0, sigma_0, y_a, P, y_w, y_0, y_1;

      lstream >> Z >> N >> E_th >> E_max >> E_0 >> sigma_0 >> y_a >> P >> y_w >>
          y_0 >> y_1;

      unsigned int iZ = Z - 1;
      unsigned int iN = N - 1;
      if (Z > _data_B.size()) {
        _data_B.resize(Z);
      }
      if (N > _data_B[iZ].size()) {
        _data_B[iZ].resize(N);
      }
      _data_B[iZ][iN].resize(VERNERDATA_B_NUMELEMENTS);
      _data_B[iZ][iN][VERNERDATA_B_E_0] = E_0;
      _data_B[iZ][iN][VERNERDATA_B_sigma_0] = sigma_0;
      _data_B[iZ][iN][VERNERDATA_B_y_a] = y_a;
      _data_B[iZ][iN][VERNERDATA_B_P] = P;
      _data_B[iZ][iN][VERNERDATA_B_y_w_squared] = y_w * y_w;
      _data_B[iZ][iN][VERNERDATA_B_y_0] = y_0;
      _data_B[iZ][iN][VERNERDATA_B_y_1_squared] = y_1 * y_1;
    }
  }

  std::ifstream fileC(VERNERCROSSSECTIONSDATALOCATION_C);
  std::string lineC;
  while (getline(fileC, lineC)) {
    if (lineC[0] != '#') {
      std::istringstream lstream(lineC);

      unsigned int N, Ninn, Ntot;

      lstream >> N >> Ninn >> Ntot;

      unsigned int iN = N - 1;
      if (N > _data_C.size()) {
        _data_C.resize(N);
      }
      _data_C[iN].resize(VERNERDATA_C_NUMELEMENTS);
      _data_C[iN][VERNERDATA_C_Ninn] = Ninn;
      _data_C[iN][VERNERDATA_C_Ntot] = Ntot;
    }
  }
}

/**
 * @brief C++ version of Verner's phfit2 routine, completely rewritten based on
 * Verner's papers and his original code.
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
                                                     double e) const {
  // convert Hz to eV
  e /= (1.5091902e33 * 1.60217662e-19);

  // assertions. Most compilers optimize these out if optimization flags are
  // given.
  assert(nz > 0 && nz <= 30);
  assert(ne > 0 && ne <= nz);

  unsigned int iZ = nz - 1;
  unsigned int iN = ne - 1;
  unsigned int in = is - 1;

  // if the energy is lower than the ionization threshold energy, the cross
  // section is trivially 0
  if (e < _data_A[iZ][iN][in][VERNERDATA_A_E_th]) {
    return 0.;
  }

  // now figure out which fitting formula to use:
  //  - the fit in the range from E_th to E_inn, the inner shell photoionization
  //    cross section jump (_data_B)
  //  - the fit to the inner shell cross sections (_data_A)
  unsigned int nout = _data_C[iN][VERNERDATA_C_Ntot];
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

  unsigned int nint = _data_C[iN][VERNERDATA_C_Ninn];
  double einn;
  if (nz == 15 || nz == 17 || nz == 19 || (nz > 20 && nz != 26)) {
    einn = 0.;
  } else {
    if (ne < 3) {
      einn = 1.e30;
    } else {
      einn = _data_A[iZ][iN][nint - 1][VERNERDATA_A_E_th];
    }
  }

  if (is < nout && is > nint && e < einn) {
    return 0.;
  }

  if (is <= nint || e >= einn) {
    double E_0 = _data_A[iZ][iN][in][VERNERDATA_A_E_0];
    double sigma_0 = _data_A[iZ][iN][in][VERNERDATA_A_sigma_0];
    double y_a = _data_A[iZ][iN][in][VERNERDATA_A_y_a];
    double P = _data_A[iZ][iN][in][VERNERDATA_A_P];
    double y_w_squared = _data_A[iZ][iN][in][VERNERDATA_A_y_w_squared];
    unsigned int l = _data_A[iZ][iN][in][VERNERDATA_A_l];

    double y = e / E_0;
    double ym1 = y - 1.;
    double Fy = (ym1 * ym1 + y_w_squared) * std::pow(y, 0.5 * P - 5.5 - l) *
                std::pow(1. + std::sqrt(y / y_a), -P);
    return 1.e-22 * sigma_0 * Fy;
  } else {
    double E_0 = _data_B[iZ][iN][VERNERDATA_B_E_0];
    double y_0 = _data_B[iZ][iN][VERNERDATA_B_y_0];
    double y_1_squared = _data_B[iZ][iN][VERNERDATA_B_y_1_squared];
    double y_w_squared = _data_B[iZ][iN][VERNERDATA_B_y_w_squared];
    double P = _data_B[iZ][iN][VERNERDATA_B_P];
    double y_a = _data_B[iZ][iN][VERNERDATA_B_y_a];
    double sigma_0 = _data_B[iZ][iN][VERNERDATA_B_sigma_0];

    double x = e / E_0 - y_0;
    double y = std::sqrt(x * x + y_1_squared);
    double xm1 = x - 1.;
    double Fy = (xm1 * xm1 + y_w_squared) * std::pow(y, 0.5 * P - 5.5) *
                std::pow(1. + std::sqrt(y / y_a), -P);
    return 1.e-22 * sigma_0 * Fy;
  }
}

/**
 * @brief Get the photoionization cross section of the given ion.
 *
 * @param ion IonName of a valid ion.
 * @param energy Photon energy (in Hz).
 * @return Photoionization cross section for the given ion and for the given
 * photon energy (in m^2).
 */
double VernerCrossSections::get_cross_section(IonName ion,
                                              double energy) const {
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
    cmac_error("Unknown ion: %i", ion);
  }
  return 0.;
}
