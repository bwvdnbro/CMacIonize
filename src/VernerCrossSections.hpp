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
 * @file VernerCrossSections.hpp
 *
 * @brief Verner photoionization cross sections: header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VERNERCROSSSECTIONS_HPP
#define VERNERCROSSSECTIONS_HPP

#include "CrossSections.hpp"

#include <vector>

/**
 * @brief Enum used to identify the columns in verner_A.txt.
 */
enum VernerDataA {
  VERNERDATA_A_l = 0,
  VERNERDATA_A_E_th,
  VERNERDATA_A_E_0,
  VERNERDATA_A_sigma_0,
  VERNERDATA_A_y_a,
  VERNERDATA_A_P,
  VERNERDATA_A_y_w_squared,
  /*! @brief Counter. Should always be the last element. */
  VERNERDATA_A_NUMELEMENTS
};

/**
 * @brief Enum used to identify the columns in verner_B.txt.
 */
enum VernerDataB {
  VERNERDATA_B_E_0 = 0,
  VERNERDATA_B_sigma_0,
  VERNERDATA_B_y_a,
  VERNERDATA_B_P,
  VERNERDATA_B_y_w_squared,
  VERNERDATA_B_y_0,
  VERNERDATA_B_y_1_squared,
  /*! @brief Counter. Should always be the last element. */
  VERNERDATA_B_NUMELEMENTS
};

/**
 * @brief Enum used to identify the columns in verner_C.txt.
 */
enum VernerDataC {
  VERNERDATA_C_Ninn = 0,
  VERNERDATA_C_Ntot,
  /*! @brief Counter. Should always be the last element. */
  VERNERDATA_C_NUMELEMENTS
};

/**
 * @brief CrossSections implementation for Verner's cross sections.
 */
class VernerCrossSections : public CrossSections {
private:
  /*! @brief Verner data A: inner shell cross sections. */
  std::vector< std::vector< std::vector< std::vector< double > > > > _data_A;

  /*! @brief Verner data B: better fits for low energies. */
  std::vector< std::vector< std::vector< double > > > _data_B;

  /*! @brief Verner data C: inner and outer shell numbers used to decide which
   *  fit to use. */
  std::vector< std::vector< double > > _data_C;

public:
  VernerCrossSections();

  double get_cross_section_verner(unsigned char nz, unsigned char ne,
                                  unsigned char is, double e);

  virtual double get_cross_section(IonName ion, double energy);
};

#endif // CROSSSECTIONS_HPP
