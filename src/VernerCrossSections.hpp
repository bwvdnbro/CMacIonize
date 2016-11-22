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

/**
 * @brief CrossSections implementation for Verner's cross sections.
 */
class VernerCrossSections : public CrossSections {
private:
  /*! @brief L array from Verner's script. */
  unsigned char _L[7];

  /*! @brief NINN array from Verner's script. */
  unsigned char _NINN[30];

  /*! @brief NTOT array from Verner's script. */
  unsigned char _NTOT[30];

  /*! @brief PH1 array from Verner's script. */
  double _PH1[6][30][30][7];

  /*! @brief PH2 array from Verner's script. */
  double _PH2[7][30][30];

public:
  VernerCrossSections();

  double get_cross_section_verner(unsigned char nz, unsigned char ne,
                                  unsigned char is, double e);

  virtual double get_cross_section(IonName ion, double energy);
};

#endif // CROSSSECTIONS_HPP
