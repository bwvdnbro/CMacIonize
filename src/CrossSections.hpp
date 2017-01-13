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
 * @file CrossSections.hpp
 *
 * @brief General interface for photoionization cross sections.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CROSSSECTIONS_HPP
#define CROSSSECTIONS_HPP

#include "ElementNames.hpp"

/**
 * @brief General interface for photoionization cross sections.
 */
class CrossSections {
public:
  /**
   * @brief Get the photoionization cross section for the given ion at the
   * given photon energy.
   *
   * @param ion IonName for a valid ion.
   * @param energy Photon frequency (in Hz).
   * @return Photoionization cross section (in m^-2).
   */
  virtual double get_cross_section(IonName ion, double energy) const = 0;
};

#endif // CROSSSECTIONS_HPP
