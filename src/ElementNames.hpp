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
 * @file ElementNames.hpp
 *
 * @brief Names of elements used in cross section and recombination rate
 * routines.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ELEMENTNAMES_HPP
#define ELEMENTNAMES_HPP

/**
 * @brief Names of supported elements.
 *
 * These are the atoms/ions that can be ionized by radiation in our program. A
 * suffix 'p' denotes an ion, and is followed by a number denoting the
 * ionization state. Note that the total number of ionization states of a single
 * element will be one more than listed below, since the highest ionized state
 * is not mentioned (as it can not be ionized).
 */
enum ElementName {
  /*! @brief Neutral hydrogen. */
  ELEMENT_H = 0,
  /*! @brief Neutral helium. */
  ELEMENT_He,
  /*! @brief Ionized carbon. */
  ELEMENT_Cp1,
  /*! @brief Double ionized carbon. */
  ELEMENT_Cp2,
  /*! @brief Neutral nytrogen. */
  ELEMENT_N,
  /*! @brief Ionized nytrogen. */
  ELEMENT_Np1,
  /*! @brief Double ionized nytrogen. */
  ELEMENT_Np2,
  /*! @brief Neutral oxygen. */
  ELEMENT_O,
  /*! @brief Ionized oxygen. */
  ELEMENT_Op1,
  /*! @brief Neutral neon. */
  ELEMENT_Ne,
  /*! @brief Ionized neon. */
  ELEMENT_Nep1,
  /*! @brief Ionized sulfur. */
  ELEMENT_Sp1,
  /*! @brief Double ionized sulfur. */
  ELEMENT_Sp2,
  /*! @brief Triple ionized sulfur. */
  ELEMENT_Sp3,
  /*! @brief Number of supported elements. */
  NUMBER_OF_ELEMENTS
};

#endif // ELEMENTNAMES_HPP
