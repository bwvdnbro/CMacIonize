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

#include "Error.hpp"

#include <string>

/**
 * @brief Names of supported atoms.
 *
 * These names are currently only used by Abundances, that is the reason
 * hydrogen is missing from the list. It should be perfectly fine to add it, but
 * why should we if it is not used?
 */
enum ElementName {
  /*! @brief Helium. */
  ELEMENT_He = 0,
  /*! @brief Carbon. */
  ELEMENT_C,
  /*! @brief Nitrogen. */
  ELEMENT_N,
  /*! @brief Oxygen. */
  ELEMENT_O,
  /*! @brief Neon. */
  ELEMENT_Ne,
  /*! @brief Sulphur. */
  ELEMENT_S,
  /*! @brief Atom number counter. Add new atoms above this element! */
  NUMBER_OF_ELEMENTNAMES
};

/**
 * @brief Names of supported ions.
 *
 * These are the atoms/ions that can be ionized by radiation in our program. A
 * suffix 'p' denotes an ion, and is followed by a number denoting the
 * ionization state. Note that the total number of ionization states of a single
 * ion will be one more than listed below, since the highest ionized state
 * is not mentioned (as it can not be ionized).
 */
enum IonName {
  /*! @brief Neutral hydrogen. */
  ION_H_n = 0,
  /*! @brief Neutral helium. */
  ION_He_n,
  /*! @brief Ionized carbon. */
  ION_C_p1,
  /*! @brief Double ionized carbon. */
  ION_C_p2,
  /*! @brief Neutral nytrogen. */
  ION_N_n,
  /*! @brief Ionized nytrogen. */
  ION_N_p1,
  /*! @brief Double ionized nytrogen. */
  ION_N_p2,
  /*! @brief Neutral oxygen. */
  ION_O_n,
  /*! @brief Ionized oxygen. */
  ION_O_p1,
  /*! @brief Neutral neon. */
  ION_Ne_n,
  /*! @brief Ionized neon. */
  ION_Ne_p1,
  /*! @brief Ionized sulphur. */
  ION_S_p1,
  /*! @brief Double ionized sulphur. */
  ION_S_p2,
  /*! @brief Triple ionized sulphur. */
  ION_S_p3,
  /*! @brief Number of supported elements. */
  NUMBER_OF_IONNAMES
};

/**
 * @brief Get the name of the given element.
 *
 * @param element ElementName for a valid element.
 * @return Textual representation of the element name.
 */
static inline std::string get_element_name(int element) {
  switch (element) {
  case ELEMENT_He:
    return "He";
  case ELEMENT_C:
    return "C";
  case ELEMENT_N:
    return "N";
  case ELEMENT_O:
    return "O";
  case ELEMENT_Ne:
    return "Ne";
  case ELEMENT_S:
    return "S";
  default:
    cmac_error("Unknown element: %i!", element);
    return "";
  }
}

/**
 * @brief Get the name of the given ion.
 *
 * @param ion IonName for a valid ion.
 * @return Textual representation of the ion name.
 */
static inline std::string get_ion_name(int ion) {
  switch (ion) {

  case ION_H_n:
    return "H";

  case ION_He_n:
    return "He";

  case ION_C_p1:
    return "C+";
  case ION_C_p2:
    return "C++";

  case ION_N_n:
    return "N";
  case ION_N_p1:
    return "N+";
  case ION_N_p2:
    return "N++";

  case ION_O_n:
    return "O";
  case ION_O_p1:
    return "O+";

  case ION_Ne_n:
    return "Ne";
  case ION_Ne_p1:
    return "Ne+";

  case ION_S_p1:
    return "S+";
  case ION_S_p2:
    return "S++";
  case ION_S_p3:
    return "S+++";

  default:
    cmac_error("Unknown ion: %i!", ion);
    return "";
  }
}

#endif // ELEMENTNAMES_HPP
