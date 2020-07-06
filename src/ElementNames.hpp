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

#include <cinttypes>
#include <string>

// activate additional metals
#ifndef HAVE_HYDROGEN_ONLY
#define HAS_HELIUM
#define HAS_CARBON
#define HAS_NITROGEN
#define HAS_OXYGEN
#define HAS_NEON
#define HAS_SULPHUR
#endif

// activate new additional metals
#ifdef ADDITIONAL_COOLANTS
#define HAS_ARGON
#endif

/**
 * @brief Names of supported atoms.
 *
 * These names are currently only used by Abundances, that is the reason
 * hydrogen is missing from the list. It should be perfectly fine to add it, but
 * why should we if it is not used?
 */
enum ElementName {
#ifdef HAS_HELIUM
  /*! @brief Helium. */
  ELEMENT_He = 0,
#else
  /*! @brief Hydrogen. */
  ELEMENT_H = 0,
#endif
#ifdef HAS_CARBON
  /*! @brief Carbon. */
  ELEMENT_C,
#endif
#ifdef HAS_NITROGEN
  /*! @brief Nitrogen. */
  ELEMENT_N,
#endif
#ifdef HAS_OXYGEN
  /*! @brief Oxygen. */
  ELEMENT_O,
#endif
#ifdef HAS_NEON
  /*! @brief Neon. */
  ELEMENT_Ne,
#endif
#ifdef HAS_SULPHUR
  /*! @brief Sulphur. */
  ELEMENT_S,
#endif
#ifdef HAS_ARGON
  ELEMENT_Ar,
#endif
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
  ION_H_n,
#ifdef HAS_HELIUM
  /*! @brief Neutral helium. */
  ION_He_n,
#endif
#ifdef HAS_CARBON
  /*! @brief Ionized carbon. */
  ION_C_p1,
  /*! @brief Double ionized carbon. */
  ION_C_p2,
#endif
#ifdef HAS_NITROGEN
  /*! @brief Neutral nytrogen. */
  ION_N_n,
  /*! @brief Ionized nytrogen. */
  ION_N_p1,
  /*! @brief Double ionized nytrogen. */
  ION_N_p2,
#endif
#ifdef HAS_OXYGEN
  /*! @brief Neutral oxygen. */
  ION_O_n,
  /*! @brief Ionized oxygen. */
  ION_O_p1,
#endif
#ifdef HAS_NEON
  /*! @brief Neutral neon. */
  ION_Ne_n,
  /*! @brief Ionized neon. */
  ION_Ne_p1,
#endif
#ifdef HAS_SULPHUR
  /*! @brief Ionized sulphur. */
  ION_S_p1,
  /*! @brief Double ionized sulphur. */
  ION_S_p2,
  /*! @brief Triple ionized sulphur. */
  ION_S_p3,
#endif
#ifdef HAS_ARGON
  /*! @brief Neutral argon. */
  ION_Ar_n,
  /*! @brief Ionized argon. */
  ION_Ar_p1,
  /*! @brief Double ionized argon. */
  ION_Ar_p2,
  /*! @brief Triple ionized argon. */
  ION_Ar_p3,
#endif
  /*! @brief Number of supported elements. */
  NUMBER_OF_IONNAMES
};

/**
 * @brief Get the name of the given element.
 *
 * @param element ElementName for a valid element.
 * @return Textual representation of the element name.
 */
static inline std::string get_element_name(const int_fast32_t element) {

  switch (element) {
#ifdef HAS_HELIUM
  case ELEMENT_He:
    return "He";
#else
  case ELEMENT_H:
    return "H";
#endif
#ifdef HAS_CARBON
  case ELEMENT_C:
    return "C";
#endif
#ifdef HAS_NITROGEN
  case ELEMENT_N:
    return "N";
#endif
#ifdef HAS_OXYGEN
  case ELEMENT_O:
    return "O";
#endif
#ifdef HAS_NEON
  case ELEMENT_Ne:
    return "Ne";
#endif
#ifdef HAS_SULPHUR
  case ELEMENT_S:
    return "S";
#endif
#ifdef HAS_ARGON
  case ELEMENT_Ar:
    return "Ar";
#endif
  default:
    cmac_error("Unknown element: %" PRIiFAST32 "!", element);
    return "";
  }
}

/**
 * @brief Get the name of the given ion.
 *
 * @param ion IonName for a valid ion.
 * @return Textual representation of the ion name.
 */
static inline std::string get_ion_name(const int_fast32_t ion) {

  switch (ion) {

  case ION_H_n:
    return "H";

#ifdef HAS_HELIUM
  case ION_He_n:
    return "He";
#endif

#ifdef HAS_CARBON
  case ION_C_p1:
    return "C+";
  case ION_C_p2:
    return "C++";
#endif

#ifdef HAS_NITROGEN
  case ION_N_n:
    return "N";
  case ION_N_p1:
    return "N+";
  case ION_N_p2:
    return "N++";
#endif

#ifdef HAS_OXYGEN
  case ION_O_n:
    return "O";
  case ION_O_p1:
    return "O+";
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    return "Ne";
  case ION_Ne_p1:
    return "Ne+";
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1:
    return "S+";
  case ION_S_p2:
    return "S++";
  case ION_S_p3:
    return "S+++";
#endif

#ifdef HAS_ARGON
  case ION_Ar_n:
    return "Ar";
  case ION_Ar_p1:
    return "Ar+";
  case ION_Ar_p2:
    return "Ar++";
  case ION_Ar_p3:
    return "Ar+++";
#endif

  default:
    cmac_error("Unknown ion: %" PRIiFAST32 "!", ion);
    return "";
  }
}

/**
 * @brief Get the element corresponding to the given ion.
 *
 * @param ion IonName.
 * @return Corresponding ElementName.
 */
static inline int_fast32_t get_element(const int_fast32_t ion) {
  switch (ion) {

  case ION_H_n:
    cmac_error("Hydrogen ions don't have a corresponding element!");
    return -1;

#ifdef HAS_HELIUM
  case ION_He_n:
    return ELEMENT_He;
#endif

#ifdef HAS_CARBON
  case ION_C_p1:
    return ELEMENT_C;
  case ION_C_p2:
    return ELEMENT_C;
#endif

#ifdef HAS_NITROGEN
  case ION_N_n:
    return ELEMENT_N;
  case ION_N_p1:
    return ELEMENT_N;
  case ION_N_p2:
    return ELEMENT_N;
#endif

#ifdef HAS_OXYGEN
  case ION_O_n:
    return ELEMENT_O;
  case ION_O_p1:
    return ELEMENT_O;
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    return ELEMENT_Ne;
  case ION_Ne_p1:
    return ELEMENT_Ne;
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1:
    return ELEMENT_S;
  case ION_S_p2:
    return ELEMENT_S;
  case ION_S_p3:
    return ELEMENT_S;
#endif

#ifdef HAS_ARGON
  case ION_Ar_n:
    return ELEMENT_Ar;
  case ION_Ar_p1:
    return ELEMENT_Ar;
  case ION_Ar_p2:
    return ELEMENT_Ar;
  case ION_Ar_p3:
    return ELEMENT_Ar;
#endif

  default:
    cmac_error("Unknown ion: %" PRIiFAST32 "!", ion);
    return -1;
  }
}

#endif // ELEMENTNAMES_HPP
