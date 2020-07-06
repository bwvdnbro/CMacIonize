/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ElementData.hpp
 *
 * @brief Atomic data for the elements tracked in the code.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ELEMENTDATA_HPP
#define ELEMENTDATA_HPP

#include "ElementNames.hpp"

/**
 * @brief Get the ionization energy of the given ion.
 *
 * The returned value is the energy required to ionize the corresponding ion.
 *
 * @param ion IonName for a valid ion.
 * @return Ionization energy for that ion (in Hz).
 */
static inline double get_ionization_energy(const int_fast32_t ion) {

  switch (ion) {

  case ION_H_n:
    return 3.28810279e+15;

#ifdef HAS_HELIUM
  case ION_He_n:
    return 5.94523574e+15;
#endif

#ifdef HAS_CARBON
  case ION_C_p1:
    return 5.89588678e+15;
  case ION_C_p2:
    return 1.15792700e+16;
#endif

#ifdef HAS_NITROGEN
  case ION_N_n:
    return 3.51435505e+15;
  case ION_N_p1:
    return 7.15759434e+15;
  case ION_N_p2:
    return 1.14732262e+16;
#endif

#ifdef HAS_OXYGEN
  case ION_O_n:
    return 3.29284691e+15;
  case ION_O_p1:
    return 8.49136314e+15;
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    return 5.21432028e+15;
  case ION_Ne_p1:
    return 9.90492110e+15;
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1:
    return 5.64310422e+15;
  case ION_S_p2:
    return 8.41222200e+15;
  case ION_S_p3:
    return 1.14182796e+16;
#endif

#ifdef HAS_ARGON
  case ION_Ar_n:
    return Ar 3.81067612e+15;
  case ION_Ar_p1:
    return 6.68085421e+15;
  case ION_Ar_p2:
    return 9.85093200e+15;
  case ION_Ar_p3:
    return 1.44620580e+16;
#endif

  default:
    cmac_error("Unknown ion: %" PRIiFAST32 "!", ion);
    return 0.;
  }
}

#endif // ELEMENTDATA_HPP
