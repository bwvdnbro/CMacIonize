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
 * @file EmissivityValues.hpp
 *
 * @brief Emissivity values in a cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EMISSIVITYVALUES_HPP
#define EMISSIVITYVALUES_HPP

#include "Error.hpp"

#include <cinttypes>
#include <string>

/**
 * @brief Names for emission lines.
 */
enum EmissionLine {
  EMISSIONLINE_HAlpha = 0,
  EMISSIONLINE_HBeta,
  EMISSIONLINE_HII,
  EMISSIONLINE_BALMER_JUMP_LOW,
  EMISSIONLINE_BALMER_JUMP_HIGH,
  EMISSIONLINE_OI_6300,
  EMISSIONLINE_OI_6364,
  EMISSIONLINE_OII_3727,
  EMISSIONLINE_OIII_5007,
  EMISSIONLINE_OIII_4959,
  EMISSIONLINE_OIII_4363,
  EMISSIONLINE_OIII_52mu,
  EMISSIONLINE_OIII_88mu,
  EMISSIONLINE_NII_5755,
  EMISSIONLINE_NII_6548,
  EMISSIONLINE_NII_6584,
  EMISSIONLINE_NeIII_3869,
  EMISSIONLINE_NeIII_3968,
  EMISSIONLINE_SII_6725,
  EMISSIONLINE_SII_4072,
  EMISSIONLINE_SIII_9405,
  EMISSIONLINE_SIII_6312,
  EMISSIONLINE_SIII_19mu,
  EMISSIONLINE_SIII_33mu,
  EMISSIONLINE_avg_T,
  EMISSIONLINE_avg_T_count,
  EMISSIONLINE_avg_nH_nHe,
  EMISSIONLINE_avg_nH_nHe_count,
  EMISSIONLINE_NeII_12mu,
  EMISSIONLINE_NIII_57mu,
  EMISSIONLINE_NeIII_15mu,
  EMISSIONLINE_NII_122mu,
  EMISSIONLINE_CII_158mu,
  EMISSIONLINE_CII_2325,
  EMISSIONLINE_CIII_1908,
  EMISSIONLINE_OII_7325,
  EMISSIONLINE_SIV_10mu,
  EMISSIONLINE_HeI_5876,
  EMISSIONLINE_Hrec_s,
  EMISSIONLINE_WFC2_F439W,
  EMISSIONLINE_WFC2_F555W,
  EMISSIONLINE_WFC2_F675W,
  NUMBER_OF_EMISSIONLINES
};

/**
 * @brief Emissivity values in a cell.
 */
class EmissivityValues {
private:
  /*! @brief Emissivity values (in J m^-3s^-1). */
  double _emissivities[NUMBER_OF_EMISSIONLINES];

public:
  /**
   * @brief (Empty) constructor.
   */
  inline EmissivityValues() {
    for (int_fast32_t i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
      _emissivities[i] = 0.;
    }
  }

  /**
   * @brief Get the emissivity of a given line.
   *
   * @param line Valid EmissionLine.
   * @return Emissivity of that line (in J m^-3s^-1).
   */
  inline double get_emissivity(const int_fast32_t line) const {
    return _emissivities[line];
  }

  /**
   * @brief Set the emissivity of the given line.
   *
   * @param line Valid EmissionLine.
   * @param emissivity Emissivity of that line (in J m^-3s^-1).
   */
  inline void set_emissivity(const int_fast32_t line, const double emissivity) {
    _emissivities[line] = emissivity;
  }

  /**
   * @brief Get a std::string representation of the given EmissionLine name.
   *
   * @param line Valid EmissionLine.
   * @return std::string containing the name of the EmissionLine.
   */
  static inline std::string get_name(const int_fast32_t line) {
    switch (line) {
    case EMISSIONLINE_HAlpha:
      return "Halpha";
    case EMISSIONLINE_HBeta:
      return "Hbeta";
    case EMISSIONLINE_HII:
      return "HII";
    case EMISSIONLINE_BALMER_JUMP_LOW:
      return "BaLow";
    case EMISSIONLINE_BALMER_JUMP_HIGH:
      return "BaHigh";
    case EMISSIONLINE_OI_6300:
      return "OI_6300";
    case EMISSIONLINE_OI_6364:
      return "OI_6364";
    case EMISSIONLINE_OII_3727:
      return "OII_3727";
    case EMISSIONLINE_OIII_5007:
      return "OIII_5007";
    case EMISSIONLINE_OIII_4959:
      return "OIII_4959";
    case EMISSIONLINE_OIII_4363:
      return "OIII_4363";
    case EMISSIONLINE_OIII_52mu:
      return "OIII_52mu";
    case EMISSIONLINE_OIII_88mu:
      return "OIII_88mu";
    case EMISSIONLINE_NII_5755:
      return "NII_5755";
    case EMISSIONLINE_NII_6548:
      return "NII_6548";
    case EMISSIONLINE_NII_6584:
      return "NII_6584";
    case EMISSIONLINE_NeIII_3869:
      return "NeIII_3869";
    case EMISSIONLINE_NeIII_3968:
      return "NeIII_3968";
    case EMISSIONLINE_SII_6725:
      return "SII_6725";
    case EMISSIONLINE_SII_4072:
      return "SII_4072";
    case EMISSIONLINE_SIII_9405:
      return "SIII_9405";
    case EMISSIONLINE_SIII_6312:
      return "SIII_6213";
    case EMISSIONLINE_SIII_19mu:
      return "SIII_19mu";
    case EMISSIONLINE_SIII_33mu:
      return "SIII_33mu";
    case EMISSIONLINE_avg_T:
      return "avg_T";
    case EMISSIONLINE_avg_T_count:
      return "avg_T_count";
    case EMISSIONLINE_avg_nH_nHe:
      return "avg_nH_nHe";
    case EMISSIONLINE_avg_nH_nHe_count:
      return "avg_nH_nHe_count";
    case EMISSIONLINE_NeII_12mu:
      return "NeII_12mu";
    case EMISSIONLINE_NIII_57mu:
      return "NIII_57mu";
    case EMISSIONLINE_NeIII_15mu:
      return "NeIII_15mu";
    case EMISSIONLINE_NII_122mu:
      return "NII_122mu";
    case EMISSIONLINE_CII_158mu:
      return "CII_158mu";
    case EMISSIONLINE_CII_2325:
      return "CII_2325";
    case EMISSIONLINE_CIII_1908:
      return "CIII_1908";
    case EMISSIONLINE_OII_7325:
      return "OII_7325";
    case EMISSIONLINE_SIV_10mu:
      return "SIV_10mu";
    case EMISSIONLINE_HeI_5876:
      return "HeI_5876";
    case EMISSIONLINE_Hrec_s:
      return "Hrec_s";
    case EMISSIONLINE_WFC2_F439W:
      return "WFC2_F439W";
    case EMISSIONLINE_WFC2_F555W:
      return "WFC2_F555W";
    case EMISSIONLINE_WFC2_F675W:
      return "WFC2_F675W";
    default:
      cmac_error("Unknown EmissionLine: %" PRIiFAST32 "!", line);
      return "";
    }
  }

  /**
   * @brief Get the central wavelength for the given emission line.
   *
   * @param line EmissionLine.
   * @return Corresponding wavelength, in m (or 0 if the line is a pseudo-line).
   */
  static inline double get_central_wavelength(const int_fast32_t line) {

    switch (line) {
    case EMISSIONLINE_HAlpha:
      // NIST (via wikepedia)
      return 656.279e-9;
    case EMISSIONLINE_HBeta:
      // NIST (via wikipedia)
      return 486.135e-9;
    case EMISSIONLINE_HII:
      // pseudo-line (not sure what it is exactly)
      return 0.;
    case EMISSIONLINE_BALMER_JUMP_LOW:
      // pseudo-line (not sure what it is exactly)
      return 0.;
    case EMISSIONLINE_BALMER_JUMP_HIGH:
      // pseudo-line (not sure what it is exactly)
      return 0.;
    case EMISSIONLINE_OI_6300:
      // Osterbrock & Ferland (2006), Table 3.14
      return 6300.3e-10;
    case EMISSIONLINE_OI_6364:
      // Osterbrock & Ferland (2006), Table 3.14
      return 6363.8e-10;
    case EMISSIONLINE_OII_3727:
      // Osterbrock & Ferland (2006), Table 3.13
      // rate average of the 3726.0 and 3728.8 angstrom transitions
      return 3726.5e-10;
    case EMISSIONLINE_OIII_5007:
      // Osterbrock & Ferland (2006), Table 3.12
      return 5006.9e-10;
    case EMISSIONLINE_OIII_4959:
      // Osterbrock & Ferland (2006), Table 3.12
      return 4958.9e-10;
    case EMISSIONLINE_OIII_4363:
      // Osterbrock & Ferland (2006), Table 3.12
      return 4363.2e-10;
    case EMISSIONLINE_OIII_52mu:
      // Osterbrock & Ferland (2006), Table 3.12
      return 51.814e-6;
    case EMISSIONLINE_OIII_88mu:
      // Osterbrock & Ferland (2006), Table 3.12
      return 88.356e-6;
    case EMISSIONLINE_NII_5755:
      // Osterbrock & Ferland (2006), Table 3.12
      return 5754.6e-10;
    case EMISSIONLINE_NII_6548:
      // Osterbrock & Ferland (2006), Table 3.12
      return 6548.0e-10;
    case EMISSIONLINE_NII_6584:
      // Osterbrock & Ferland (2006), Table 3.12
      return 6583.4e-10;
    case EMISSIONLINE_NeIII_3869:
      // Osterbrock & Ferland (2006), Table 3.14
      return 3868.8e-10;
    case EMISSIONLINE_NeIII_3968:
      // Osterbrock & Ferland (2006), Table 3.14
      return 3967.5e-10;
    case EMISSIONLINE_SII_6725:
      // Osterbrock & Ferland (2006), Table 3.13
      // rate average of the 6716.5 and 6730.8 angstrom transition
      return 6727.5e-10;
    case EMISSIONLINE_SII_4072:
      // Osterbrock & Ferland (2006), Table 3.13
      // rate average of the 4068.6 and 4076.4 angstrom transition
      return 4070.8e-10;
    case EMISSIONLINE_SIII_9405:
      // Osterbrock & Ferland (2006), Table 3.12
      // rate average of the 9068.9 and 9531.0 angstrom transition
      return 9403.3e-10;
    case EMISSIONLINE_SIII_6312:
      // Osterbrock & Ferland (2006), Table 3.12
      return 6312.0e-10;
    case EMISSIONLINE_SIII_19mu:
      // Osterbrock & Ferland (2006), Table 3.12
      return 18.713e-6;
    case EMISSIONLINE_SIII_33mu:
      // Osterbrock & Ferland (2006), Table 3.12
      return 33.47e-6;
    case EMISSIONLINE_avg_T:
      // pseudo-line
      return 0.;
    case EMISSIONLINE_avg_T_count:
      // pseudo-line
      return 0.;
    case EMISSIONLINE_avg_nH_nHe:
      // pseudo-line
      return 0.;
    case EMISSIONLINE_avg_nH_nHe_count:
      // pseudo-line
      return 0.;
    case EMISSIONLINE_NeII_12mu:
      // Osterbrock & Ferland (2006), Table 3.11
      return 12.814e-6;
    case EMISSIONLINE_NIII_57mu:
      // Osterbrock & Ferland (2006), Table 3.9
      return 57.343e-6;
    case EMISSIONLINE_NeIII_15mu:
      // Osterbrock & Ferland (2006), Table 3.14
      return 15.555e-6;
    case EMISSIONLINE_NII_122mu:
      // Osterbrock & Ferland (2006), Table 3.12
      return 121.89e-6;
    case EMISSIONLINE_CII_158mu:
      // Osterbrock & Ferland (2006), Table 3.9
      return 157.74e-6;
    case EMISSIONLINE_CII_2325:
      // Osterbrock & Ferland (2006), Table 3.9
      // rate average of the 2325.4, 2326.9, 2328.1, 2323.6 and 2324.7 angstrom
      // transitions
      return 2326.2e-10;
    case EMISSIONLINE_CIII_1908:
      // Osterbrock & Ferland (2006), Table 3.8
      // rate average of the 1906.7 and 1908.7 angstrom transitions
      return 1908.7e-10;
    case EMISSIONLINE_OII_7325:
      // Osterbrock & Ferland (2006), Table 3.13
      // rate average of the 7318.8, 7319.9, 7329.6 and 7330.7 angstrom
      // transitions
      return 7324.5e-10;
    case EMISSIONLINE_SIV_10mu:
      // Osterbrock & Ferland (2006), Table 3.10
      return 10.514e-6;
    case EMISSIONLINE_HeI_5876:
      // Osterbrock & Ferland (2006), Table 4.6
      // no detailed wavelength info available
      return 5876e-10;
    case EMISSIONLINE_Hrec_s:
      // pseudo-line
      return 0.;
    case EMISSIONLINE_WFC2_F439W:
      // pseudo-line
      return 0.;
    case EMISSIONLINE_WFC2_F555W:
      // pseudo-line
      return 0.;
    case EMISSIONLINE_WFC2_F675W:
      // pseudo-line
      return 0.;
    default:
      cmac_error("Unknown EmissionLine: %" PRIiFAST32 "!", line);
      return 0.;
    }
  }
};

#endif // EMISSIVITYVALUES_HPP
