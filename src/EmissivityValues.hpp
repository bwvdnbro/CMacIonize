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
  EMISSIONLINE_OII_3727,
  EMISSIONLINE_OIII_5007,
  EMISSIONLINE_OIII_4363,
  EMISSIONLINE_OIII_88mu,
  EMISSIONLINE_NII_5755,
  EMISSIONLINE_NII_6584,
  EMISSIONLINE_NeIII_3869,
  EMISSIONLINE_SII_6725,
  EMISSIONLINE_SII_4072,
  EMISSIONLINE_SIII_9405,
  EMISSIONLINE_SIII_6312,
  EMISSIONLINE_SIII_19mu,
  EMISSIONLINE_avg_T,
  EMISSIONLINE_avg_T_count,
  EMISSIONLINE_avg_nH_nHe,
  EMISSIONLINE_avg_nH_nHe_count,
  EMISSIONLINE_NeII_12mu,
  EMISSIONLINE_NIII_57mu,
  EMISSIONLINE_NeIII_15mu,
  EMISSIONLINE_NII_122mu,
  EMISSIONLINE_CII_2325,
  EMISSIONLINE_CIII_1908,
  EMISSIONLINE_OII_7325,
  EMISSIONLINE_SIV_10mu,
  EMISSIONLINE_HeI_5876,
  EMISSIONLINE_Hrec_s,
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
    for (int i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
      _emissivities[i] = 0.;
    }
  }

  /**
   * @brief Get the emissivity of a given line.
   *
   * @param line Valid EmissionLine.
   * @return Emissivity of that line (in J m^-3s^-1).
   */
  inline double get_emissivity(EmissionLine line) const {
    return _emissivities[line];
  }

  /**
   * @brief Set the emissivity of the given line.
   *
   * @param line Valid EmissionLine.
   * @param emissivity Emissivity of that line (in J m^-3s^-1).
   */
  inline void set_emissivity(EmissionLine line, double emissivity) {
    _emissivities[line] = emissivity;
  }

  /**
   * @brief Get a std::string representation of the given EmissionLine name.
   *
   * @param line Valid EmissionLine.
   * @return std::string containing the name of the EmissionLine.
   */
  static inline std::string get_name(EmissionLine line) {
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
    case EMISSIONLINE_OII_3727:
      return "OII_3727";
    case EMISSIONLINE_OIII_5007:
      return "OIII_5007";
    case EMISSIONLINE_OIII_4363:
      return "OIII_4363";
    case EMISSIONLINE_OIII_88mu:
      return "OIII_88mu";
    case EMISSIONLINE_NII_5755:
      return "NII_5755";
    case EMISSIONLINE_NII_6584:
      return "NII_6584";
    case EMISSIONLINE_NeIII_3869:
      return "NeIII_3869";
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
    default:
      cmac_error("Unknown EmissionLine: %i!", line);
      return "";
    }
  }
};

#endif // EMISSIVITYVALUES_HPP
