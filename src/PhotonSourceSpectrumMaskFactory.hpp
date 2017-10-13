/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file PhotonSourceSpectrumMaskFactory.hpp
 *
 * @brief Factory for PhotonSourceSpectrumMask implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCESPECTRUMMASKFACTORY_HPP
#define PHOTONSOURCESPECTRUMMASKFACTORY_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceSpectrumMask.hpp"

/**
 * @brief Factory for PhotonSourceSpectrumMask implementations.
 */
class PhotonSourceSpectrumMaskFactory {
public:
  /**
   * @brief Generate a PhotonSourceSpectrumMask of the type selected in the
   * given ParameterFile.
   *
   * @param role Role the PhotonSourceSpectrum fulfils in the simulation.
   * Parameters are read from the mask block in the corresponding block of the
   * parameter file.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created PhotonSourceSpectrumMask instance.
   * Memory management for the pointer should be handled by the calling routine.
   */
  inline static PhotonSourceSpectrumMask *
  generate(std::string role, ParameterFile &params, Log *log = nullptr) {

    std::string type = params.get_value< std::string >(
        role + ":PhotonSourceSpectrumMask:type", "None");
    if (log) {
      log->write_status("Requested PhotonSourceSpectrumMask for ", role, ": ",
                        type);
    }

    if (type == "None") {
      cmac_error("No types implemented yet!");
      return nullptr;
    } else {
      cmac_error("Unknown PhotonSourceSpectrumMask type: \"%s\"!",
                 type.c_str());
      return nullptr;
    }
  }
};

#endif // PHOTONSOURCESPECTRUMMASKFACTORY_HPP
