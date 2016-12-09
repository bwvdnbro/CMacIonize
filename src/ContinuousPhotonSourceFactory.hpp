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
 * @file ContinuousPhotonSourceFactory.hpp
 *
 * @brief Factory class for ContinuousPhotonSource instances (currently there
 * is only one instance, and the factory is only used to either generate it, or
 * return a nullptr).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CONTINUOUSPHOTONSOURCEFACTORY_HPP
#define CONTINUOUSPHOTONSOURCEFACTORY_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "IsotropicContinuousPhotonSource.hpp"

/**
 * @brief Factory class for ContinuousPhotonSource instances (currently there
 * is only one instance, and the factory is only used to either generate it, or
 * return a nullptr).
 */
class ContinuousPhotonSourceFactory {
public:
  /**
   * @brief Generate an IsotropicContinuousPhotonSource, or return a nullptr.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created IsotropicContinuousPhotonSource
   * instance. Memory management for the pointer needs to be handled by the
   * calling routine.
   */
  inline static IsotropicContinuousPhotonSource *generate(ParameterFile &params,
                                                          Log *log = nullptr) {
    std::string type =
        params.get_value< std::string >("continuousphotonsource.type", "None");
    if (log) {
      log->write_info("Requested ContinuousPhotonSource type: ", type, ".");
    }
    if (type == "Isotropic") {
      return new IsotropicContinuousPhotonSource(params, log);
    } else if (type == "None") {
      return nullptr;
    } else {
      cmac_error("Unknown ContinuousPhotonSource type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // CONTINUOUSPHOTONSOURCEFACTORY_HPP
