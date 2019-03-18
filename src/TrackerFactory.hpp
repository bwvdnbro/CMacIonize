/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file TrackerFactory.hpp
 *
 * @brief Factory for Tracker instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TRACKERFACTORY_HPP
#define TRACKERFACTORY_HPP

#include "Tracker.hpp"
#include "YAMLDictionary.hpp"

// implementations
#include "SpectrumTracker.hpp"

/**
 * @brief Factory for Tracker instances.
 */
class TrackerFactory {
public:
  /**
   * @brief Generate a Tracker instance of the type found in the given
   * YAMLDictionary.
   *
   * Supported types are:
   *  - Spectrum: Spectrum tracker
   *
   * @param name Name of the corresponding block in the dictionary.
   * @param blocks YAMLDictionary to read from.
   * @return Pointer to a newly created Tracker instance. Memory management for
   * this pointer should be done by the calling routine.
   */
  inline static Tracker *generate(const std::string name,
                                  YAMLDictionary &blocks) {

    const std::string type = blocks.get_value< std::string >(name + "type");

    if (type == "Spectrum") {
      return new SpectrumTracker(name, blocks);
    } else {
      cmac_error("Unknown Tracker type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

#endif // TRACKERFACTORY_HPP
