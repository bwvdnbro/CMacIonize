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
 * @file FrequencyBinsFactory.hpp
 *
 * @brief Factory for FrequencyBins instances.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FREQUENCYBINSFACTORY_HPP
#define FREQUENCYBINSFACTORY_HPP

#include "FrequencyBins.hpp"
#include "YAMLDictionary.hpp"

// implementations
#include "LevelFrequencyBins.hpp"
#include "LinearFrequencyBins.hpp"

/**
 * @brief Factory for FrequencyBins instances.
 */
class FrequencyBinsFactory {
public:
  /**
   * @brief Generate a FrequencyBins instance based on the parameters in the
   * given YAMLDictionary.
   *
   * Supported types are (default: Linear):
   *  - Level: The bottom frequency for each bin corresponds to the ionization
   *    energy for one of the ions in the simulation.
   *  - Linear: Fixed number of linear bins between a given lower and upper
   *    limit.
   *
   * @param name Name of the block in the dictionary that contains the
   * parameters for this FrequencyBins object.
   * @param blocks YAMLDictionary to read from.
   * @return Pointer to a newly created FrequencyBins instance. Memory
   * management for this pointer should be done by the calling routine.
   */
  inline static FrequencyBins *generate(const std::string name,
                                        YAMLDictionary &blocks) {

    const std::string type =
        blocks.get_value< std::string >(name + "FrequencyBins:type", "Linear");

    if (type == "Level") {
      return new LevelFrequencyBins();
    } else if (type == "Linear") {
      return new LinearFrequencyBins(name, blocks);
    } else {
      cmac_error("Unknown FrequencyBins type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // FREQUENCYBINSFACTORY_HPP
