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
 * @file SpectrumTrackerManager.hpp
 *
 * @brief Class that manages SpectrumTrackers in the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPECTRUMTRACKERMANAGER_HPP
#define SPECTRUMTRACKERMANAGER_HPP

#include "DensityGrid.hpp"
#include "ParameterFile.hpp"
#include "SpectrumTracker.hpp"
#include "YAMLDictionary.hpp"

#include <fstream>
#include <vector>

/**
 * @brief Class that manages SpectrumTrackers in the grid.
 */
class SpectrumTrackerManager {
private:
  /*! @brief List of tracker positions. */
  std::vector< CoordinateVector<> > _tracker_positions;

  /*! @brief List of trackers. */
  std::vector< SpectrumTracker * > _trackers;

  /*! @brief List of output file names. */
  std::vector< std::string > _output_names;

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the file that contains the positions of the
   * trackers.
   */
  SpectrumTrackerManager(const std::string filename) {

    std::ifstream file(filename);

    if (!file) {
      cmac_error("Error while opening file \"%s\"!", filename.c_str());
    }

    YAMLDictionary blocks(file);
    const uint_fast32_t number_of_trackers =
        blocks.get_value< uint_fast32_t >("number of trackers");
    _tracker_positions.resize(number_of_trackers);
    _trackers.resize(number_of_trackers, nullptr);
    _output_names.resize(number_of_trackers);
    for (uint_fast32_t i = 0; i < number_of_trackers; ++i) {
      std::stringstream blockname;
      blockname << "tracker[" << i << "]:";

      _tracker_positions[i] = blocks.get_physical_vector< QUANTITY_LENGTH >(
          blockname.str() + "position");

      const uint_fast32_t nbin = blocks.get_value< uint_fast32_t >(
          blockname.str() + "number of bins", 100);
      _trackers[i] = new SpectrumTracker(nbin);

      std::stringstream default_name;
      default_name << "SpectrumTracker" << i << ".txt";
      _output_names[i] = blocks.get_value< std::string >(
          blockname.str() + "output name", default_name.str());
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  SpectrumTrackerManager(ParameterFile &params)
      : SpectrumTrackerManager(params.get_value< std::string >(
            "SpectrumTrackerManager:filename")) {}

  /**
   * @brief Destructor.
   */
  ~SpectrumTrackerManager() {
    for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
      delete _trackers[i];
    }
  }

  /**
   * @brief Add the trackers to the given DensityGrid.
   *
   * @param grid Grid to add trackers to.
   */
  inline void add_trackers(DensityGrid &grid) {
    for (uint_fast32_t i = 0; i < _tracker_positions.size(); ++i) {
      grid.get_cell(_tracker_positions[i])
          .get_ionization_variables()
          .add_tracker(_trackers[i]);
    }
  }

  /**
   * @brief Output the tracker information.
   */
  inline void output_trackers() const {
    for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
      _trackers[i]->output_spectrum(_output_names[i]);
    }
  }
};

#endif // SPECTRUMTRACKERMANAGER_HPP
