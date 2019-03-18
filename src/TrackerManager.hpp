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
 * @file TrackerManager.hpp
 *
 * @brief Class that manages Trackers in the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TRACKERMANAGER_HPP
#define TRACKERMANAGER_HPP

#include "DensityGrid.hpp"
#include "ParameterFile.hpp"
#include "TrackerFactory.hpp"
#include "YAMLDictionary.hpp"

#include <fstream>
#include <vector>

/**
 * @brief Class that manages Trackers in the grid.
 */
class TrackerManager {
private:
  /*! @brief List of tracker positions. */
  std::vector< CoordinateVector<> > _tracker_positions;

  /*! @brief List of trackers. */
  std::vector< Tracker * > _trackers;

  /*! @brief List of output file names. */
  std::vector< std::string > _output_names;

  /*! @brief Number of photon packets to use during the tracking step. */
  const uint_fast64_t _number_of_photons;

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the file that contains the positions of the
   * trackers.
   * @param number_of_photons Number of photon packets to use during the
   * tracking step.
   */
  TrackerManager(const std::string filename,
                 const uint_fast64_t number_of_photons)
      : _number_of_photons(number_of_photons) {

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

      _trackers[i] = TrackerFactory::generate(blockname.str(), blocks);

      std::stringstream default_name;
      default_name << "SpectrumTracker" << i << ".txt";
      _output_names[i] = blocks.get_value< std::string >(
          blockname.str() + "output name", default_name.str());
    }

    std::ofstream ofile(filename + ".used-values");
    blocks.print_contents(ofile, true);
    ofile.close();
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read from the parameter file:
   *  - filename: Name of the file that contains the tracker positions
   *    (required)
   *  - minimum number of photon packets: Minimum number of photon packets to
   *    use during the spectrum tracking step (default: 0, meaning we do not
   *    use a different number for the spectrum tracking step)
   *
   * @param params ParameterFile to read from.
   */
  TrackerManager(ParameterFile &params)
      : TrackerManager(
            params.get_filename("SpectrumTrackerManager:filename"),
            params.get_value< uint_fast64_t >(
                "SpectrumTrackerManager:minimum number of photon packets", 0)) {
  }

  /**
   * @brief Destructor.
   */
  ~TrackerManager() {
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
      _trackers[i]->output_tracker(_output_names[i]);
    }
  }

  /**
   * @brief Get the number of photon packets to use during the tracking step.
   *
   * @return Number of photons to use during the tracking step.
   */
  inline uint_fast64_t get_number_of_photons() const {
    return _number_of_photons;
  }
};

#endif // TRACKERMANAGER_HPP
