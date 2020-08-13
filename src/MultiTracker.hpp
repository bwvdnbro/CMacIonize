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
 * @file MultiTracker.hpp
 *
 * @brief Tracker that allows attaching multiple trackers to one cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MULTITRACKER_HPP
#define MULTITRACKER_HPP

#include "Tracker.hpp"

#include <sstream>
#include <vector>

class Photon;
class YAMLDictionary;

/**
 * @brief Tracker that allows attaching multiple trackers to one cell.
 */
class MultiTracker : public Tracker {
private:
  /*! @brief Trackers attached to this tracker. */
  std::vector< Tracker * > _trackers;

  /*! @brief Output file names for the different trackers attached to this
   *  tracker. */
  std::vector< std::string > _output_names;

public:
  /**
   * @brief Empty constructor.
   */
  MultiTracker() {}

  MultiTracker(const std::string name, YAMLDictionary &blocks);

  virtual ~MultiTracker();

  /**
   * @brief Add the given tracker to the list of stored trackers.
   *
   * @param tracker Tracker to add.
   * @param name Name for this tracker (default: tracker X, where X is the index
   * of this tracker in the list).
   */
  inline void add_tracker(Tracker *tracker, std::string name = "") {

    if (name == "") {
      std::stringstream namestr;
      namestr << "tracker " << _trackers.size();
      name = namestr.str();
    }
    _trackers.push_back(tracker);
    _output_names.push_back(name);
  }

  /**
   * @brief Erase the internally stored trackers to make sure they are not
   * deleted twice.
   */
  inline void erase_trackers() {
    for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
      _trackers[i] = nullptr;
    }
  }

  virtual Tracker *duplicate();
  virtual void merge(Tracker *tracker);

  virtual void count_photon(const Photon &photon);
  virtual void count_photon(const PhotonPacket &photon,
                            const double *absorption);

  virtual void output_tracker(const std::string filename) const;

  /**
   * @brief Is this tracker a MultiTracker?
   *
   * @return True.
   */
  virtual bool is_multi_tracker() const { return true; }
};

#endif // MULTITRACKER_HPP
