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
  MultiTracker(const std::string name, YAMLDictionary &blocks);

  virtual ~MultiTracker();

  virtual void count_photon(const Photon &photon);

  virtual void output_tracker(const std::string filename) const;
};

#endif // MULTITRACKER_HPP
