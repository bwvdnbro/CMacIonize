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
 * @file Tracker.hpp
 *
 * @brief General interface for trackers that record photon properties.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TRACKER_HPP
#define TRACKER_HPP

#include "Photon.hpp"

#include <string>

/**
 * @brief General interface for trackers that record photon properties.
 */
class Tracker {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~Tracker() {}

  /**
   * @brief Add the contribution of the given photon packet to the tracker.
   *
   * @param photon Photon to add.
   */
  virtual void count_photon(const Photon &photon) = 0;

  /**
   * @brief Output the tracker data to the file with the given name.
   *
   * @param filename Name of the output file.
   */
  virtual void output_tracker(const std::string filename) const = 0;

  /**
   * @brief Describe the tracker in the given output stream, appending the given
   * prefix to each line of output.
   *
   * @param prefix Prefix to add to each output line.
   * @param stream std::ostream to write to.
   */
  virtual void describe(const std::string prefix, std::ostream &stream) const {}
};

#endif // TRACKER_HPP
