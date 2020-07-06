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

#include "Configuration.hpp"

#ifdef HAVE_HDF5
#include "HDF5Tools.hpp"
#endif

#include <string>

class Cell;
class Photon;
class PhotonPacket;

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
   * @brief Normalize the tracker for use in a cell with the given size.
   *
   * @param cell Cell the tracker is attached to.
   */
  virtual void normalize_for_cell(const Cell &cell) {}

  /**
   * @brief Normalize the tracker with the appropriate ionizing luminosity per
   * unit photon packet weight.
   *
   * @param luminosity_per_weight Ionizing luminosity per unit photon packet
   * weight (in s^-1).
   */
  virtual void normalize(const double luminosity_per_weight) {}

  /**
   * @brief Make a duplicate of the current tracker.
   *
   * @return Pointer to a new duplicate of the tracker.
   */
  virtual Tracker *duplicate() = 0;

  /**
   * @brief Add the contribution from the given duplicate tracker to this
   * tracker.
   *
   * @param tracker Duplicate tracker (created using Tracker::duplicate()).
   */
  virtual void merge(Tracker *tracker) = 0;

  /**
   * @brief Add the contribution of the given photon packet to the tracker.
   *
   * @param photon Photon to add.
   */
  virtual void count_photon(const Photon &photon) = 0;

  /**
   * @brief Add the contribution of the given photon packet to the tracker.
   *
   * @param photon Photon to add.
   * @param absorption Absorption counters within the cell for this photon
   * (in m^-1).
   */
  virtual void count_photon(const PhotonPacket &photon,
                            const double *absorption) = 0;

  /**
   * @brief Output the tracker data to the file with the given name.
   *
   * @param filename Name of the output file.
   */
  virtual void output_tracker(const std::string filename) const = 0;

#ifdef HAVE_HDF5
  /**
   * @brief Does the given tracker belong to the same group as this tracker?
   *
   * @param tracker Other tracker.
   * @return True if both trackers belong to the same group.
   */
  virtual bool same_group(const Tracker *tracker) const { return false; }

  /**
   * @brief Create the header and shared datasets for an HDF5 group containing
   * one or multiple trackers of this type.
   *
   * @param group HDF5Group to write to.
   * @param group_size Number of trackers in the group.
   */
  virtual void create_group(const HDF5Tools::HDF5Group group,
                            const uint_fast32_t group_size) {
    cmac_error("Function has not been implemented for this tracker type!");
  }

  /**
   * @brief Append the tracker to the given group.
   *
   * We assume the group was already properly initialised using create_group().
   *
   * @param group HDF5Group to write to.
   * @param group_index Index of this particular tracker within the group.
   */
  virtual void append_to_group(const HDF5Tools::HDF5Group group,
                               const uint_fast32_t group_index) {
    cmac_error("Function has not been implemented for this tracker type!");
  }
#endif

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
