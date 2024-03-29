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
 * @file MultiTracker.cpp
 *
 * @brief MultiTracker implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "MultiTracker.hpp"
#include "TrackerFactory.hpp"

/**
 * @brief YAMLDictionary constructor.
 *
 * @param name Name of the block in the dictionary that contains additional
 * parameters for the spectrum tracker.
 * @param blocks YAMLDictionary that contains additional parameters.
 */
MultiTracker::MultiTracker(const std::string name, YAMLDictionary &blocks) {

  const uint_fast32_t number_of_trackers =
      blocks.get_value< uint_fast32_t >(name + "number of trackers");

  _trackers.resize(number_of_trackers, nullptr);
  _output_names.resize(number_of_trackers);

  for (uint_fast32_t i = 0; i < number_of_trackers; ++i) {
    std::stringstream blockname;
    blockname << name << "tracker[" << i << "]:";
    _trackers[i] = TrackerFactory::generate(blockname.str(), blocks);
    _output_names[i] =
        blocks.get_value< std::string >(blockname.str() + "output name", "");
  }
}

/**
 * @brief Virtual destructor.
 */
MultiTracker::~MultiTracker() {
  for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
    delete _trackers[i];
  }
}

/**
 * @brief Make a duplicate of the current tracker.
 *
 * @return Pointer to a new duplicate of the tracker.
 */
Tracker *MultiTracker::duplicate() {
  MultiTracker *copy = new MultiTracker();
  copy->_trackers.resize(_trackers.size(), nullptr);
  copy->_output_names.resize(_trackers.size());
  for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
    copy->_trackers[i] = _trackers[i]->duplicate();
    copy->_output_names[i] = _output_names[i];
  }
  return copy;
}

/**
 * @brief Add the contribution from the given duplicate tracker to this
 * tracker.
 *
 * @param tracker Duplicate tracker (created using Tracker::duplicate()).
 */
void MultiTracker::merge(Tracker *tracker) {
  const MultiTracker *other = reinterpret_cast< MultiTracker * >(tracker);
  for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
    _trackers[i]->merge(other->_trackers[i]);
  }
}

/**
 * @brief Add the contribution of the given photon packet to the bins.
 *
 * @param photon Photon to add.
 */
void MultiTracker::count_photon(const Photon &photon) {

  for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
    _trackers[i]->count_photon(photon);
  }
}

/**
 * @brief Add the contribution of the given photon packet to the bins.
 *
 * @param photon Photon to add.
 * @param absorption Absorption counters within the cell for this photon
 * (in m^-1).
 */
void MultiTracker::count_photon(const PhotonPacket &photon,
                                const double *absorption) {

  for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
    _trackers[i]->count_photon(photon, absorption);
  }
}

/**
 * @brief Output the trackers to the file with the given name.
 *
 * We output one file for every tracker in this tracker + an additional file
 * with the names of all output files.
 *
 * @param filename Name of the output file.
 */
void MultiTracker::output_tracker(const std::string filename) const {

  std::ofstream ofile(filename);
  for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
    std::string this_filename = _output_names[i];
    if (this_filename == "") {
      std::stringstream this_filename_stream;
      this_filename_stream << filename << "." << i << ".txt";
      this_filename = this_filename_stream.str();
    }
    _trackers[i]->output_tracker(this_filename);

    ofile << "tracker[" << i << "]:\n";
    ofile << "  output name: " << this_filename << "\n";
    _trackers[i]->describe("  ", ofile);
  }
}
