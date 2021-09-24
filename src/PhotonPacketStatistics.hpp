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
 * @file PhotonPacketStatistics.hpp
 *
 * @brief Statistical information about the photon packets and re-emission
 * events.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONPACKETSTATISTICS_HPP
#define PHOTONPACKETSTATISTICS_HPP

#include "AtomicValue.hpp"
#include "ParameterFile.hpp"
#include "PhotonPacket.hpp"
#include <fstream>
#include <vector>

/**
 * @brief Statistical information about the photon packets and re-emission
 * events.
 */
class PhotonPacketStatistics {
private:
  /**
   * @brief variable to store the number of photons on bins according to number
   * of scatterings.
   */
  std::vector< AtomicValue< uint_fast64_t > > _scatter_histogram;

public:
  /**
   * @brief constructor
   *
   * @param max_scatter maximum number of scatters recorded by the statistics
   */
  inline PhotonPacketStatistics(uint_fast32_t max_scatter)
      : _scatter_histogram(max_scatter + 2) {}
  /**
   * @brief parameter file constructor
   *
   * These are the parameters that are used by this function:
   *   - maximum number of scatters; (default: 5)
   * @param params reference to the parameter file
   */
  inline PhotonPacketStatistics(ParameterFile &params)
      : PhotonPacketStatistics(params.get_value< uint_fast32_t >(
            "PhotonPacketStatistics:maximum number of scatters", 5)) {}
  /**
   * @brief function that implements photon absorption termination
   *
   * @param packet photon packet retrieved to be terminated
   */
  inline void absorb_photon(const PhotonPacket &packet) {
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
  }
  /**
   * @brief Function that implements photon escape termination
   *
   * @param packet Photon packet.
   */
  inline void escape_photon(const PhotonPacket &packet) {
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
  }
  /**
   * @brief function that outputs re-emission statistics of photons
   */
  inline void print_stats() {
    std::ofstream output_stats("photon_statistics.txt");
    output_stats << "# Scattering statistics for photons\n";
    output_stats << "# Nscatter\t BinCount  \n";
    for (uint_fast32_t i = 0; i < _scatter_histogram.size(); i++) {
      output_stats << i << "\t" << _scatter_histogram[i].value() << "\n";
    }
  }
};

#endif // PHOTONPACKETSTATISTICS_HPP
