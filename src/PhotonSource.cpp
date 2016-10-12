/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file PhotonSource.cpp
 *
 * @brief Photon source: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "PhotonSource.hpp"
#include "CrossSections.hpp"
#include "ElementNames.hpp"
#include "Error.hpp"
#include "PhotonSourceDistribution.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "Utilities.hpp"
#include <cmath>
using namespace std;

/**
 * @brief Constructor.
 *
 * @param distribution PhotonSourceDistribution giving the positions of the
 * discrete photon sources.
 * @param spectrum PhotonSourceSpectrum for the discrete photon sources.
 * @param cross_sections Cross sections for photoionization.
 * @param number_of_photons Total number of photons emitted from all discrete
 * photon sources together.
 */
PhotonSource::PhotonSource(PhotonSourceDistribution &distribution,
                           PhotonSourceSpectrum &spectrum,
                           CrossSections &cross_sections,
                           unsigned int number_of_photons)
    : _number_of_photons(number_of_photons), _spectrum(spectrum),
      _cross_sections(cross_sections) {
  _positions.resize(distribution.get_number_of_sources());
  _weights.resize(distribution.get_number_of_sources());
  for (unsigned int i = 0; i < _positions.size(); ++i) {
    _positions[i] = distribution.get_position(i);
    _weights[i] = distribution.get_weight(i);
  }

  _active_source_index = 0;
  _active_photon_index = 0;
  _active_number_of_photons = _number_of_photons * _weights[0];
}

/**
 * @brief Set the number of photons that should be emitted from this source
 * during the next iteration.
 *
 * This also resets the internal counters.
 *
 * @param number_of_photons Number of photons during the next iteration.
 */
void PhotonSource::set_number_of_photons(unsigned int number_of_photons) {
  _number_of_photons = number_of_photons;
  _active_source_index = 0;
  _active_photon_index = 0;
  _active_number_of_photons = _number_of_photons * _weights[0];
}

/**
 * @brief Get a photon with a random direction and energy, originating at one
 * of the discrete sources.
 *
 * @return Photon.
 */
Photon PhotonSource::get_random_photon() {
  if (_active_source_index == _positions.size()) {
    error("No more photons available!");
  }

  CoordinateVector<> position = _positions[_active_source_index];

  CoordinateVector<> direction = get_random_direction();

  double energy = _spectrum.get_random_frequency();
  double xsecH = _cross_sections.get_cross_section(ELEMENT_H, energy);
  double xsecHe = _cross_sections.get_cross_section(ELEMENT_He, energy);

  ++_active_photon_index;
  if (_active_photon_index == _active_number_of_photons) {
    _active_photon_index = 0;
    ++_active_source_index;
    if (_active_source_index < _positions.size()) {
      _active_number_of_photons =
          _number_of_photons * _weights[_active_source_index];
    }
  }

  return Photon(position, direction, energy, xsecH, xsecHe);
}
