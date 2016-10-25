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
#include "DensityValues.hpp"
#include "ElementNames.hpp"
#include "Error.hpp"
#include "Log.hpp"
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
 * @param log Log to write logging info to.
 */
PhotonSource::PhotonSource(PhotonSourceDistribution &distribution,
                           PhotonSourceSpectrum &spectrum,
                           CrossSections &cross_sections, Log *log)
    : _number_of_photons(0), _spectrum(spectrum),
      _cross_sections(cross_sections), _HLyc_spectrum(cross_sections),
      _HeLyc_spectrum(cross_sections), _log(log) {
  _positions.resize(distribution.get_number_of_sources());
  _weights.resize(distribution.get_number_of_sources());
  for (unsigned int i = 0; i < _positions.size(); ++i) {
    _positions[i] = distribution.get_position(i);
    _weights[i] = distribution.get_weight(i);
  }
  _total_luminosity = distribution.get_total_luminosity();

  if (_log) {
    _log->write_status("Constructed PhotonSource with ", _positions.size(),
                       " positions and weights.");
  }

  _active_source_index = 0;
  _active_photon_index = 0;
  _active_number_of_photons = 0;
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
  _active_number_of_photons = round(_number_of_photons * _weights[0]);

  if (_log) {
    _log->write_status("Number of photons for PhotonSource reset to ",
                       _number_of_photons, ".");
  }
}

/**
 * @brief Get a photon with a random direction and energy, originating at one
 * of the discrete sources.
 *
 * @return Photon.
 */
Photon PhotonSource::get_random_photon() {
  if (_active_source_index == _positions.size()) {
    if (_log) {
      _log->write_warning("Too many photons requested!");
      --_active_source_index;
    }
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
          round(_number_of_photons * _weights[_active_source_index]);
    }
  }

  return Photon(position, direction, energy, xsecH, xsecHe);
}

/**
 * @brief Get the total luminosity of all sources together.
 *
 * @return Total luminosity (in s^-1).
 */
double PhotonSource::get_total_luminosity() { return _total_luminosity; }

/**
 * @brief Reemit the given Photon.
 *
 * This routine randomly chooses if the photon is absorbed by hydrogen or
 * helium, and then reemits it at a new random frequency and in a new random
 * direction.
 *
 * @param photon Photon to reemit.
 * @param cell DensityValues of the cell in which the Photon currently resides.
 * @param helium_abundance Helium abundance.
 * @return True if the photon is re-emitted as an ionizing photon, false if it
 * leaves the system.
 */
bool PhotonSource::reemit(Photon &photon, DensityValues &cell,
                          double helium_abundance) {
  double new_frequency = 0.;
  double pHabs = 1. / (1. +
                       cell.get_neutral_fraction_He() * helium_abundance *
                           photon.get_helium_cross_section() /
                           cell.get_neutral_fraction_H() /
                           photon.get_hydrogen_cross_section());

  double x = Utilities::random_double();
  if (x <= pHabs) {
    // photon absorbed by hydrogen
    x = Utilities::random_double();
    if (x <= cell.get_pHion()) {
      // sample new frequency from H Ly c
      _HLyc_spectrum.set_temperature(cell.get_temperature());
      new_frequency = _HLyc_spectrum.get_random_frequency();
      photon.set_type(PHOTONTYPE_DIFFUSE_HI);
    } else {
      // photon escapes
      photon.set_type(PHOTONTYPE_ABSORBED);
      return false;
    }
  } else {
    // photon absorbed by helium
    x = Utilities::random_double();
    if (x <= cell.get_pHe_em(0)) {
      // sample new frequency from He Ly c
      _HeLyc_spectrum.set_temperature(cell.get_temperature());
      new_frequency = _HeLyc_spectrum.get_random_frequency();
      photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
    } else if (x <= cell.get_pHe_em(1)) {
      // new frequency is 19.8eV (no idea why)
      new_frequency = 19.8 / 13.6;
      photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
    } else if (x <= cell.get_pHe_em(2)) {
      x = Utilities::random_double();
      if (x < 0.56) {
        // sample new frequency from H-ionizing part of He 2-photon continuum
        new_frequency = _He2pc_spectrum.get_random_frequency();
        photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
      } else {
        // photon escapes
        photon.set_type(PHOTONTYPE_ABSORBED);
        return false;
      }
    } else if (x <= cell.get_pHe_em(3)) {
      // HeI Ly-alpha, is either absorbed on the spot or converted to HeI
      // 2-photon continuum
      double pHots = 1. / (1. +
                           77. * cell.get_neutral_fraction_He() /
                               sqrt(cell.get_temperature()) /
                               cell.get_neutral_fraction_H());
      x = Utilities::random_double();
      if (x < pHots) {
        // absorbed on the spot
        x = Utilities::random_double();
        if (x <= cell.get_pHion()) {
          // H Ly c, like above
          _HLyc_spectrum.set_temperature(cell.get_temperature());
          new_frequency = _HLyc_spectrum.get_random_frequency();
          photon.set_type(PHOTONTYPE_DIFFUSE_HI);
        } else {
          // photon escapes
          photon.set_type(PHOTONTYPE_ABSORBED);
          return false;
        }
      } else {
        // He 2-photon continuum
        x = Utilities::random_double();
        if (x < 0.56) {
          // sample like above
          new_frequency = _He2pc_spectrum.get_random_frequency();
          photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
        } else {
          // photon escapes
          photon.set_type(PHOTONTYPE_ABSORBED);
          return false;
        }
      }
    } else {
      // not in Kenny's code, since the probabilities above are forced to sum
      // to 1.
      // the code below is hence never executed
      photon.set_type(PHOTONTYPE_ABSORBED);
      return false;
    }
  }

  photon.set_energy(new_frequency);

  CoordinateVector<> direction = get_random_direction();
  photon.set_direction(direction);

  double xsecH = _cross_sections.get_cross_section(ELEMENT_H, new_frequency);
  double xsecHe = _cross_sections.get_cross_section(ELEMENT_He, new_frequency);
  photon.set_hydrogen_cross_section(xsecH);
  photon.set_helium_cross_section(xsecHe);

  return true;
}
