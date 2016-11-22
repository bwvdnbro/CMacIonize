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
#include "Abundances.hpp"
#include "CrossSections.hpp"
#include "DensityValues.hpp"
#include "ElementNames.hpp"
#include "Error.hpp"
#include "IsotropicContinuousPhotonSource.hpp"
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
 * @param discrete_spectrum PhotonSourceSpectrum for the discrete photon
 * sources.
 * @param continuous_source IsotropicContinuousPhotonSource.
 * @param continuous_spectrum PhotonSourceSpectrum for the continuous photon
 * source.
 * @param abundances Abundances of the elements in the ISM.
 * @param cross_sections Cross sections for photoionization.
 * @param random_generator RandomGenerator used to generate random numbers.
 * @param log Log to write logging info to.
 */
PhotonSource::PhotonSource(PhotonSourceDistribution *distribution,
                           PhotonSourceSpectrum *discrete_spectrum,
                           IsotropicContinuousPhotonSource *continuous_source,
                           PhotonSourceSpectrum *continuous_spectrum,
                           Abundances &abundances,
                           CrossSections &cross_sections,
                           RandomGenerator &random_generator, Log *log)
    : _discrete_number_of_photons(0), _discrete_spectrum(discrete_spectrum),
      _continuous_source(continuous_source),
      _continuous_spectrum(continuous_spectrum), _abundances(abundances),
      _cross_sections(cross_sections), _random_generator(random_generator),
      _HLyc_spectrum(cross_sections, random_generator),
      _HeLyc_spectrum(cross_sections, random_generator),
      _He2pc_spectrum(random_generator), _log(log) {

  double luminosity_discrete = 0.;
  double luminosity_continuous = 0.;
  if (distribution != nullptr) {
    _discrete_positions.resize(distribution->get_number_of_sources());
    _discrete_weights.resize(distribution->get_number_of_sources());
    for (unsigned int i = 0; i < _discrete_positions.size(); ++i) {
      _discrete_positions[i] = distribution->get_position(i);
      _discrete_weights[i] = distribution->get_weight(i);
    }
    luminosity_discrete = distribution->get_total_luminosity();

    if (_log) {
      _log->write_status("Constructed PhotonSource with ",
                         _discrete_positions.size(), " positions and weights.");
    }

    _discrete_active_source_index = 0;
    _discrete_active_photon_index = 0;
    _discrete_active_number_of_photons = 0;
  }

  if (continuous_source != nullptr) {
    _continuous_active_number_of_photons = 0;
    luminosity_continuous = continuous_source->get_total_luminosity();
  }

  _total_luminosity = luminosity_discrete + luminosity_continuous;
  _discrete_fraction = luminosity_discrete / _total_luminosity;

  if (_log) {
    _log->write_status(_discrete_fraction * 100.,
                       "% of the photons is emitted by discrete sources.");
  }
}

/**
 * @brief Set the number of photons that should be emitted from this source
 * during the next iteration.
 *
 * This also resets the internal counters.
 *
 * @param number_of_photons Number of photons during the next iteration.
 * @return Actual number of photons that was set (since the weights might not
 * sum nicely to 1).
 */
unsigned int
PhotonSource::set_number_of_photons(unsigned int number_of_photons) {
  _discrete_number_of_photons = 0;
  _continuous_number_of_photons = 0;

  // make sure we have at least 10 photons per discrete source
  if (number_of_photons * _discrete_fraction < 10 * _discrete_weights.size()) {
    number_of_photons = 10 * _discrete_weights.size() / _discrete_fraction;
  }

  // make sure we have at least 100 photons for the continuous sources
  if (_discrete_fraction < 1. &&
      number_of_photons * (1. - _discrete_fraction) < 100) {
    number_of_photons = 100 / (1. - _discrete_fraction);
  }

  while (number_of_photons !=
         _discrete_number_of_photons + _continuous_number_of_photons) {
    _discrete_number_of_photons = 0;
    _continuous_number_of_photons =
        std::round((1. - _discrete_fraction) * number_of_photons);
    for (unsigned int i = 0; i < _discrete_weights.size(); ++i) {
      _discrete_number_of_photons += std::round(
          number_of_photons * _discrete_fraction * _discrete_weights[i]);
    }
    number_of_photons =
        _discrete_number_of_photons + _continuous_number_of_photons;
  }

  _discrete_number_of_photons =
      std::round(number_of_photons * _discrete_fraction);

  _discrete_active_source_index = 0;
  _discrete_active_photon_index = 0;
  if (_discrete_weights.size() > 0) {
    _discrete_active_number_of_photons =
        std::round(_discrete_number_of_photons * _discrete_weights[0]);
  } else {
    _discrete_active_number_of_photons = 0;
  }

  _continuous_number_of_photons =
      std::round((1. - _discrete_fraction) * number_of_photons);
  if (_discrete_number_of_photons == 0) {
    _continuous_active_number_of_photons = 0;
  } else {
    // disable continuous sources and first do discrete sources
    _continuous_active_number_of_photons = _continuous_number_of_photons;
  }

  if (_log) {
    _log->write_info("Number of photons for PhotonSource reset to ",
                     _discrete_number_of_photons, " discrete photons and ",
                     _continuous_number_of_photons, " continuous photons.");
  }

  return _discrete_number_of_photons + _continuous_number_of_photons;
}

/**
 * @brief Set the cross sections of the given photon with the given energy.
 *
 * @param photon Photon.
 * @param energy Energy of the photon (in Hz).
 */
void PhotonSource::set_cross_sections(Photon &photon, double energy) {
  for (int i = 0; i < NUMBER_OF_ELEMENTS; ++i) {
    ElementName element = static_cast< ElementName >(i);
    photon.set_cross_section(
        element, _cross_sections.get_cross_section(element, energy));
  }
  photon.set_cross_section_He_corr(_abundances.get_abundance(ATOM_HELIUM) *
                                   photon.get_cross_section(ELEMENT_He));
}

/**
 * @brief Get a photon with a random direction and energy, originating at one
 * of the discrete sources.
 *
 * @return Photon.
 */
Photon PhotonSource::get_random_photon() {

  CoordinateVector<> position, direction;
  double energy;
  // check if we have a continuous or a discrete source photon
  if (_continuous_active_number_of_photons < _continuous_number_of_photons) {
    std::pair< CoordinateVector<>, CoordinateVector<> > posdir =
        _continuous_source->get_random_incoming_direction();
    position = posdir.first;
    direction = posdir.second;
    energy = _continuous_spectrum->get_random_frequency();
  } else {
    position = _discrete_positions[_discrete_active_source_index];
    direction = get_random_direction();
    energy = _discrete_spectrum->get_random_frequency();
  }

  Photon photon(position, direction, energy);
  set_cross_sections(photon, energy);

  update_indices();

  return photon;
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
 * @return True if the photon is re-emitted as an ionizing photon, false if it
 * leaves the system.
 */
bool PhotonSource::reemit(Photon &photon, DensityValues &cell) {
  double new_frequency = 0.;
  double helium_abundance = _abundances.get_abundance(ATOM_HELIUM);
  double pHabs = 1. / (1. +
                       cell.get_ionic_fraction(ELEMENT_He) * helium_abundance *
                           photon.get_cross_section(ELEMENT_He) /
                           cell.get_ionic_fraction(ELEMENT_H) /
                           photon.get_cross_section(ELEMENT_H));

  double x = _random_generator.get_uniform_random_double();
  if (x <= pHabs) {
    // photon absorbed by hydrogen
    x = _random_generator.get_uniform_random_double();
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
    x = _random_generator.get_uniform_random_double();
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
      x = _random_generator.get_uniform_random_double();
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
                           77. * cell.get_ionic_fraction(ELEMENT_He) /
                               sqrt(cell.get_temperature()) /
                               cell.get_ionic_fraction(ELEMENT_H));
      x = _random_generator.get_uniform_random_double();
      if (x < pHots) {
        // absorbed on the spot
        x = _random_generator.get_uniform_random_double();
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
        x = _random_generator.get_uniform_random_double();
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

  set_cross_sections(photon, new_frequency);

  return true;
}
