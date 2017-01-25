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
#include "ContinuousPhotonSource.hpp"
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
 * @param discrete_spectrum PhotonSourceSpectrum for the discrete photon
 * sources.
 * @param continuous_source IsotropicContinuousPhotonSource.
 * @param continuous_spectrum PhotonSourceSpectrum for the continuous photon
 * source.
 * @param abundances Abundances of the elements in the ISM.
 * @param cross_sections Cross sections for photoionization.
 * @param log Log to write logging info to.
 */
PhotonSource::PhotonSource(PhotonSourceDistribution *distribution,
                           PhotonSourceSpectrum *discrete_spectrum,
                           ContinuousPhotonSource *continuous_source,
                           PhotonSourceSpectrum *continuous_spectrum,
                           Abundances &abundances,
                           CrossSections &cross_sections, Log *log)
    : _discrete_spectrum(discrete_spectrum),
      _continuous_source(continuous_source),
      _continuous_spectrum(continuous_spectrum), _abundances(abundances),
      _cross_sections(cross_sections), _HLyc_spectrum(cross_sections),
      _HeLyc_spectrum(cross_sections), _log(log) {

  double discrete_luminosity = 0.;
  double continuous_luminosity = 0.;
  if (distribution != nullptr) {
    _discrete_positions.resize(distribution->get_number_of_sources());
    _discrete_probabilities.resize(distribution->get_number_of_sources());
    for (unsigned int i = 0; i < _discrete_positions.size(); ++i) {
      _discrete_positions[i] = distribution->get_position(i);
      if (i > 0) {
        _discrete_probabilities[i] =
            _discrete_probabilities[i - 1] + distribution->get_weight(i);
      } else {
        _discrete_probabilities[i] = distribution->get_weight(i);
      }
    }
    if (_discrete_probabilities.size() > 0) {
      if (std::abs(_discrete_probabilities.back() - 1.) > 1.e-9) {
        cmac_error("Discrete source weights do not sum to 1.0 (%g)!",
                   _discrete_probabilities.back());
      } else {
        _discrete_probabilities.back() = 1.;
      }
    }
    discrete_luminosity = distribution->get_total_luminosity();

    if (_log) {
      _log->write_status("Constructed PhotonSource with ",
                         _discrete_positions.size(), " positions and weights.");
    }
  }

  if (continuous_source != nullptr) {
    continuous_luminosity = continuous_source->get_total_surface_area() *
                            continuous_spectrum->get_total_flux();
  }

  _total_luminosity = discrete_luminosity + continuous_luminosity;

  if (_total_luminosity == 0.) {
    if (_log) {
      _log->write_error("Total luminosity of all sources is zero! Not doing "
                        "radiative transfer, as there is no radiation to "
                        "propagate.");
      cmac_error("Total luminosity is zero!");
    }
  }

  double discrete_fraction = discrete_luminosity / _total_luminosity;

  if (discrete_luminosity > 0.) {
    if (continuous_luminosity > 0.) {
      _continuous_probability = 0.5;
    } else {
      _continuous_probability = 0.;
    }
    _discrete_photon_weight = 1.;
    _continuous_photon_weight = (1. - _continuous_probability) *
                                continuous_luminosity /
                                _continuous_probability / discrete_luminosity;
  } else {
    _continuous_probability = 1.;
    _discrete_photon_weight = 0.;
    _continuous_photon_weight = 1.;
  }

  if (_log) {
    _log->write_status("Total luminosity of discrete sources: ",
                       discrete_luminosity, " s^-1.");
    _log->write_status("Total luminosity of continuous sources: ",
                       continuous_luminosity, " s^-1.");
    _log->write_status(
        discrete_fraction * 100.,
        "% of the ionizing radiation is emitted by discrete sources.");
  }
}

/**
 * @brief Set the cross sections of the given photon with the given energy.
 *
 * @param photon Photon.
 * @param energy Energy of the photon (in Hz).
 */
void PhotonSource::set_cross_sections(Photon &photon, double energy) const {
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    IonName ion = static_cast< IonName >(i);
    photon.set_cross_section(ion,
                             _cross_sections.get_cross_section(ion, energy));
  }
  photon.set_cross_section_He_corr(_abundances.get_abundance(ELEMENT_He) *
                                   photon.get_cross_section(ION_He_n));
}

/**
 * @brief Get a photon with a random direction and energy, originating at one
 * of the discrete sources.
 *
 * @param random_generator RandomGenerator to use.
 * @return Photon.
 */
Photon
PhotonSource::get_random_photon(RandomGenerator &random_generator) const {

  CoordinateVector<> position, direction;
  double energy;
  double weight;

  double x = random_generator.get_uniform_random_double();
  if (x > _continuous_probability) {
    // discrete photon
    x = random_generator.get_uniform_random_double();
    unsigned int i = 0;
    while (x > _discrete_probabilities[i]) {
      ++i;
    }
    position = _discrete_positions[i];
    direction = get_random_direction(random_generator);
    energy = _discrete_spectrum->get_random_frequency(random_generator);
    weight = _discrete_photon_weight;
  } else {
    // continuous photon
    std::pair< CoordinateVector<>, CoordinateVector<> > posdir =
        _continuous_source->get_random_incoming_direction(random_generator);
    position = posdir.first;
    direction = posdir.second;
    energy = _continuous_spectrum->get_random_frequency(random_generator);
    weight = _continuous_photon_weight;
  }

  Photon photon(position, direction, energy);
  set_cross_sections(photon, energy);

  photon.set_weight(weight);

  return photon;
}

/**
 * @brief Get the total luminosity of all sources together.
 *
 * @return Total luminosity (in s^-1).
 */
double PhotonSource::get_total_luminosity() const { return _total_luminosity; }

/**
 * @brief Reemit the given Photon.
 *
 * This routine randomly chooses if the photon is absorbed by hydrogen or
 * helium, and then reemits it at a new random frequency and in a new random
 * direction.
 *
 * @param photon Photon to reemit.
 * @param cell DensityValues of the cell in which the Photon currently resides.
 * @param random_generator RandomGenerator to use.
 * @return True if the photon is re-emitted as an ionizing photon, false if it
 * leaves the system.
 */
bool PhotonSource::reemit(Photon &photon, DensityValues &cell,
                          RandomGenerator &random_generator) const {
  double new_frequency = 0.;
  double helium_abundance = _abundances.get_abundance(ELEMENT_He);
  double pHabs = 1. / (1. +
                       cell.get_ionic_fraction(ION_He_n) * helium_abundance *
                           photon.get_cross_section(ION_He_n) /
                           cell.get_ionic_fraction(ION_H_n) /
                           photon.get_cross_section(ION_H_n));

  double x = random_generator.get_uniform_random_double();
  if (x <= pHabs) {
    // photon absorbed by hydrogen
    x = random_generator.get_uniform_random_double();
    if (x <= cell.get_pHion()) {
      // sample new frequency from H Ly c
      new_frequency = _HLyc_spectrum.get_random_frequency(
          random_generator, cell.get_temperature());
      photon.set_type(PHOTONTYPE_DIFFUSE_HI);
    } else {
      // photon escapes
      photon.set_type(PHOTONTYPE_ABSORBED);
      return false;
    }
  } else {
    // photon absorbed by helium
    x = random_generator.get_uniform_random_double();
    if (x <= cell.get_pHe_em(0)) {
      // sample new frequency from He Ly c
      new_frequency = _HeLyc_spectrum.get_random_frequency(
          random_generator, cell.get_temperature());
      photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
    } else if (x <= cell.get_pHe_em(1)) {
      // new frequency is 19.8eV
      new_frequency = 4.788e15;
      photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
    } else if (x <= cell.get_pHe_em(2)) {
      x = random_generator.get_uniform_random_double();
      if (x < 0.56) {
        // sample new frequency from H-ionizing part of He 2-photon continuum
        new_frequency = _He2pc_spectrum.get_random_frequency(
            random_generator, cell.get_temperature());
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
                           77. * cell.get_ionic_fraction(ION_He_n) /
                               sqrt(cell.get_temperature()) /
                               cell.get_ionic_fraction(ION_H_n));
      x = random_generator.get_uniform_random_double();
      if (x < pHots) {
        // absorbed on the spot
        x = random_generator.get_uniform_random_double();
        if (x <= cell.get_pHion()) {
          // H Ly c, like above
          new_frequency = _HLyc_spectrum.get_random_frequency(
              random_generator, cell.get_temperature());
          photon.set_type(PHOTONTYPE_DIFFUSE_HI);
        } else {
          // photon escapes
          photon.set_type(PHOTONTYPE_ABSORBED);
          return false;
        }
      } else {
        // He 2-photon continuum
        x = random_generator.get_uniform_random_double();
        if (x < 0.56) {
          // sample like above
          new_frequency =
              _He2pc_spectrum.get_random_frequency(random_generator);
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

  CoordinateVector<> direction = get_random_direction(random_generator);
  photon.set_direction(direction);

  set_cross_sections(photon, new_frequency);

  return true;
}
