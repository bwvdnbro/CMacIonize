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
      _HeLyc_spectrum(cross_sections), _scattering_hgg(0.44),
      _scattering_g2(_scattering_hgg * _scattering_hgg),
      _scattering_omg2(1. - _scattering_g2),
      _scattering_thgg(2. * _scattering_hgg),
      _scattering_omhgg(1. - _scattering_hgg),
      _scattering_od2hgg(0.5 / _scattering_hgg),
      _scattering_opg2(1. + _scattering_g2), _scattering_pl(0.43),
      _scattering_sc(1.), _scattering_pc(0.), _log(log) {

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
  if (x >= _continuous_probability) {
    cmac_assert(_discrete_probabilities.size() > 0);
    cmac_assert(_discrete_probabilities.back() == 1.);
    // discrete photon
    x = random_generator.get_uniform_random_double();
    unsigned int i = 0;
    while (x > _discrete_probabilities[i]) {
      ++i;
      cmac_assert(i < _discrete_probabilities.size());
    }
    position = _discrete_positions[i];
    direction = get_random_direction(random_generator);
    energy = _discrete_spectrum->get_random_frequency(random_generator);
    weight = _discrete_photon_weight;
  } else {
    cmac_assert(_continuous_source != nullptr);
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
 * @param ionization_variables IonizationVariables of the cell that contains the
 * current location of the Photon.
 * @param random_generator RandomGenerator to use.
 * @return True if the photon is re-emitted as an ionizing photon, false if it
 * leaves the system.
 */
bool PhotonSource::reemit(Photon &photon,
                          const IonizationVariables &ionization_variables,
                          RandomGenerator &random_generator) const {

  double new_frequency = 0.;
  double helium_abundance = _abundances.get_abundance(ELEMENT_He);
  double pHabs =
      1. / (1. +
            ionization_variables.get_ionic_fraction(ION_He_n) *
                helium_abundance * photon.get_cross_section(ION_He_n) /
                ionization_variables.get_ionic_fraction(ION_H_n) /
                photon.get_cross_section(ION_H_n));

  double x = random_generator.get_uniform_random_double();
  if (x <= pHabs) {
    // photon absorbed by hydrogen
    x = random_generator.get_uniform_random_double();
    if (x <= ionization_variables.get_reemission_probability(
                 REEMISSIONPROBABILITY_HYDROGEN)) {
      // sample new frequency from H Ly c
      new_frequency = _HLyc_spectrum.get_random_frequency(
          random_generator, ionization_variables.get_temperature());
      photon.set_type(PHOTONTYPE_DIFFUSE_HI);
    } else {
      // photon escapes
      photon.set_type(PHOTONTYPE_ABSORBED);
      return false;
    }
  } else {
    // photon absorbed by helium
    x = random_generator.get_uniform_random_double();
    if (x <= ionization_variables.get_reemission_probability(
                 REEMISSIONPROBABILITY_HELIUM_LYC)) {
      // sample new frequency from He Ly c
      new_frequency = _HeLyc_spectrum.get_random_frequency(
          random_generator, ionization_variables.get_temperature());
      photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
    } else if (x <= ionization_variables.get_reemission_probability(
                        REEMISSIONPROBABILITY_HELIUM_NPEEV)) {
      // new frequency is 19.8eV
      new_frequency = 4.788e15;
      photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
    } else if (x <= ionization_variables.get_reemission_probability(
                        REEMISSIONPROBABILITY_HELIUM_TPC)) {
      x = random_generator.get_uniform_random_double();
      if (x < 0.56) {
        // sample new frequency from H-ionizing part of He 2-photon continuum
        new_frequency = _He2pc_spectrum.get_random_frequency(
            random_generator, ionization_variables.get_temperature());
        photon.set_type(PHOTONTYPE_DIFFUSE_HeI);
      } else {
        // photon escapes
        photon.set_type(PHOTONTYPE_ABSORBED);
        return false;
      }
    } else if (x <= ionization_variables.get_reemission_probability(
                        REEMISSIONPROBABILITY_HELIUM_LYA)) {
      // HeI Ly-alpha, is either absorbed on the spot or converted to HeI
      // 2-photon continuum
      double pHots =
          1. / (1. +
                77. * ionization_variables.get_ionic_fraction(ION_He_n) /
                    sqrt(ionization_variables.get_temperature()) /
                    ionization_variables.get_ionic_fraction(ION_H_n));
      x = random_generator.get_uniform_random_double();
      if (x < pHots) {
        // absorbed on the spot
        x = random_generator.get_uniform_random_double();
        if (x <= ionization_variables.get_reemission_probability(
                     REEMISSIONPROBABILITY_HYDROGEN)) {
          // H Ly c, like above
          new_frequency = _HLyc_spectrum.get_random_frequency(
              random_generator, ionization_variables.get_temperature());
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

/**
 * @brief Scatter the given Photon through an interaction with dust in the ISM.
 *
 * This routine will set a new direction for the Photon, and will use and update
 * its Stokes parameters.
 *
 * @param photon Photon to scatter.
 * @param random_generator RandomGenerator used to generate pseudo random
 * numbers.
 */
void PhotonSource::scatter(Photon &photon,
                           RandomGenerator &random_generator) const {

  // Stokes parameters
  double fi_old, fq_old, fu_old, fv_old;
  photon.get_stokes_parameters(fi_old, fq_old, fu_old, fv_old);

  const double fi = 1.;
  const double fq = fq_old / fi_old;
  const double fu = fu_old / fi_old;
  const double fv = fv_old / fi_old;

  const double term =
      _scattering_omg2 /
      (_scattering_omhgg +
       _scattering_thgg * random_generator.get_uniform_random_double());
  const double bmu = std::min(
      1., std::max(-1., _scattering_od2hgg * (_scattering_opg2 - term * term)));

  double fi_new, fq_new, fu_new, fv_new;
  CoordinateVector<> direction_new;
  double sin_theta_new, cos_theta_new, phi_new, sin_phi_new, cos_phi_new;

  if (std::abs(bmu) == 1.) {

    fi_new = fi;
    fq_new = fq;
    if (bmu == -1.) {
      fu_new = -fu;
    } else {
      fu_new = fu;
    }
    fv_new = fv;

    direction_new = photon.get_direction();
    photon.get_direction_parameters(sin_theta_new, cos_theta_new, phi_new,
                                    sin_phi_new, cos_phi_new);

  } else {

    const double weight = fi_old;

    double sin_theta_old, cos_theta_old, phi_old, sin_phi_old, cos_phi_old;
    photon.get_direction_parameters(sin_theta_old, cos_theta_old, phi_old,
                                    sin_phi_old, cos_phi_old);

    const double cosb2 = bmu * bmu;
    const double b = cosb2 - 1.;

    // dustmat
    const double p1 = _scattering_omg2 *
                      std::pow(_scattering_opg2 - _scattering_thgg * bmu, -1.5);
    const double p2 = -_scattering_pl * p1 * (1. - cosb2) / (1. + cosb2);
    const double p3 = 2. * p1 * bmu / (1. + cosb2);
    const double phi_temp = std::acos(bmu) * 180. / M_PI;
    const double f = 3.13 * phi_temp * std::exp(-7. * phi_temp / 180.);
    const double f2 = (phi_temp + _scattering_sc * f) * M_PI / 180.;
    const double c = std::cos(f2);
    const double c2 = c * c;
    const double p4 = -_scattering_pc * p1 * (1. - c2) / (1. + c2);

    const double a = p1;
    const double sinbt = std::sqrt(-b);
    const double ri1 = 2. * M_PI * random_generator.get_uniform_random_double();

    double a11, a12, a13, a21, a22, a23, a24, a31, a32, a33, a34, a42, a43, a44;

    if (ri1 > M_PI) {

      const double ri3 = 2. * M_PI - ri1;
      const double cosi3 = std::cos(ri3);
      const double sini3 = std::sin(ri3);
      const double sin2i3 = 2. * sini3 * cosi3;
      const double cos2i3 = 2. * cosi3 * cosi3 - 1.;

      cos_theta_new = cos_theta_old * bmu + sin_theta_old * sinbt * cosi3;
      double sini2, cosi2;
      if (std::abs(cos_theta_new) < 1.) {
        sin_theta_new = std::abs(std::sqrt(1. - cos_theta_new * cos_theta_new));
        sini2 = sini3 * sin_theta_old / sin_theta_new;
        const double bott = sin_theta_new * sinbt;
        cosi2 = cos_theta_old / bott - cos_theta_new * bmu / bott;
      } else {
        sin_theta_new = 0.;
        sini2 = 0.;
        if (cos_theta_new >= 1.) {
          cosi2 = -1.;
        } else {
          cosi2 = 1.;
        }
      }

      const double cosdph =
          std::min(1., std::max(-1., -cosi2 * cosi3 + sini2 * sini3 * bmu));
      phi_new = phi_old + std::acos(cosdph);
      if (phi_new > 2. * M_PI) {
        phi_new -= 2. * M_PI;
      }
      if (phi_new < 0.) {
        phi_new += 2. * M_PI;
      }

      const double sin2i2 = 2. * sini2 * cosi2;
      const double cos2i2 = 2. * cosi2 * cosi2 - 1.;
      const double sin2 = sin2i2 * sin2i3;
      const double cos2 = cos2i2 * cos2i3;
      const double sin2cos1 = sin2i2 * cos2i3;
      const double cos2sin1 = cos2i2 * sin2i3;

      a11 = p1;
      a12 = p2 * cos2i3;
      a13 = p2 * sin2i3;

      a21 = p2 * cos2i2;
      a22 = p1 * cos2 - p3 * sin2;
      a23 = p1 * cos2sin1 + p3 * sin2cos1;
      a24 = -p4 * sin2i2;

      a31 = -p2 * sin2i2;
      a32 = -p1 * sin2cos1 - p3 * cos2sin1;
      a33 = -p1 * sin2 + p3 * cos2;
      a34 = -p4 * cos2i2;

      a42 = -p4 * sin2i3;
      a43 = p4 * cos2i3;
      a44 = p3;

    } else {

      const double cosi1 = std::cos(ri1);
      const double sini1 = std::sin(ri1);
      const double sin2i1 = 2. * sini1 * cosi1;
      const double cos2i1 = 2. * cosi1 * cosi1 - 1.;

      cos_theta_new = cos_theta_old * bmu + sin_theta_old * sinbt * cosi1;
      double sini2, cosi2;
      if (std::abs(cos_theta_new) < 1.) {
        sin_theta_new = std::abs(std::sqrt(1. - cos_theta_new * cos_theta_new));
        sini2 = sini1 * sin_theta_old / sin_theta_new;
        const double bott = sin_theta_new * sinbt;
        cosi2 = cos_theta_old / bott - cos_theta_new * bmu / bott;
      } else {
        sin_theta_new = 0.;
        sini2 = 0.;
        if (cos_theta_new >= 1.) {
          cosi2 = -1.;
        } else {
          cosi2 = 1.;
        }
      }

      const double cosdph =
          std::min(1., std::max(-1., -cosi2 * cosi1 + sini2 * sini1 * bmu));
      phi_new = phi_old - std::acos(cosdph);
      if (phi_new > 2. * M_PI) {
        phi_new -= 2. * M_PI;
      }
      if (phi_new < 0.) {
        phi_new += 2. * M_PI;
      }

      const double sin2i2 = 2. * sini2 * cosi2;
      const double cos2i2 = 2. * cosi2 * cosi2 - 1.;
      const double sin2 = sin2i2 * sin2i1;
      const double cos2 = cos2i2 * cos2i1;
      const double sin2cos1 = sin2i2 * cos2i1;
      const double cos2sin1 = cos2i2 * sin2i1;

      a11 = p1;
      a12 = p2 * cos2i1;
      a13 = -p2 * sin2i1;

      a21 = p2 * cos2i2;
      a22 = p1 * cos2 - p3 * sin2;
      a23 = -p1 * cos2sin1 - p3 * sin2cos1;
      a24 = p4 * sin2i2;

      a31 = p2 * sin2i2;
      a32 = p1 * sin2cos1 + p3 * cos2sin1;
      a33 = -p1 * sin2 + p3 * cos2;
      a34 = -p4 * cos2i2;

      a42 = p4 * sin2i1;
      a43 = p4 * cos2i1;
      a44 = p3;
    }

    const double si = (a11 * fi + a12 * fq + a13 * fu) / a;
    const double sq = (a21 * fi + a22 * fq + a23 * fu + a24 * fv) / a;
    const double su = (a31 * fi + a32 * fq + a33 * fu + a34 * fv) / a;
    const double sv = (a42 * fq + a43 * fu + a44 * fv) / a;

    fi_new = si * weight;
    fq_new = sq * weight;
    fu_new = su * weight;
    fv_new = sv * weight;

    cos_phi_new = std::cos(phi_new);
    sin_phi_new = std::sin(phi_new);
    direction_new[0] = sin_theta_new * cos_phi_new;
    direction_new[1] = sin_theta_new * sin_phi_new;
    direction_new[2] = cos_theta_new;
  }

  photon.set_stokes_parameters(fi_new, fq_new, fu_new, fv_new);
  photon.set_direction(direction_new);
  photon.set_direction_parameters(sin_theta_new, cos_theta_new, phi_new,
                                  sin_phi_new, cos_phi_new);
}

/**
 * @brief Update the Stokes parameters for the given Photon for a scattering
 * interaction that scatters it towards an observer in the given direction.
 *
 * This routine will set a new direction for the Photon, and will use and update
 * its Stokes parameters.
 *
 * @param photon Photon to scatter.
 * @param direction Direction in which the photon is observed.
 * @param sint @f$\sin(\theta{})@f$ of the observation direction.
 * @param cost @f$\cos(\theta{})@f$ of the observation direction.
 * @param phi @f$\phi{}@f$ of the observation direction (in radians).
 * @param sinp @f$\sin(\phi{})@f$ of the observation direction.
 * @param cosp @f$\cos(\phi{})@f$ of the observation direction.
 * @return Weight of the scattered photon.
 */
double PhotonSource::scatter_towards(Photon &photon,
                                     const CoordinateVector<> direction,
                                     double sint, double cost, double phi,
                                     double sinp, double cosp) const {

  double cos_theta_old, sin_theta_old, phi_old, sin_phi_old, cos_phi_old;
  const CoordinateVector<> direction_old = photon.get_direction();
  photon.get_direction_parameters(sin_theta_old, cos_theta_old, phi_old,
                                  sin_phi_old, cos_phi_old);

  const double calpha =
      CoordinateVector<>::dot_product(direction, direction_old);
  const double cos_theta = cost;
  double sin_theta = sint;

  // Stokes parameters
  double fi_old, fq_old, fu_old, fv_old;
  photon.get_stokes_parameters(fi_old, fq_old, fu_old, fv_old);

  const double fi = 1.;
  const double fq = fq_old / fi_old;
  const double fu = fu_old / fi_old;
  const double fv = fv_old / fi_old;

  const double bmu = calpha;

  double fi_new, fq_new, fu_new, fv_new;

  if (std::abs(bmu) == 1.) {

    fi_new = fi;
    fq_new = fq;
    if (bmu == -1.) {
      fu_new = -fu;
    } else {
      fu_new = fu;
    }
    fv_new = fv;

  } else {

    const double weight = fi_old;

    const double cosb2 = bmu * bmu;
    const double b = cosb2 - 1.;

    // dustmat
    const double p1 = _scattering_omg2 *
                      std::pow(_scattering_opg2 - _scattering_thgg * bmu, -1.5);
    const double p2 = -_scattering_pl * p1 * (1. - cosb2) / (1. + cosb2);
    const double p3 = 2. * p1 * bmu / (1. + cosb2);
    const double phi_temp = std::acos(bmu) * 180. / M_PI;
    const double f = 3.13 * phi_temp * std::exp(-7. * phi_temp / 180.);
    const double f2 = (phi_temp + _scattering_sc * f) * M_PI / 180.;
    const double c = std::cos(f2);
    const double c2 = c * c;
    const double p4 = -_scattering_pc * p1 * (1. - c2) / (1. + c2);

    const double a = p1;
    const double sinbt = std::sqrt(-b);

    double ri1;
    if (sin_theta_old == 0.) {
      ri1 = M_PI;
    } else {
      const double cosi1 =
          (cos_theta - cos_theta_old * bmu) / (sin_theta_old * sinbt);
      const double sini1 = std::sin(phi_old - phi - M_PI) * sin_theta / sinbt;
      ri1 = std::atan2(sini1, cosi1) + M_PI;
    }

    double a11, a12, a13, a21, a22, a23, a24, a31, a32, a33, a34, a42, a43, a44;

    if (ri1 > M_PI) {

      const double ri3 = 2. * M_PI - ri1;
      const double cosi3 = std::cos(ri3);
      const double sini3 = std::sin(ri3);
      const double sin2i3 = 2. * sini3 * cosi3;
      const double cos2i3 = 2. * cosi3 * cosi3 - 1.;

      double sini2, cosi2;
      if (std::abs(cos_theta) < 1.) {
        sini2 = sini3 * sin_theta_old / sin_theta;
        const double bott = sin_theta * sinbt;
        cosi2 = cos_theta_old / bott - cos_theta * bmu / bott;
      } else {
        sini2 = 0.;
        if (cos_theta >= 1.) {
          cosi2 = -1.;
        } else {
          cosi2 = 1.;
        }
      }

      const double sin2i2 = 2. * sini2 * cosi2;
      const double cos2i2 = 2. * cosi2 * cosi2 - 1.;
      const double sin2 = sin2i2 * sin2i3;
      const double cos2 = cos2i2 * cos2i3;
      const double sin2cos1 = sin2i2 * cos2i3;
      const double cos2sin1 = cos2i2 * sin2i3;

      a11 = p1;
      a12 = p2 * cos2i3;
      a13 = p2 * sin2i3;

      a21 = p2 * cos2i2;
      a22 = p1 * cos2 - p3 * sin2;
      a23 = p1 * cos2sin1 + p3 * sin2cos1;
      a24 = -p4 * sin2i2;

      a31 = -p2 * sin2i2;
      a32 = -p1 * sin2cos1 - p3 * cos2sin1;
      a33 = -p1 * sin2 + p3 * cos2;
      a34 = -p4 * cos2i2;

      a42 = -p4 * sin2i3;
      a43 = p4 * cos2i3;
      a44 = p3;

    } else {

      const double cosi1 = std::cos(ri1);
      const double sini1 = std::sin(ri1);
      const double sin2i1 = 2. * sini1 * cosi1;
      const double cos2i1 = 2. * cosi1 * cosi1 - 1.;

      double sini2, cosi2;
      if (std::abs(cos_theta) < 1.) {
        sini2 = sini1 * sin_theta_old / sin_theta;
        const double bott = sin_theta * sinbt;
        cosi2 = cos_theta_old / bott - cos_theta * bmu / bott;
      } else {
        sini2 = 0.;
        if (cos_theta >= 1.) {
          cosi2 = -1.;
        } else {
          cosi2 = 1.;
        }
      }

      const double sin2i2 = 2. * sini2 * cosi2;
      const double cos2i2 = 2. * cosi2 * cosi2 - 1.;
      const double sin2 = sin2i2 * sin2i1;
      const double cos2 = cos2i2 * cos2i1;
      const double sin2cos1 = sin2i2 * cos2i1;
      const double cos2sin1 = cos2i2 * sin2i1;

      a11 = p1;
      a12 = p2 * cos2i1;
      a13 = -p2 * sin2i1;

      a21 = p2 * cos2i2;
      a22 = p1 * cos2 - p3 * sin2;
      a23 = -p1 * cos2sin1 - p3 * sin2cos1;
      a24 = p4 * sin2i2;

      a31 = p2 * sin2i2;
      a32 = p1 * sin2cos1 + p3 * cos2sin1;
      a33 = -p1 * sin2 + p3 * cos2;
      a34 = -p4 * cos2i2;

      a42 = p4 * sin2i1;
      a43 = p4 * cos2i1;
      a44 = p3;
    }

    const double si = (a11 * fi + a12 * fq + a13 * fu) / a;
    const double sq = (a21 * fi + a22 * fq + a23 * fu + a24 * fv) / a;
    const double su = (a31 * fi + a32 * fq + a33 * fu + a34 * fv) / a;
    const double sv = (a42 * fq + a43 * fu + a44 * fv) / a;

    fi_new = si * weight;
    fq_new = sq * weight;
    fu_new = su * weight;
    fv_new = sv * weight;
  }

  photon.set_direction(direction);
  photon.set_direction_parameters(sint, cost, phi, sinp, cosp);
  photon.set_stokes_parameters(fi_new, fq_new, fu_new, fv_new);

  return 0.25 * _scattering_omg2 *
         std::pow(_scattering_opg2 - _scattering_thgg * calpha, -1.5) / M_PI;
}
