/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DustScattering.cpp
 *
 * @brief DustScattering implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DustScattering.hpp"
#include "Photon.hpp"
#include "RandomGenerator.hpp"
#include <cmath>

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
void DustScattering::scatter(Photon &photon,
                             RandomGenerator &random_generator) const {

  // Stokes parameters
  double I_old, Q_old, U_old, V_old;
  photon.get_stokes_parameters(I_old, Q_old, U_old, V_old);

  // first randomly sample a scattering angle
  // we randomly sample from the Henyey-Greenstein phase function using
  // Witt (1977), equation (19) to get the cosine of the scattering angle, where
  // the scattering angle is measured in the plane formed by the old and new
  // photon path
  const double term =
      _scattering_omg2 /
      (_scattering_omhgg +
       _scattering_thgg * random_generator.get_uniform_random_double());
  const double cos_scattering_angle = std::min(
      1., std::max(-1., _scattering_od2hgg * (_scattering_opg2 - term * term)));

  double I_new, Q_new, U_new, V_new;
  CoordinateVector<> direction_new;
  double sin_theta_new, cos_theta_new, phi_new, sin_phi_new, cos_phi_new;

  if (std::abs(cos_scattering_angle) == 1.) {

    // if we have forward/backward scattering, the polarization change is
    // trivial

    direction_new = photon.get_direction();
    photon.get_direction_parameters(sin_theta_new, cos_theta_new, phi_new,
                                    sin_phi_new, cos_phi_new);

    I_new = I_old;
    Q_new = Q_old;
    if (cos_scattering_angle == -1.) {
      // backward scattering
      U_new = -U_old;
      // theta' = pi - theta
      // phi' = pi + phi
      // cos(theta') = -cos(theta)
      // sin(theta') = sin(theta)
      // cos(phi') = -cos(phi)
      // sin(phi') = -sin(phi)
      direction_new[0] = -direction_new[0];
      direction_new[1] = -direction_new[1];
      direction_new[2] = -direction_new[2];
      cos_theta_new = -cos_theta_new;
      sin_phi_new = -sin_phi_new;
      cos_phi_new = -cos_phi_new;
      phi_new += M_PI;
    } else {
      // forward scattering
      U_new = U_old;
    }
    V_new = V_old;

  } else {

    // general scattering angle
    // in this case, we need to draw a random phi angle, which is the angle that
    // determines the orientation of the scattering plane

    // normalize the Stokes parameters
    // note that this is not necessary, as we multiply the weight back in in the
    // end...
    const double weight = I_old;
    const double weight_inv = 1. / weight;

    const double I = 1.;
    const double Q = Q_old * weight_inv;
    const double U = U_old * weight_inv;
    const double V = V_old * weight_inv;

    double sin_theta_old, cos_theta_old, phi_old, sin_phi_old, cos_phi_old;
    photon.get_direction_parameters(sin_theta_old, cos_theta_old, phi_old,
                                    sin_phi_old, cos_phi_old);

    const double cos2_scattering_angle =
        cos_scattering_angle * cos_scattering_angle;

    // White (1979), equation (3)
    const double P1 =
        _scattering_omg2 *
        std::pow(_scattering_opg2 - _scattering_thgg * cos_scattering_angle,
                 -1.5);
    const double one_over_one_plus_cos2_scattering_angle =
        1. / (1. + cos2_scattering_angle);
    // White (1979), equation (4)
    const double P2 = -_scattering_pl * P1 * (1. - cos2_scattering_angle) *
                      one_over_one_plus_cos2_scattering_angle;
    // White (1979), equation (5)
    const double P3 = 2. * P1 * cos_scattering_angle *
                      one_over_one_plus_cos2_scattering_angle;
    // White (1979), equation (6)
    const double scattering_angle = std::acos(cos_scattering_angle);
    const double cos_skewed_scattering_angle =
        std::cos(scattering_angle +
                 _scattering_sc * 3.13 * scattering_angle *
                     std::exp(-7. * scattering_angle * M_1_PI));
    const double cos2_skewed_scattering_angle =
        cos_skewed_scattering_angle * cos_skewed_scattering_angle;
    const double P4 = -_scattering_pc * P1 *
                      (1. - cos2_skewed_scattering_angle) /
                      (1. + cos2_skewed_scattering_angle);

    // interestingly, the common factor P1 is divided out in the end...
    // this means that we could just set P1 to 1 and drop the common factor
    // everywhere...
    const double a_inv = 1. / P1;
    const double sin_scattering_angle =
        std::sqrt(std::max(0., 1. - cos2_scattering_angle));
    const double random_phi1 =
        2. * M_PI * random_generator.get_uniform_random_double();

    // rotation matrix elements
    double A00, A01, A02, A10, A11, A12, A13, A20, A21, A22, A23, A31, A32, A33;

    if (random_phi1 > M_PI) {

      const double random_phi3 = 2. * M_PI - random_phi1;
      const double cos_random_phi3 = std::cos(random_phi3);
      const double sin_random_phi3 = std::sin(random_phi3);
      const double sin_2random_phi3 = 2. * sin_random_phi3 * cos_random_phi3;
      const double cos_2random_phi3 =
          2. * cos_random_phi3 * cos_random_phi3 - 1.;

      // Yusef-Zadeh, Morris & White (1984), equation (16)
      cos_theta_new = cos_theta_old * cos_scattering_angle +
                      sin_theta_old * sin_scattering_angle * cos_random_phi3;
      double sini2, cosi2;
      if (std::abs(cos_theta_new) < 1.) {
        sin_theta_new = std::abs(std::sqrt(1. - cos_theta_new * cos_theta_new));
        sini2 = sin_random_phi3 * sin_theta_old / sin_theta_new;
        cosi2 = (cos_theta_old - cos_theta_new * cos_scattering_angle) /
                (sin_theta_new * sin_scattering_angle);
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
          std::min(1., std::max(-1., -cosi2 * cos_random_phi3 +
                                         sini2 * sin_random_phi3 *
                                             cos_scattering_angle));
      phi_new = phi_old + std::acos(cosdph);
      if (phi_new > 2. * M_PI) {
        phi_new -= 2. * M_PI;
      }
      if (phi_new < 0.) {
        phi_new += 2. * M_PI;
      }

      const double sin2i2 = 2. * sini2 * cosi2;
      const double cos2i2 = 2. * cosi2 * cosi2 - 1.;
      const double sin2 = sin2i2 * sin_2random_phi3;
      const double cos2 = cos2i2 * cos_2random_phi3;
      const double sin2cos1 = sin2i2 * cos_2random_phi3;
      const double cos2sin1 = cos2i2 * sin_2random_phi3;

      // double matrix multiplication: Code & Whitney (1995), equation (2)
      A00 = P1;
      A01 = P2 * cos_2random_phi3;
      A02 = P2 * sin_2random_phi3;

      A10 = P2 * cos2i2;
      A11 = P1 * cos2 - P3 * sin2;
      A12 = P1 * cos2sin1 + P3 * sin2cos1;
      A13 = -P4 * sin2i2;

      A20 = -P2 * sin2i2;
      A21 = -P1 * sin2cos1 - P3 * cos2sin1;
      A22 = -P1 * sin2 + P3 * cos2;
      A23 = -P4 * cos2i2;

      A31 = -P4 * sin_2random_phi3;
      A32 = P4 * cos_2random_phi3;
      A33 = P3;

    } else {

      const double cosi1 = std::cos(random_phi1);
      const double sini1 = std::sin(random_phi1);
      const double sin2i1 = 2. * sini1 * cosi1;
      const double cos2i1 = 2. * cosi1 * cosi1 - 1.;

      cos_theta_new = cos_theta_old * cos_scattering_angle +
                      sin_theta_old * sin_scattering_angle * cosi1;
      double sini2, cosi2;
      if (std::abs(cos_theta_new) < 1.) {
        sin_theta_new = std::abs(std::sqrt(1. - cos_theta_new * cos_theta_new));
        sini2 = sini1 * sin_theta_old / sin_theta_new;
        cosi2 = (cos_theta_old - cos_theta_new * cos_scattering_angle) /
                (sin_theta_new * sin_scattering_angle);
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
          std::min(1., std::max(-1., -cosi2 * cosi1 +
                                         sini2 * sini1 * cos_scattering_angle));
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

      // double matrix multiplication: Code & Whitney (1995), equation (2)
      A00 = P1;
      A01 = P2 * cos2i1;
      A02 = -P2 * sin2i1;

      A10 = P2 * cos2i2;
      A11 = P1 * cos2 - P3 * sin2;
      A12 = -P1 * cos2sin1 - P3 * sin2cos1;
      A13 = P4 * sin2i2;

      A20 = P2 * sin2i2;
      A21 = P1 * sin2cos1 + P3 * cos2sin1;
      A22 = -P1 * sin2 + P3 * cos2;
      A23 = -P4 * cos2i2;

      A31 = P4 * sin2i1;
      A32 = P4 * cos2i1;
      A33 = P3;
    }

    const double si = (A00 * I + A01 * Q + A02 * U) * a_inv;
    const double sq = (A10 * I + A11 * Q + A12 * U + A13 * V) * a_inv;
    const double su = (A20 * I + A21 * Q + A22 * U + A23 * V) * a_inv;
    const double sv = (A31 * Q + A32 * U + A33 * V) * a_inv;

    I_new = si * weight;
    Q_new = sq * weight;
    U_new = su * weight;
    V_new = sv * weight;

    cos_phi_new = std::cos(phi_new);
    sin_phi_new = std::sin(phi_new);
    direction_new[0] = sin_theta_new * cos_phi_new;
    direction_new[1] = sin_theta_new * sin_phi_new;
    direction_new[2] = cos_theta_new;
  }

  photon.set_stokes_parameters(I_new, Q_new, U_new, V_new);
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
double DustScattering::scatter_towards(Photon &photon,
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

  const double bmu = calpha;

  double fi_new, fq_new, fu_new, fv_new;

  if (std::abs(bmu) == 1.) {

    fi_new = fi_old;
    fq_new = fq_old;
    if (bmu == -1.) {
      fu_new = -fu_old;
    } else {
      fu_new = fu_old;
    }
    fv_new = fv_old;

  } else {

    const double weight = fi_old;
    const double weight_inv = 1. / weight;

    const double fi = 1.;
    const double fq = fq_old * weight_inv;
    const double fu = fu_old * weight_inv;
    const double fv = fv_old * weight_inv;

    const double cosb2 = bmu * bmu;
    const double b = cosb2 - 1.;

    // White (1979), equation (3)
    const double p1 = _scattering_omg2 *
                      std::pow(_scattering_opg2 - _scattering_thgg * bmu, -1.5);
    const double opcosb2inv = 1. / (1. + cosb2);
    // White (1979), equation (4)
    const double p2 = -_scattering_pl * p1 * (1. - cosb2) * opcosb2inv;
    // White (1979), equation (5)
    const double p3 = 2. * p1 * bmu * opcosb2inv;
    // White (1979), equation (6)
    const double phi_temp = std::acos(bmu) * 180. * M_1_PI;
    const double f = 3.13 * phi_temp * std::exp(-7. * phi_temp / 180.);
    const double f2 = (phi_temp + _scattering_sc * f) * M_PI / 180.;
    const double c = std::cos(f2);
    const double c2 = c * c;
    const double p4 = -_scattering_pc * p1 * (1. - c2) / (1. + c2);

    const double a_inv = 1. / p1;
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

      // double matrix multiplication: Code & Whitney (1995), equation (2)
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

      // double matrix multiplication: Code & Whitney (1995), equation (2)
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

    const double si = (a11 * fi + a12 * fq + a13 * fu) * a_inv;
    const double sq = (a21 * fi + a22 * fq + a23 * fu + a24 * fv) * a_inv;
    const double su = (a31 * fi + a32 * fq + a33 * fu + a34 * fv) * a_inv;
    const double sv = (a42 * fq + a43 * fu + a44 * fv) * a_inv;

    fi_new = si * weight;
    fq_new = sq * weight;
    fu_new = su * weight;
    fv_new = sv * weight;
  }

  photon.set_direction(direction);
  photon.set_direction_parameters(sint, cost, phi, sinp, cosp);
  photon.set_stokes_parameters(fi_new, fq_new, fu_new, fv_new);

  // the weight of the scattering event is given by the relevant Henyey-
  // Greenstein phase function, White (1979), equation (3)
  return 0.25 * _scattering_omg2 *
         std::pow(_scattering_opg2 - _scattering_thgg * calpha, -1.5) * M_1_PI;
}
