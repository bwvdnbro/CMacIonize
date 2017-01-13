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
 * @file testIsotropicContinuousPhotonSource.cpp
 *
 * @brief Unit test for the IsotropicContinuousPhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "IsotropicContinuousPhotonSource.hpp"
#include "RandomGenerator.hpp"
#include <fstream>

/*! @brief Number of angular bins. */
#define NUMANGLEBIN 16

/*! @brief Number of directional bins. */
#define NUMDIRBIN 100

/**
 * @brief Get the intersection point of the line through the given point and
 * with the given direction with the sphere surrounding the entire box.
 *
 * @param f Point on the line.
 * @param d Direction of the line.
 * @return Intersection point.
 */
CoordinateVector<> get_sphere_position(CoordinateVector<> f,
                                       CoordinateVector<> d) {
  const double R = 2.;
  double R2 = R * R;
  double ddotf = d.x() * f.x() + d.y() * f.y() + d.z() * f.z();
  double f2 = f.norm2();

  double disc = 4. * ddotf * ddotf - 4. * (f2 - R2);
  double l1 = -ddotf + 0.5 * std::sqrt(disc);
  double l2 = -ddotf - 0.5 * std::sqrt(disc);
  double l;

  assert_condition(l1 * l2 < 0.);

  if (l1 < 0.) {
    l = l1;
  } else {
    l = l2;
  }

  CoordinateVector<> i = f + l * d;
  assert_values_equal(i.norm(), R);
  return i;
}

/**
 * @brief Unit test for the IsotropicContinuousPhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box box(CoordinateVector<>(-0.5), CoordinateVector<>(1.));
  RandomGenerator random_generator(44);
  IsotropicContinuousPhotonSource source(box);

  // to see the angular dependence, this number should be a factor 100 larger
  // but then the unit test takes too long to complete
  const unsigned int numphoton = 1000000;

  unsigned int bins[6][101];
  for (unsigned int i = 0; i < 6; ++i) {
    for (unsigned int j = 0; j < 101; ++j) {
      bins[i][j] = 0;
    }
  }
  double sbins[NUMANGLEBIN][NUMDIRBIN + 1];
  for (unsigned int i = 0; i < NUMANGLEBIN; ++i) {
    for (unsigned int j = 0; j < NUMDIRBIN + 1; ++j) {
      sbins[i][j] = 0.;
    }
  }

  std::ofstream sfile("isotropic_sphere.txt");

  CoordinateVector<> average_position;
  CoordinateVector<> average_direction;
  for (unsigned int i = 0; i < numphoton; ++i) {
    std::pair< CoordinateVector<>, CoordinateVector<> > posdir =
        source.get_random_incoming_direction(random_generator);
    average_position += posdir.first;
    average_direction += posdir.second;

    CoordinateVector<> spos = get_sphere_position(posdir.first, posdir.second);

    // we make angular bins for all rays in a thin band around the equator
    // plane
    if (std::abs(spos.z()) < 0.2) {
      // convert the position to a position on the unit sphere
      spos *= 0.5;
      assert_values_equal(spos.norm(), 1.);
      // get the phi angle (assuming x = cos(phi) and y = sin(phi))
      double phi = std::atan2(spos.y(), spos.x());
      // phi lies in the range [-pi, pi]
      // find the index of phi in the NUMANGLEBIN bins, starting with bin
      // [ -pi-0.5*(2.*pi)/NUMANGLEBIN, -pi+0.5*(2.*pi)/NUMANGLEBIN [
      int ibin = 0.5 * NUMANGLEBIN * (phi + M_PI + M_PI / NUMANGLEBIN) / M_PI;
      assert_condition(ibin >= 0);
      assert_condition(ibin <= NUMANGLEBIN);
      // rebin to correct for offset
      if (ibin >= NUMANGLEBIN) {
        ibin -= NUMANGLEBIN;
      }

      // we now find the outgoing angle, relative to the bin angle
      double phi_ref = -M_PI + (2. * M_PI * ibin) / NUMANGLEBIN;

      double phi_out = std::atan2(posdir.second.y(), posdir.second.x());
      // + pi, since we expect the outgoing angle to be opposite to the bin
      // angle
      phi_out -= phi_ref + M_PI;
      if (phi_out < -M_PI) {
        phi_out += 2. * M_PI;
      }
      if (phi_out >= M_PI) {
        phi_out -= 2. * M_PI;
      }
      // we only bin between -1. and 1.
      int jbin = 0.5 * NUMDIRBIN * (phi_out + 1.);
      assert_condition(jbin >= 0);
      assert_condition(jbin < NUMDIRBIN);
      sbins[ibin][jbin] += 1.;
    }

    for (unsigned int ip = 0; ip < 3; ++ip) {
      unsigned int ibin = 100 * (posdir.first[ip] + 0.5);
      ++bins[ip][ibin];
    }
    for (unsigned int id = 0; id < 3; ++id) {
      unsigned int ibin = 50 * (posdir.second[id] + 1.);
      ++bins[id + 3][ibin];
    }
  }
  average_position /= numphoton;
  average_direction /= numphoton;

  std::ofstream ofile("isotropic_bins.txt");
  unsigned int binsum[6] = {0};
  for (unsigned int i = 0; i < 101; ++i) {
    ofile << i;
    for (unsigned int j = 0; j < 6; ++j) {
      ofile << "\t" << bins[j][i];
      binsum[j] += bins[j][i];
    }
    ofile << "\n";
  }

  for (unsigned int i = 0; i < NUMANGLEBIN; ++i) {
    assert_condition(sbins[i][NUMDIRBIN] == 0);
    for (unsigned int j = 0; j < NUMDIRBIN; ++j) {
      sbins[i][NUMDIRBIN] += sbins[i][j];
    }
    sbins[i][NUMDIRBIN] = 1. / sbins[i][NUMDIRBIN];
    for (unsigned int j = 0; j < NUMDIRBIN; ++j) {
      sbins[i][j] *= sbins[i][NUMDIRBIN];
    }
  }

  for (unsigned int i = 0; i < NUMDIRBIN + 1; ++i) {
    sfile << i;
    for (unsigned int j = 0; j < NUMANGLEBIN; ++j) {
      sfile << "\t" << sbins[j][i];
    }
    sfile << "\n";
  }

  return 0;
}
