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
 * @file PhantomSnapshotDensityFunction.cpp
 *
 * @brief PhantomSnapshotDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "PhantomSnapshotDensityFunction.hpp"
#include "DensityValues.hpp"
#include "Log.hpp"
#include "Octree.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"

#include <cfloat>
#include <fstream>
#include <map>

/**
 * @brief Cubic spline kernel.
 *
 * As in Price, 2007, Publications of the Astronomical Society of Australia, 24,
 * 159 (equation 5).
 *
 * @param q Distance between the kernel center and the evaluation point, in
 * units of the smoothing length.
 * @param h Smoothing length.
 * @return Value of the kernel.
 */
double PhantomSnapshotDensityFunction::kernel(const double q, const double h) {
  if (q < 1.) {
    double q2 = q * q;
    double h2 = h * h;
    double h3 = h2 * h;
    return (1. - 1.5 * q2 + 0.75 * q2 * q) / (M_PI * h3);
  } else if (q < 2.) {
    double c = 2. - q;
    double c2 = c * c;
    double h2 = h * h;
    double h3 = h * h2;
    return 0.25 * c2 * c / (M_PI * h3);
  } else {
    return 0.;
  }
}

/**
 * @brief Function that gives the 3D integral of the kernel of
 * a particle for a given vertex of a cell face.
 *
 * @param phi Azimuthal angle of the vertex.
 * @param r0 Distance from the particle to the face of the cell.
 * @param R_0 Distance from the orthogonal projection of the particle
 * onto the face of the cell to a side of the face (containing the vertex).
 * @param h The kernel smoothing length of the particle.
 * @return The integral of the kernel for the given vertex.
 */

double PhantomSnapshotDensityFunction::full_integral(double phi, double r0,
                                                     double R_0, double h) {

  double B1, B2, B3, mu, a, logs, u;
  double full_int;
  double a2, cosp, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3, tanp;
  double r2, R, linedist2, phi1, phi2, h2;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double D2, D3;

  B1 = 0.0;
  B2 = 0.0;
  B3 = 0.0;

  if (r0 == 0.0)
    return 0.0;
  if (R_0 == 0.0)
    return 0.0;
  if (phi == 0.0)
    return 0.0;

  h2 = h * h;
  r03 = r0 * r0 * r0;
  r0h2 = r0 / h * r0 / h;
  r0h3 = r0h2 * r0 / h;
  r0h_2 = h / r0 * h / r0;
  r0h_3 = r0h_2 * h / r0;

  // Setting up the B1, B2, B3 constants of integration.

  if (r0 >= 2.0 * h) {
    B3 = h2 * h / 4.;
  } else if (r0 > h) {
    B3 = r03 / 4. *
         (-4. / 3. + (r0 / h) - 0.3 * r0h2 + 1. / 30. * r0h3 -
          1. / 15. * r0h_3 + 8. / 5. * r0h_2);
    B2 =
        r03 / 4. *
        (-4. / 3. + (r0 / h) - 0.3 * r0h2 + 1. / 30. * r0h3 - 1. / 15. * r0h_3);
  } else {
    B3 = r03 / 4. * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3 + 7. / 5. * r0h_2);
    B2 = r03 / 4. * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3 - 1. / 5. * r0h_2);
    B1 = r03 / 4. * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3);
  }

  a = R_0 / r0;
  a2 = a * a;

  linedist2 = r0 * r0 + R_0 * R_0;
  R = R_0 / cos(phi);
  r2 = r0 * r0 + R * R;

  full_int = 0.0;
  D2 = 0.0;
  D3 = 0.0;

  if (linedist2 <= h2) {
    ////// phi1 business /////
    phi1 = acos(R_0 / sqrt(h * h - r0 * r0));

    cosp = cos(phi1);
    cosp2 = cosp * cosp;
    mu = cosp / a / sqrt(1. + cosp2 / a2);

    tanp = tan(phi1);

    I0 = phi1;
    I_2 = phi1 + a2 * tanp;
    I_4 = phi1 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sin(phi1) * sqrt(1 - mu * mu);
    logs = log(1 + u) - log(1 - u);
    I1 = atan(u / a);

    I_1 = a / 2. * logs + I1;
    I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
    I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. *
                    ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) +
                     3. * logs);

    D2 = -1. / 6. * I_2 + 0.25 * (r0 / h) * I_3 - 0.15 * r0h2 * I_4 +
         1. / 30. * r0h3 * I_5 - 1. / 60. * r0h_3 * I1 + (B1 - B2) / r03 * I0;

    ////// phi2 business /////
    phi2 = acos(R_0 / sqrt(4.0 * h * h - r0 * r0));

    cosp = cos(phi2);
    cosp2 = cosp * cosp;
    mu = cosp / a / sqrt(1. + cosp2 / a2);

    tanp = tan(phi2);

    I0 = phi2;
    I_2 = phi2 + a2 * tanp;
    I_4 = phi2 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sin(phi2) * sqrt(1 - mu * mu);
    logs = log(1 + u) - log(1 - u);
    I1 = atan(u / a);

    I_1 = a / 2. * logs + I1;
    I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
    I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. *
                    ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) +
                     3. * logs);

    D3 = 1. / 3. * I_2 - 0.25 * (r0 / h) * I_3 + 3. / 40. * r0h2 * I_4 -
         1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 +
         D2;
  } else if (linedist2 <= 4.0 * h2) {
    ////// phi2 business /////
    phi2 = acos(R_0 / sqrt(4.0 * h * h - r0 * r0));

    cosp = cos(phi2);
    cosp2 = cosp * cosp;
    mu = cosp / a / sqrt(1. + cosp2 / a2);

    tanp = tan(phi2);

    I0 = phi2;
    I_2 = phi2 + a2 * tanp;
    I_4 = phi2 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sin(phi2) * sqrt(1 - mu * mu);
    logs = log(1 + u) - log(1 - u);
    I1 = atan(u / a);

    I_1 = a / 2. * logs + I1;
    I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
    I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. *
                    ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) +
                     3. * logs);

    D3 = 1. / 3. * I_2 - 0.25 * (r0 / h) * I_3 + 3. / 40. * r0h2 * I_4 -
         1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 +
         D2;
  }

  //////////////////////////////////
  // Calculating I_n expressions. //
  //////////////////////////////////

  cosp = cos(phi);
  cosp2 = cosp * cosp;
  mu = cosp / a / sqrt(1. + cosp2 / a2);

  tanp = tan(phi);

  I0 = phi;
  I_2 = phi + a2 * tanp;
  I_4 = phi + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

  u = sin(phi) * sqrt(1 - mu * mu);
  logs = log(1 + u) - log(1 - u);
  I1 = atan(u / a);

  I_1 = a / 2. * logs + I1;
  I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
  I_5 = I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

  // Calculating the integral expression.

  if (r2 < h2) {
    full_int = r0h3 / M_PI *
               (1. / 6. * I_2 - 3. / 40. * r0h2 * I_4 + 1. / 40. * r0h3 * I_5 +
                B1 / r03 * I0);
  } else if (r2 < 4.0 * h2) {
    full_int = r0h3 / M_PI *
               (0.25 * (4. / 3. * I_2 - (r0 / h) * I_3 + 0.3 * r0h2 * I_4 -
                        1. / 30. * r0h3 * I_5 + 1. / 15. * r0h_3 * I1) +
                B2 / r03 * I0 + D2);
  } else {
    full_int = r0h3 / M_PI * (-0.25 * r0h_3 * I1 + B3 / r03 * I0 + D3);
  }

  return full_int;
}

/**
 * @brief Function that calculates the mass contribution of a particle
 * towards the total mass of a cell.
 *
 * @param cell Geometrical information about the cell.
 * @param particle The particle position.
 * @param h The kernel smoothing length of the particle.
 * @return The mass contribution of the particle to the cell.
 */

double PhantomSnapshotDensityFunction::mass_contribution(
    const Cell &cell, const CoordinateVector<> particle, const double h) {

  double M, Msum;

  Msum = 0.;
  M = 0.;

  std::vector< Face > face_vector = cell.get_faces();

  // Loop over each face of a cell.
  for (size_t i = 0; i < face_vector.size(); i++) {

    CoordinateVector<> vert_position1;
    CoordinateVector<> projected_particle;
    double r0 = 0.;
    double ar0 = 0.;
    double s2 = 0.;

    // Loop over the vertices of each face.
    for (Face::Vertices j = face_vector[i].first_vertex();
         j != face_vector[i].last_vertex(); ++j) {
      if (j == face_vector[i].first_vertex()) { // Calculating the distance from
                                                // particle to each face of a
                                                // cell
        // http://mathinsight.org/distance_point_plane
        // http://mathinsight.org/forming_planes
        Face::Vertices j_twin = j;
        vert_position1 = j_twin.get_position();
        CoordinateVector<> vert_position2 = (++j_twin).get_position();
        CoordinateVector<> vert_position3 = (++j_twin).get_position();

        const double A = (vert_position2[1] - vert_position1[1]) *
                             (vert_position3[2] - vert_position1[2]) -
                         (vert_position3[1] - vert_position1[1]) *
                             (vert_position2[2] - vert_position1[2]);
        const double B = (vert_position2[2] - vert_position1[2]) *
                             (vert_position3[0] - vert_position1[0]) -
                         (vert_position3[2] - vert_position1[2]) *
                             (vert_position2[0] - vert_position1[0]);
        const double C = (vert_position2[0] - vert_position1[0]) *
                             (vert_position3[1] - vert_position1[1]) -
                         (vert_position3[0] - vert_position1[0]) *
                             (vert_position2[1] - vert_position1[1]);
        const double D = -A * vert_position1[0] - B * vert_position1[1] -
                         C * vert_position1[2];

        const double norm = sqrt(A * A + B * B + C * C);
        r0 = (A * particle[0] + B * particle[1] + C * particle[2] + D) / norm;
        ar0 = fabs(r0);

        // Calculate of the orthogonal projection of the particle position onto
        // the face.
        projected_particle[0] = particle[0] - r0 * A / norm;
        projected_particle[1] = particle[1] - r0 * B / norm;
        projected_particle[2] = particle[2] - r0 * C / norm;

        // s2 contains information about the orientation of the face vertices.
        s2 = vert_position1[0] * (vert_position2[1] * vert_position3[2] -
                                  vert_position2[2] * vert_position3[1]) +
             vert_position1[1] * (vert_position2[2] * vert_position3[0] -
                                  vert_position2[0] * vert_position3[2]) +
             vert_position1[2] * (vert_position2[0] * vert_position3[1] -
                                  vert_position2[1] * vert_position3[0]);
      }

      Face::Vertices j_twin = j;
      const CoordinateVector<> vert_position2 = j_twin.get_position();
      CoordinateVector<> vert_position3;

      if (++j_twin == face_vector[i].last_vertex()) {
        vert_position3 = vert_position1;
      } else {
        vert_position3 = j_twin.get_position();
      }

      const double r23 = (vert_position2 - vert_position3).norm();
      const double r12 = (projected_particle - vert_position2).norm();
      const double r13 = (projected_particle - vert_position3).norm();
      const double cosa = ((vert_position3[0] - vert_position2[0]) *
                               (projected_particle[0] - vert_position2[0]) +
                           (vert_position3[1] - vert_position2[1]) *
                               (projected_particle[1] - vert_position2[1]) +
                           (vert_position3[2] - vert_position2[2]) *
                               (projected_particle[2] - vert_position2[2])) /
                          r12 / r23;

      double R_0 = 0.;
      double phi1 = 0.;
      double phi2 = 0.;

      if (fabs(cosa) < 1.0) {
        R_0 = r12 * sqrt(1 - cosa * cosa);
      } else {
        if (fabs(cosa) - 1.0 < 0.00001) {
          R_0 = 0.0;
        } else {
          printf("Error: cosa > 1: %g\n", cosa);
        }
      }

      const double s1 =
          projected_particle[0] * (vert_position2[1] * vert_position3[2] -
                                   vert_position2[2] * vert_position3[1]) +
          projected_particle[1] * (vert_position2[2] * vert_position3[0] -
                                   vert_position2[0] * vert_position3[2]) +
          projected_particle[2] * (vert_position2[0] * vert_position3[1] -
                                   vert_position2[1] * vert_position3[0]);

      if (R_0 < r12) {
        phi1 = acos(R_0 / r12);
      } else {
        if ((R_0 - r12) / h < 0.00001) {
          phi1 = 0.0;
        } else {
          printf("Error: R0 > r12: %g\n", R_0 - r12);
        }
      }
      if (R_0 < r13) {
        phi2 = acos(R_0 / r13);
      } else {
        if ((R_0 - r13) / h < 0.00001) {
          phi2 = 0.0;
        } else {
          printf("Error: R0 > r13: %g\n", R_0 - r13);
        }
      }

      // Find out if the vertex integral will contribute positively or
      // negatively to the cell mass.
      if (s1 * s2 * r0 <= 0) {
        M = -1.;
      } else {
        M = 1.;
      }

      // Calculate the vertex integral.
      if ((r12 * sin(phi1) >= r23) || (r13 * sin(phi2) >= r23)) {
        if (phi1 >= phi2) {
          M = M * (full_integral(phi1, ar0, R_0, h) -
                   full_integral(phi2, ar0, R_0, h));
        } else {
          M = M * (full_integral(phi2, ar0, R_0, h) -
                   full_integral(phi1, ar0, R_0, h));
        }
      } else {
        M = M * (full_integral(phi1, ar0, R_0, h) +
                 full_integral(phi2, ar0, R_0, h));
      }
      Msum = Msum + M;
    }
  }

  return Msum;
}

/**
 * @brief Constructor.
 *
 * Reads the file and stores the particle variables in internal arrays. Does not
 * construct the Octree, this is done in initialize().
 *
 * @param filename Name of the file to read.
 * @param initial_temperature Initial temperature of the gas (in K).
 * @param use_new_algorithm Use the new mapping algorithm?
 * @param use_periodic_box Use a periodic box for the density mapping?
 * @param binary_dump Dump the particle positions and densities to a binary
 * file?
 * @param binary_dump_name Name of the file in which the particle positions and
 * densities are dumped.
 * @param log Log to write logging info to.
 */
PhantomSnapshotDensityFunction::PhantomSnapshotDensityFunction(
    std::string filename, double initial_temperature,
    const bool use_new_algorithm, const bool use_periodic_box,
    const bool binary_dump, const std::string binary_dump_name, Log *log)
    : _octree(nullptr), _initial_temperature(initial_temperature),
      _use_new_algorithm(use_new_algorithm),
      _use_periodic_box(use_periodic_box), _log(log) {

  if (binary_dump && binary_dump_name == "") {
    cmac_error("No valid name provided for binary dump file!");
  }

  std::ifstream file(filename, std::ios::binary | std::ios::in);

  if (!file) {
    cmac_error("Unable to open file \"%s\"!", filename.c_str());
  }

  // skip the first block: it contains garbage
  skip_block(file);
  // read the second block: it contains the file identity
  // we only support file identities starting with 'F'
  // there is currently no support for untagged files ('T')
  std::string fileident;
  read_block(file, fileident);
  if (fileident[0] != 'F') {
    cmac_error("Unsupported Phantom snapshot format: %s!", fileident.c_str());
  }
  bool tagged = true;
  if (fileident[1] != 'T') {
    tagged = false;
  }
  // for now, we assume a tagged file
  cmac_assert(tagged);

  // read header blocks
  std::map< std::string, int32_t > ints = read_dict< int32_t >(file, tagged);
  std::map< std::string, int8_t > int8s = read_dict< int8_t >(file, tagged);
  std::map< std::string, int16_t > int16s = read_dict< int16_t >(file, tagged);
  std::map< std::string, int32_t > int32s = read_dict< int32_t >(file, tagged);
  std::map< std::string, int64_t > int64s = read_dict< int64_t >(file, tagged);
  std::map< std::string, double > reals = read_dict< double >(file, tagged);
  std::map< std::string, float > real4s = read_dict< float >(file, tagged);
  std::map< std::string, double > real8s = read_dict< double >(file, tagged);

  if (_log) {
    _log->write_info("Phantom header:");
    _log->write_info("ints:");
    for (auto it = ints.begin(); it != ints.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
    _log->write_info("int8s:");
    for (auto it = int8s.begin(); it != int8s.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
    _log->write_info("int16s:");
    for (auto it = int16s.begin(); it != int16s.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
    _log->write_info("int32s:");
    for (auto it = int32s.begin(); it != int32s.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
    _log->write_info("int64s:");
    for (auto it = int64s.begin(); it != int64s.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
    _log->write_info("reals:");
    for (auto it = reals.begin(); it != reals.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
    _log->write_info("real4s:");
    for (auto it = real4s.begin(); it != real4s.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
    _log->write_info("real8s:");
    for (auto it = real8s.begin(); it != real8s.end(); ++it) {
      _log->write_info(it->first, ": ", it->second);
    }
  }

  int_fast32_t numblocks = ints["nblocks"];
  // for now, we assume all data is stored in a single block
  cmac_assert(numblocks == 1);
  int32_t narraylengths;
  read_block(file, narraylengths);
  narraylengths /= numblocks;
  cmac_assert(narraylengths > 1);
  cmac_assert(narraylengths < 4);

  // create temporary vectors to store particle data
  const uint_fast32_t numpart = int64s["npartoftype"];
  std::vector< double > x(numpart), y(numpart), z(numpart);
  std::vector< float > h(numpart);

  // read the number of variables in the single block for each array
  std::vector< int64_t > varnumber(narraylengths);
  std::vector< std::vector< int32_t > > varnums(narraylengths,
                                                std::vector< int32_t >(8));
  for (int_fast32_t iarray = 0; iarray < narraylengths; ++iarray) {

    read_block(file, varnumber[iarray], varnums[iarray]);
  }

  // now read the data
  // we only read the coordinates x, y, z and the smoothing length h
  for (int_fast32_t iarray = 0; iarray < narraylengths; ++iarray) {
    for (uint_fast8_t idata = 0; idata < 8; ++idata) {
      for (int_fast32_t i = 0; i < varnums[iarray][idata]; ++i) {
        std::string tag;
        read_block(file, tag);
        if (tag == "x") {
          read_block(file, x);
        } else if (tag == "y") {
          read_block(file, y);
        } else if (tag == "z") {
          read_block(file, z);
        } else if (tag == "h") {
          read_block(file, h);
        } else {
          skip_block(file);
        }
      }
    }
  }

  // the particle mass is constant in Phantom
  // the mass unit is given in CGS, we convert to SI
  const double pmass = reals["massoftype"] * real8s["umass"] * 0.001;
  // get the length unit in CGS and convert to SI
  const double unit_length_in_SI = real8s["udist"] * 0.01;

  // initialise variables to store box dimensions in
  _partbox.get_anchor()[0] = DBL_MAX;
  _partbox.get_anchor()[1] = DBL_MAX;
  _partbox.get_anchor()[2] = DBL_MAX;
  _partbox.get_sides()[0] = -DBL_MAX;
  _partbox.get_sides()[1] = -DBL_MAX;
  _partbox.get_sides()[2] = -DBL_MAX;

  // now add the particle data to the right internal vectors
  // we also convert units
  _positions.resize(numpart);
  _smoothing_lengths.resize(numpart);
  _masses.resize(numpart);
  for (uint_fast32_t ipart = 0; ipart < numpart; ++ipart) {
    const CoordinateVector<> p(x[ipart] * unit_length_in_SI,
                               y[ipart] * unit_length_in_SI,
                               z[ipart] * unit_length_in_SI);
    _positions[ipart] = p;
    _smoothing_lengths[ipart] = h[ipart] * unit_length_in_SI;
    _masses[ipart] = pmass;

    // keep track of the min and max position for all particles
    _partbox.get_anchor() = CoordinateVector<>::min(_partbox.get_anchor(), p);
    _partbox.get_sides() = CoordinateVector<>::max(_partbox.get_sides(), p);
  }

  if (binary_dump) {
    std::ofstream dump_file(binary_dump_name);
    for (uint_fast32_t i = 0; i < _positions.size(); ++i) {
      const double xi = _positions[i].x();
      const double yi = _positions[i].y();
      const double zi = _positions[i].z();
      const double hi = _smoothing_lengths[i];
      const double mi = _masses[i];
      // we approximate the particle density using the mass and smoothing length
      const double rhoi = 0.75 * M_1_PI * mi / (hi * hi * hi);
      dump_file.write(reinterpret_cast< const char * >(&xi), sizeof(double));
      dump_file.write(reinterpret_cast< const char * >(&yi), sizeof(double));
      dump_file.write(reinterpret_cast< const char * >(&zi), sizeof(double));
      dump_file.write(reinterpret_cast< const char * >(&rhoi), sizeof(double));
    }
  }

  // convert max positions to box sides and add safety margins
  _partbox.get_sides() -= _partbox.get_anchor();
  // add some margin to the box
  _partbox.get_anchor() -= 0.01 * _partbox.get_sides();
  _partbox.get_sides() *= 1.02;

  if (_use_periodic_box) {
    _partbox.get_anchor()[0] = reals["xmin"] * unit_length_in_SI;
    _partbox.get_anchor()[1] = reals["ymin"] * unit_length_in_SI;
    _partbox.get_anchor()[2] = reals["zmin"] * unit_length_in_SI;
    _partbox.get_sides()[0] =
        (reals["xmax"] - reals["xmin"]) * unit_length_in_SI;
    _partbox.get_sides()[1] =
        (reals["ymax"] - reals["ymin"]) * unit_length_in_SI;
    _partbox.get_sides()[2] =
        (reals["zmax"] - reals["zmin"]) * unit_length_in_SI;
  }

  if (_log) {
    _log->write_status("Snapshot contains ", _positions.size(),
                       " gas particles.");
    _log->write_status("Total particle mass is ", pmass * _masses.size(),
                       "kg.");
    const std::string periodic = (_use_periodic_box) ? " periodic" : "";
    _log->write_status(
        "Will create octree in", periodic, " box with anchor [",
        _partbox.get_anchor().x(), " m, ", _partbox.get_anchor().y(), " m, ",
        _partbox.get_anchor().z(), " m] and sides [", _partbox.get_sides().x(),
        " m, ", _partbox.get_sides().y(), " m, ", _partbox.get_sides().z(),
        " m]...");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name fo the snapshot file (required)
 *  - initial temperature: Initial temperature of the gas (default: 8000. K)
 *  - use new algorithm: Use the new algorithm of Maya Petkova? (default: false)
 *  - use periodic box: Use a periodic box for the density mapping? (default:
 *    false)
 *  - binary dump: Produce a binary dump of the particle positions and
 *    densities (default: false)?
 *  - binary dump name: Name of the file in which to dump the particle positions
 *    and densities.
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
PhantomSnapshotDensityFunction::PhantomSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : PhantomSnapshotDensityFunction(
          params.get_filename("DensityFunction:filename"),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "DensityFunction:initial temperature", "8000. K"),
          params.get_value< bool >("DensityFunction:use new algorithm", false),
          params.get_value< bool >("DensityFunction:use periodic box", false),
          params.get_value< bool >("DensityFunction:binary dump", false),
          params.get_value< std::string >("DensityFunction:binary dump name",
                                          ""),
          log) {}

/**
 * @brief Destructor.
 *
 * Clean up the octree.
 */
PhantomSnapshotDensityFunction::~PhantomSnapshotDensityFunction() {
  delete _octree;
}

/**
 * @brief This routine constructs the internal Octree that is used for neighbour
 * finding.
 */
void PhantomSnapshotDensityFunction::initialize() {

  _octree = new Octree(_positions, _partbox, _use_periodic_box);
  std::vector< double > h2s = _smoothing_lengths;
  for (uint_fast32_t i = 0; i < _smoothing_lengths.size(); ++i) {
    h2s[i] *= 2.;
  }
  _octree->set_auxiliaries(h2s, Octree::max< double >);
}

/**
 * @brief Get the number of particles in the snapshot.
 *
 * @return Number of particles in the snapshot.
 */
uint_fast32_t PhantomSnapshotDensityFunction::get_number_of_particles() const {
  return _positions.size();
}

/**
 * @brief Get the position of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return CoordinateVector<> containing the position of that particle (in m).
 */
CoordinateVector<>
PhantomSnapshotDensityFunction::get_position(const uint_fast32_t index) const {
  return _positions[index];
}

/**
 * @brief Get the mass of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return Mass of the particle (in kg).
 */
double
PhantomSnapshotDensityFunction::get_mass(const uint_fast32_t index) const {
  return _masses[index];
}

/**
 * @brief Get the smoothing length of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return Smoothing length of the particle (in m).
 */
double PhantomSnapshotDensityFunction::get_smoothing_length(
    const uint_fast32_t index) const {
  return _smoothing_lengths[index];
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues PhantomSnapshotDensityFunction::
operator()(const Cell &cell) const {

  DensityValues values;
  if (_use_new_algorithm) {

    CoordinateVector<> position = cell.get_cell_midpoint();

    // Find the vertex that is furthest away from the cell midpoint.
    std::vector< Face > face_vector = cell.get_faces();
    double radius = 0.0;
    for (size_t i = 0; i < face_vector.size(); i++) {
      for (Face::Vertices j = face_vector[i].first_vertex();
           j != face_vector[i].last_vertex(); ++j) {
        double distance = j.get_position().norm();
        if (distance > radius)
          radius = distance;
      }
    }

    // Find the neighbours that are contained inside of a sphere of centre the
    // cell midpoint
    // and radius given by the distance to the furthest vertex.
    std::vector< uint_fast32_t > ngbs =
        _octree->get_ngbs_sphere(position, radius);
    const size_t numngbs = ngbs.size();

    double density = 0.;

    // Loop over all the neighbouring particles and calculate their mass
    // contributions.
    for (size_t i = 0; i < numngbs; i++) {
      const uint_fast32_t index = ngbs[i];
      const double h = _smoothing_lengths[index];
      const CoordinateVector<> particle = _positions[index];
      density += mass_contribution(cell, particle, h) * _masses[index];
    }

    // Divide the cell mass by the cell volume to get density.
    density = density / cell.get_volume();

    // convert density to particle density (assuming hydrogen only)
    values.set_number_density(density / 1.6737236e-27);
    // TODO: other quantities
    // temporary values
    values.set_temperature(_initial_temperature);
    values.set_ionic_fraction(ION_H_n, 1.e-6);
#ifdef HAS_HELIUM
    values.set_ionic_fraction(ION_He_n, 1.e-6);
#endif

  } else {

    const CoordinateVector<> position = cell.get_cell_midpoint();

    double density = 0.;
    std::vector< uint_fast32_t > ngbs = _octree->get_ngbs(position);
    const size_t numngbs = ngbs.size();
    for (size_t i = 0; i < numngbs; ++i) {
      const uint_fast32_t index = ngbs[i];
      double r;
      if (_use_periodic_box) {
        r = _partbox.periodic_distance(position, _positions[index]).norm();
      } else {
        r = (position - _positions[index]).norm();
      }
      const double h = _smoothing_lengths[index];
      const double q = r / h;
      const double m = _masses[index];
      const double splineval = m * kernel(q, h);
      density += splineval;
    }

    // convert density to particle density (assuming hydrogen only)
    values.set_number_density(density / 1.6737236e-27);
    // TODO: other quantities
    // temporary values
    values.set_temperature(_initial_temperature);
    values.set_ionic_fraction(ION_H_n, 1.e-6);
#ifdef HAS_HELIUM
    values.set_ionic_fraction(ION_He_n, 1.e-6);
#endif
  }

  return values;
}
