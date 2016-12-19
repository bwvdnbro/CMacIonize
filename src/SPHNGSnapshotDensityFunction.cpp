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
 * @file SPHNGSnapshotDensityFunction.cpp
 *
 * @brief SPHNGSnapshotDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "SPHNGSnapshotDensityFunction.hpp"
#include "DensityValues.hpp"
#include "Octree.hpp"
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
double SPHNGSnapshotDensityFunction::kernel(const double q, const double h) {
  if (q < 1.) {
    double q2 = q * q;
    double h2 = h * h;
    double h3 = h2 * h;
    return (1. - 1.5 * q2 + 0.75 * q2 * q) / M_PI / h3;
  } else if (q < 2.) {
    double c = 2. - q;
    double c2 = c * c;
    double h2 = h * h;
    double h3 = h * h2;
    return 0.25 * c2 * c / M_PI / h3;
  } else {
    return 0.;
  }
}

/**
 * @brief Constructor.
 *
 * @param filename Name of the file to read.
 */
SPHNGSnapshotDensityFunction::SPHNGSnapshotDensityFunction(
    std::string filename) {
  std::ifstream file(filename, std::ios::binary | std::ios::in);

  if (!file) {
    cmac_error("Unable to open file \"%s\"!", filename.c_str());
  }

  // read header

  // extract length of first line
  int length;
  // sizeof(int)+1, since get reads up to n-1 characters...
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  // make sure we have arrived at the next length bit
  int length2;
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  // second block
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  // third block
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  int number;
  file.read(reinterpret_cast< char * >(&number), sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  // fourth block: tags
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  // fifth block
  // Will's Fortran code:
  //  read(10) nparttot,n1,n2,nreassign,naccrete,nkill,nblocks,iyr,&
  //        idum,(iv(i),i=1,NTAB),iplanetesimals,irotpot,idragscheme
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  int data[44];
  file.read(reinterpret_cast< char * >(data), 44 * sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }
  unsigned int totnumpart = data[0];

  // skip 3 blocks
  for (unsigned int i = 0; i < 3; ++i) {
    file.read(reinterpret_cast< char * >(&length), sizeof(int));
    // skip length bytes
    file.seekg(length, std::ios_base::cur);
    file.read(reinterpret_cast< char * >(&length2), sizeof(int));
    if (length != length2) {
      cmac_error("Something went wrong!");
    }
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&number), sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  unsigned long iuniquemax;
  if (number == 1) {
    // skip block
    file.read(reinterpret_cast< char * >(&length), sizeof(int));
    // skip length bytes
    file.seekg(length, std::ios_base::cur);
    file.read(reinterpret_cast< char * >(&length2), sizeof(int));
    if (length != length2) {
      cmac_error("Something went wrong!");
    }

    // read iuniquemax
    file.read(reinterpret_cast< char * >(&length), sizeof(int));
    file.read(reinterpret_cast< char * >(&iuniquemax), sizeof(unsigned long));
    file.read(reinterpret_cast< char * >(&length2), sizeof(int));
    if (length != length2) {
      cmac_error("Something went wrong!");
    }
  } else {
    iuniquemax = totnumpart;
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&number), sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  std::vector< std::string > tags;
  for (int i = 0; i < number; ++i) {
    char c[17];
    file.read(c, 16);
    c[16] = '\0';
    unsigned int j = 16;
    while (j > 0 && c[j - 1] == ' ') {
      --j;
      c[j] = '\0';
    }
    if (j == 0) {
      // flag tag empty
      c[0] = '0';
    }
    std::string tag(c);
    tags.push_back(tag);
  }
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  std::map< std::string, double > headerdict;
  for (unsigned int i = 0; i < tags.size(); ++i) {
    double rval;
    file.read(reinterpret_cast< char * >(&rval), sizeof(double));
    headerdict[tags[i]] = rval;
  }
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&number), sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&number), sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  // skip block
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  std::vector< double > units;
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  for (int i = 0; i < number; ++i) {
    double unit;
    file.read(reinterpret_cast< char * >(&unit), sizeof(double));
    units.push_back(unit);
  }
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&number), sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  // done reading header!

  // read block
  unsigned long npart;
  unsigned int nums[8];
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&npart), sizeof(unsigned long));
  file.read(reinterpret_cast< char * >(nums), 8 * sizeof(unsigned int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  unsigned long nptmass;
  unsigned int numsink[8];
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&nptmass), sizeof(unsigned long));
  file.read(reinterpret_cast< char * >(numsink), 8 * sizeof(unsigned int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  // skip block
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  /// continue on line '  read(10) isteps(icount+1:icount+npart)' of Will's
  /// rbin.f90 file...
  std::vector< int > isteps(npart);
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&isteps[0]), npart * sizeof(int));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  if (nums[0] >= 2) {
    // skip 2 blocks
    for (unsigned int i = 0; i < 2; ++i) {
      file.read(reinterpret_cast< char * >(&length), sizeof(int));
      // skip length bytes
      file.seekg(length, std::ios_base::cur);
      file.read(reinterpret_cast< char * >(&length2), sizeof(int));
      if (length != length2) {
        cmac_error("Something went wrong!");
      }
    }
  }

  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  char *c = new char[length + 1];
  file.read(c, length);
  c[length] = '\0';
  unsigned int j = length - 1;
  while (c[j] == ' ') {
    c[j] = '\0';
    --j;
  }
  std::string tag(c);
  delete[] c;
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  if (tag != "iphase") {
    cmac_error("Wrong tag: \"%s\" (expected \"iphase\")!", tag.c_str());
  }

  std::vector< char > iphase(npart);
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(&iphase[0], npart);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  if (nums[4] >= 1) {
    // skip iunique block
    file.read(reinterpret_cast< char * >(&length), sizeof(int));
    // skip length bytes
    file.seekg(length, std::ios_base::cur);
    file.read(reinterpret_cast< char * >(&length2), sizeof(int));
    if (length != length2) {
      cmac_error("Something went wrong!");
    }
    file.read(reinterpret_cast< char * >(&length), sizeof(int));
    // skip length bytes
    file.seekg(length, std::ios_base::cur);
    file.read(reinterpret_cast< char * >(&length2), sizeof(int));
    if (length != length2) {
      cmac_error("Something went wrong!");
    }
  }

  std::vector< double > x(npart);
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&x[0]), npart * sizeof(double));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  std::vector< double > y(npart);
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&y[0]), npart * sizeof(double));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  std::vector< double > z(npart);
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&z[0]), npart * sizeof(double));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  std::vector< double > m(npart);
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&m[0]), npart * sizeof(double));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  std::vector< double > h(npart);
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  // skip length bytes
  file.seekg(length, std::ios_base::cur);
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }
  file.read(reinterpret_cast< char * >(&length), sizeof(int));
  file.read(reinterpret_cast< char * >(&h[0]), npart * sizeof(double));
  file.read(reinterpret_cast< char * >(&length2), sizeof(int));
  if (length != length2) {
    cmac_error("Something went wrong!");
  }

  double unit_length = UnitConverter::to_SI< QUANTITY_LENGTH >(units[0], "cm");
  double unit_mass = UnitConverter::to_SI< QUANTITY_MASS >(units[1], "g");

  unsigned int ngaspart = 0;
  for (unsigned int i = 0; i < iphase.size(); ++i) {
    if (iphase[i] == 0) {
      ++ngaspart;
    }
  }

  Box partbox;
  partbox.get_anchor()[0] = DBL_MAX;
  partbox.get_anchor()[1] = DBL_MAX;
  partbox.get_anchor()[2] = DBL_MAX;
  partbox.get_sides()[0] = -DBL_MAX;
  partbox.get_sides()[1] = -DBL_MAX;
  partbox.get_sides()[2] = -DBL_MAX;

  _positions.resize(ngaspart);
  _masses.resize(ngaspart);
  _smoothing_lengths.resize(ngaspart);
  unsigned int index = 0;
  for (unsigned int i = 0; i < npart; ++i) {
    if (iphase[i] == 0) {
      _positions[index][0] = x[i] * unit_length;
      _positions[index][1] = y[i] * unit_length;
      _positions[index][2] = z[i] * unit_length;
      partbox.get_anchor() =
          CoordinateVector<>::min(partbox.get_anchor(), _positions[index]);
      partbox.get_sides() =
          CoordinateVector<>::max(partbox.get_sides(), _positions[index]);
      _masses[index] = m[i] * unit_mass;
      _smoothing_lengths[index] = h[i] * unit_length;
      ++index;
    }
  }

  partbox.get_sides() -= partbox.get_anchor();
  // add some margin to the box
  partbox.get_anchor() -= 0.01 * partbox.get_sides();
  partbox.get_sides() *= 1.02;

  _octree = new Octree(_positions, partbox, false);
  _octree->set_auxiliaries(_smoothing_lengths, Octree::max< double >);
}

/**
 * @brief Destructor.
 *
 * Clean up the octree.
 */
SPHNGSnapshotDensityFunction::~SPHNGSnapshotDensityFunction() {
  delete _octree;
}

/**
 * @brief Get the position of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return CoordinateVector<> containing the position of that particle (in m).
 */
CoordinateVector<>
SPHNGSnapshotDensityFunction::get_position(unsigned int index) {

  return _positions[index];
}

/**
 * @brief Get the mass of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return Mass of the particle (in kg).
 */
double SPHNGSnapshotDensityFunction::get_mass(unsigned int index) {
  return _masses[index];
}

/**
 * @brief Get the smoothing length of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return Smoothing length of the particle (in m).
 */
double SPHNGSnapshotDensityFunction::get_smoothing_length(unsigned int index) {
  return _smoothing_lengths[index];
}

/**
 * @brief Get the DensityValues at the given position.
 *
 * @param position Coordinates of a position (in m).
 * @return DensityValues at that position (in SI units).
 */
DensityValues SPHNGSnapshotDensityFunction::
operator()(CoordinateVector<> position) const {
  DensityValues cell;

  double density = 0.;
  std::vector< unsigned int > ngbs = _octree->get_ngbs(position);
  const unsigned int numngbs = ngbs.size();
  for (unsigned int i = 0; i < numngbs; ++i) {
    unsigned int index = ngbs[i];
    double r;
    r = (position - _positions[index]).norm();
    double h = _smoothing_lengths[index];
    double q = r / h;
    double m = _masses[index];
    double splineval = m * kernel(q, h);
    density += splineval;
  }

  // convert density to particle density (assuming hydrogen only)
  cell.set_total_density(density / 1.6737236e-27);
  // TODO: other quantities

  return cell;
}
