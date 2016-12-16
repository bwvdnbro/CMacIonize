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
#include <fstream>
#include <map>

/**
 * @brief Constructor.
 *
 * @param filename Name of the file to read.
 */
SPHNGSnapshotDensityFunction::SPHNGSnapshotDensityFunction(
    std::string filename) {
  std::ifstream file(filename, std::ios::binary | std::ios::in);

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
  cmac_status("%lu", iuniquemax);

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
    unsigned int j = 15;
    while (c[j] == ' ') {
      c[j] = '\0';
      --j;
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
}

/**
 * @brief Get the DensityValues at the given position.
 *
 * @param position Coordinates of a position (in m).
 * @return DensityValues at that position (in SI units).
 */
DensityValues SPHNGSnapshotDensityFunction::
operator()(CoordinateVector<> position) const {
  return DensityValues();
}
