/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SPHNGVoronoiGeneratorDistribution.cpp
 *
 * @brief SPHNGVoronoiGeneratorDistribution implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "SPHNGVoronoiGeneratorDistribution.hpp"
#include "SPHNGSnapshotUtilities.hpp"

#include <fstream>

/**
 * @brief Constructor.
 *
 * @param filename Name of the SPHNG snapshot file to read.
 * @param log Log to write logging info to.
 */
SPHNGVoronoiGeneratorDistribution::SPHNGVoronoiGeneratorDistribution(
    const std::string filename, Log *log)
    : _current_number(0) {

  std::ifstream file(filename, std::ios::binary | std::ios::in);

  if (!file) {
    cmac_error("Unable to open file \"%s\"!", filename.c_str());
  }

  // read header

  // skip the first block: it contains garbage
  SPHNGSnapshotUtilities::skip_block(file);
  // read the second block: it contains the file identity
  // we only support file identities starting with 'F'
  // there is currently no support for untagged files ('T')
  std::string fileident;
  SPHNGSnapshotUtilities::read_block(file, fileident);
  if (fileident[0] != 'F') {
    cmac_error("Unsupported SPHNG snapshot format: %s!", fileident.c_str());
  }
  bool tagged = true;
  if (fileident[1] != 'T') {
    tagged = false;
  }

  if (log) {
    if (tagged) {
      log->write_info("File is tagged.");
    } else {
      log->write_info("File is not tagged.");
    }
  }

  // the third, fourth and fifth block contain a dictionary of particle numbers
  // the code below reads it
  std::map< std::string, uint32_t > numbers =
      SPHNGSnapshotUtilities::read_dict< uint32_t >(file, tagged);
  if (log) {
    log->write_info("Number header:");
    for (auto it = numbers.begin(); it != numbers.end(); ++it) {
      log->write_info(it->first, ": ", it->second);
    }
  }
  if (!tagged) {
    const size_t numnumbers = numbers.size();
    numbers["nparttot"] = numbers["tag"];
    if (numnumbers == 6) {
      numbers["nblocks"] = 1;
    } else {
      numbers["nblocks"] = numbers["tag6"];
    }
  }
  const uint_fast32_t numpart = numbers["nparttot"];
  const uint_fast32_t numblock = numbers["nblocks"];
  if (log) {
    log->write_info("nparttot: ", numbers["nparttot"]);
    log->write_info("nblocks: ", numbers["nblocks"]);
  }

  // skip 3 blocks
  // in the example file I got from Will, all these blocks contain a single
  // integer with value zero. They supposedly correspond to blocks that are
  // absent
  SPHNGSnapshotUtilities::skip_block(file);
  SPHNGSnapshotUtilities::skip_block(file);
  SPHNGSnapshotUtilities::skip_block(file);

  // the next three blocks are a dictionary containing the highest unique index
  // in the snapshot. If the first block contains 0, then the next 2 blocks are
  // absent.
  int32_t number;
  SPHNGSnapshotUtilities::read_block(file, number);

  if (number == 1) {
    // skip blocks
    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::skip_block(file);
  }

  // the next three blocks contain a number of double precision floating point
  // values that we don't use. They can be read in with the code below, but we
  // just skip them.
  if (log) {
    log->write_info("Header dictionary:");
    std::map< std::string, double > headerdict =
        SPHNGSnapshotUtilities::read_dict< double >(file, tagged);
    for (auto it = headerdict.begin(); it != headerdict.end(); ++it) {
      log->write_info(it->first, ": ", it->second);
    }
  } else {
    SPHNGSnapshotUtilities::skip_block(file);
    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::skip_block(file);
  }

  // the next block again corresponds to a block that is absent from Will's
  // example file
  SPHNGSnapshotUtilities::skip_block(file);

  // the next three blocks contain the units
  std::map< std::string, double > units =
      SPHNGSnapshotUtilities::read_dict< double >(file, tagged);
  if (log) {
    log->write_info("Unit block:");
    for (auto it = units.begin(); it != units.end(); ++it) {
      log->write_info(it->first, ": ", it->second);
    }
  }
  if (!tagged) {
    units["udist"] = units["tag"];
    units["umass"] = units["tag1"];
  }
  if (log) {
    log->write_info("udist: ", units["udist"]);
    log->write_info("umass: ", units["umass"]);
  }

  // the last header block is again absent from Will's file
  SPHNGSnapshotUtilities::skip_block(file);

  // done reading header!

  const double unit_length =
      UnitConverter::to_SI< QUANTITY_LENGTH >(units["udist"], "cm");

  _generator_positions.reserve(numpart);

  // read blocks
  for (uint_fast32_t iblock = 0; iblock < numblock; ++iblock) {

    if (log) {
      log->write_info("Block ", iblock);
    }

    uint64_t npart;
    std::vector< uint32_t > nums(8);
    SPHNGSnapshotUtilities::read_block(file, npart, nums);

    uint64_t nptmass;
    std::vector< uint32_t > numssink(8);
    SPHNGSnapshotUtilities::read_block(file, nptmass, numssink);

    if (log) {
      log->write_info("npart: ", npart);
      log->write_info("nums: ", nums[0], " ", nums[1], " ", nums[2], " ",
                      nums[3], " ", nums[4], " ", nums[5], " ", nums[6], " ",
                      nums[7]);
      log->write_info("nptmass: ", nptmass);
      log->write_info("numssink: ", numssink[0], " ", numssink[1], " ",
                      numssink[2], " ", numssink[3], " ", numssink[4], " ",
                      numssink[5], " ", numssink[6], " ", numssink[7]);
    }

    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }

    std::vector< int32_t > isteps(npart);
    SPHNGSnapshotUtilities::read_block(file, isteps);

    if (nums[0] >= 2) {
      // skip 2 blocks
      if (tagged) {
        SPHNGSnapshotUtilities::skip_block(file);
      }
      SPHNGSnapshotUtilities::skip_block(file);
    }

    std::string tag;
    if (tagged) {
      SPHNGSnapshotUtilities::read_block(file, tag);
    } else {
      tag = "iphase";
    }

    if (tag != "iphase") {
      cmac_error("Wrong tag: \"%s\" (expected \"iphase\")!", tag.c_str());
    }

    std::vector< int8_t > iphase(npart);
    SPHNGSnapshotUtilities::read_block(file, iphase);

    //    std::vector<uint64_t> iunique(npart);
    if (nums[4] >= 1) {
      // skip iunique block
      if (tagged) {
        SPHNGSnapshotUtilities::skip_block(file);
      }
      SPHNGSnapshotUtilities::skip_block(file);
      //      if(tagged){
      //        skip_block(file);
      //      }
      //      read_block(file, iunique);
    }

    std::vector< double > x(npart);
    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::read_block(file, x);

    std::vector< double > y(npart);
    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::read_block(file, y);

    std::vector< double > z(npart);
    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::read_block(file, z);

    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::skip_block(file);

    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::skip_block(file);

    // skip velocity and thermal energy blocks
    for (uint_fast8_t i = 0; i < 4; ++i) {
      if (tagged) {
        SPHNGSnapshotUtilities::skip_block(file);
      }
      SPHNGSnapshotUtilities::skip_block(file);
    }

    // skip density block
    if (tagged) {
      SPHNGSnapshotUtilities::skip_block(file);
    }
    SPHNGSnapshotUtilities::skip_block(file);

    // skip igrad related blocks
    for (uint_fast32_t i = 0; i < nums[6] - 1; ++i) {
      if (tagged) {
        SPHNGSnapshotUtilities::skip_block(file);
      }
      SPHNGSnapshotUtilities::skip_block(file);
    }

    // skip sink particle data
    for (uint_fast8_t i = 0; i < 10; ++i) {
      if (tagged) {
        SPHNGSnapshotUtilities::skip_block(file);
      }
      SPHNGSnapshotUtilities::skip_block(file);
    }

    for (uint_fast32_t i = 0; i < npart; ++i) {
      if (iphase[i] == 0) {
        const CoordinateVector<> position(
            x[i] * unit_length, y[i] * unit_length, z[i] * unit_length);
        _generator_positions.push_back(position);
      }
    }
  }

  if (log) {
    log->write_info("position size: ", _generator_positions.size());
    log->write_info("number of particles: ", numpart);
  }

  // done reading file
  // we would like to reduce memory usage at this point by shrinking the
  // internal vectors, using shrink_to_fit. However, the Intel compiler for some
  // reason complains about not finding this function (while still supporting
  // all other C++11 features), so we disable it here.
  // Note that enabling the code below might reduce the memory imprint,
  // especially if large snapshots are read in (if you use a suitable compiler,
  // that is).
  // BEGIN OF CODE
  // _positions.shrink_to_fit();
  // _masses.shrink_to_fit();
  // _smoothing_lengths.shrink_to_fit();
  // END OF CODE
}
