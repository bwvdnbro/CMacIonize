/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2017 Maya Petkova (map32@st-andrews.ac.uk)
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
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#include "SPHNGSnapshotDensityFunction.hpp"
#include "DensityValues.hpp"
#include "Log.hpp"
#include "Octree.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"
#include <cfloat>
#include <fstream>
#include <map>

// Maya
#include "VoronoiFace.hpp"
#include "VoronoiGrid.hpp"

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
 * Reads the file and stores the particle variables in internal arrays. Does not
 * construct the Octree, this is done in initialize().
 *
 * @param filename Name of the file to read.
 * @param initial_temperature Initial temperature of the gas (in K).
 * @param write_stats Flag indicating whether or not to write a file with
 * neighbour statistics.
 * @param stats_numbin Number of logarithmic bins used to keep track of
 * interneighbour distances.
 * @param stats_mindist Minimum interneighbour distance bin (in m; needs to be
 * non-zero, as we use logarithmic binning).
 * @param stats_maxdist Maximum interneighbour distance bin (in m).
 * @param stats_filename Name of the file with neighbour statistics that will be
 * written out.
 * @param use_new_algorithm Use the new mapping algorithm?
 * @param log Log to write logging info to.
 */
SPHNGSnapshotDensityFunction::SPHNGSnapshotDensityFunction(
    std::string filename, double initial_temperature, bool write_stats,
    unsigned int stats_numbin, double stats_mindist, double stats_maxdist,
    std::string stats_filename, bool use_new_algorithm, Log *log)
    : _use_new_algorithm(use_new_algorithm), _octree(nullptr),
      _initial_temperature(initial_temperature), _stats_numbin(stats_numbin),
      _stats_mindist(stats_mindist), _stats_maxdist(stats_maxdist),
      _stats_filename(stats_filename), _log(log) {
  std::ifstream file(filename, std::ios::binary | std::ios::in);

  if (!file) {
    cmac_error("Unable to open file \"%s\"!", filename.c_str());
  }

  // read header

  // skip the first block: it contains garbage
  skip_block(file);
  // read the second block: it contains the file identity
  // we only support file identities starting with 'F'
  // there is currently no support for untagged files ('T')
  std::string fileident;
  read_block(file, fileident);
  if (fileident[0] != 'F') {
    cmac_error("Unsupported SPHNG snapshot format: %s!", fileident.c_str());
  }
  bool tagged = true;
  if (fileident[1] != 'T') {
    tagged = false;
  }

  // the third, fourth and fifth block contain a dictionary of particle numbers
  // the code below reads it
  // since we currently don't use these values, we skip these 3 blocks instead
  std::map< std::string, unsigned int > numbers =
      read_dict< unsigned int >(file, tagged);
  if (!tagged) {
    unsigned int numnumbers = numbers.size();
    numbers["nparttot"] = numbers["tag"];
    if (numnumbers == 6) {
      numbers["nblocks"] = 1;
    } else {
      numbers["nblocks"] = numbers["tag6"];
    }
  }
  unsigned int numpart = numbers["nparttot"];
  unsigned int numblock = numbers["nblocks"];
  //  skip_block(file);
  //  if (tagged) {
  //    skip_block(file);
  //  }
  //  skip_block(file);

  // skip 3 blocks
  // in the example file I got from Will, all these blocks contain a single
  // integer with value zero. They supposedly correspond to blocks that are
  // absent
  skip_block(file);
  skip_block(file);
  skip_block(file);

  // the next three blocks are a dictionary containing the highest unique index
  // in the snapshot. If the first block contains 0, then the next 2 blocks are
  // absent.
  int number;
  read_block(file, number);

  if (number == 1) {
    // skip blocks
    if (tagged) {
      skip_block(file);
    }
    skip_block(file);
  }

  // the next three blocks contain a number of double precision floating point
  // values that we don't use. They can be read in with the code below, but we
  // just skip them.
  //  std::map< std::string, double > headerdict = read_dict< double >(file);
  skip_block(file);
  if (tagged) {
    skip_block(file);
  }
  skip_block(file);

  // the next block again corresponds to a block that is absent from Will's
  // example file
  skip_block(file);

  // the next three blocks contain the units
  std::map< std::string, double > units = read_dict< double >(file, tagged);
  if (!tagged) {
    units["udist"] = units["tag"];
    units["umass"] = units["tag1"];
  }

  // the last header block is again absent from Will's file
  skip_block(file);

  // done reading header!

  _partbox.get_anchor()[0] = DBL_MAX;
  _partbox.get_anchor()[1] = DBL_MAX;
  _partbox.get_anchor()[2] = DBL_MAX;
  _partbox.get_sides()[0] = -DBL_MAX;
  _partbox.get_sides()[1] = -DBL_MAX;
  _partbox.get_sides()[2] = -DBL_MAX;

  Box<> rawunitsbox;
  rawunitsbox.get_anchor()[0] = DBL_MAX;
  rawunitsbox.get_anchor()[1] = DBL_MAX;
  rawunitsbox.get_anchor()[2] = DBL_MAX;
  rawunitsbox.get_sides()[0] = -DBL_MAX;
  rawunitsbox.get_sides()[1] = -DBL_MAX;
  rawunitsbox.get_sides()[2] = -DBL_MAX;

  double unit_length =
      UnitConverter::to_SI< QUANTITY_LENGTH >(units["udist"], "cm");
  double unit_mass = UnitConverter::to_SI< QUANTITY_MASS >(units["umass"], "g");

  _positions.reserve(numpart);
  _masses.reserve(numpart);
  _smoothing_lengths.reserve(numpart);
  //  std::vector<uint64_t> all_iunique;
  //  all_iunique.reserve(numpart);

  // read blocks
  for (unsigned int iblock = 0; iblock < numblock; ++iblock) {
    uint64_t npart;
    std::vector< unsigned int > nums(8);
    read_block(file, npart, nums);

    uint64_t nptmass;
    std::vector< unsigned int > numssink(8);
    read_block(file, nptmass, numssink);

    if (tagged) {
      skip_block(file);
    }

    std::vector< int > isteps(npart);
    read_block(file, isteps);

    if (nums[0] >= 2) {
      // skip 2 blocks
      if (tagged) {
        skip_block(file);
      }
      skip_block(file);
    }

    std::string tag;
    if (tagged) {
      read_block(file, tag);
    } else {
      tag = "iphase";
    }

    if (tag != "iphase") {
      cmac_error("Wrong tag: \"%s\" (expected \"iphase\")!", tag.c_str());
    }

    std::vector< char > iphase(npart);
    read_block(file, iphase);

    //    std::vector<uint64_t> iunique(npart);
    if (nums[4] >= 1) {
      // skip iunique block
      if (tagged) {
        skip_block(file);
      }
      skip_block(file);
      //      if(tagged){
      //        skip_block(file);
      //      }
      //      read_block(file, iunique);
    }

    std::vector< double > x(npart);
    if (tagged) {
      skip_block(file);
    }
    read_block(file, x);

    std::vector< double > y(npart);
    if (tagged) {
      skip_block(file);
    }
    read_block(file, y);

    std::vector< double > z(npart);
    if (tagged) {
      skip_block(file);
    }
    read_block(file, z);

    std::vector< double > m(npart);
    if (tagged) {
      skip_block(file);
    }
    read_block(file, m);

    std::vector< double > h(npart);
    if (tagged) {
      skip_block(file);
    }
    read_block(file, h);

    // skip velocity, thermal energy and density blocks
    for (unsigned int i = 0; i < 5; ++i) {
      if (tagged) {
        skip_block(file);
      }
      skip_block(file);
    }

    // skip igrad related blocks
    for (unsigned int i = 0; i < nums[6] - 1; ++i) {
      if (tagged) {
        skip_block(file);
      }
      skip_block(file);
    }

    // skip sink particle data
    for (unsigned int i = 0; i < 10; ++i) {
      if (tagged) {
        skip_block(file);
      }
      skip_block(file);
    }

    for (unsigned int i = 0; i < npart; ++i) {
      if (iphase[i] == 0) {
        CoordinateVector<> rawunitsposition(x[i], y[i], z[i]);
        rawunitsbox.get_anchor() =
            CoordinateVector<>::min(rawunitsbox.get_anchor(), rawunitsposition);
        rawunitsbox.get_sides() =
            CoordinateVector<>::max(rawunitsbox.get_sides(), rawunitsposition);
        CoordinateVector<> position(x[i] * unit_length, y[i] * unit_length,
                                    z[i] * unit_length);
        _positions.push_back(position);
        _partbox.get_anchor() =
            CoordinateVector<>::min(_partbox.get_anchor(), position);
        _partbox.get_sides() =
            CoordinateVector<>::max(_partbox.get_sides(), position);
        _masses.push_back(m[i] * unit_mass);
        _smoothing_lengths.push_back(h[i] * unit_length);
        //        all_iunique.push_back(iunique[i]);
      }
    }
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

  _partbox.get_sides() -= _partbox.get_anchor();
  // add some margin to the box
  _partbox.get_anchor() -= 0.01 * _partbox.get_sides();
  _partbox.get_sides() *= 1.02;

  if (_log) {
    _log->write_status("Snapshot contains ", _positions.size(),
                       " gas particles.");
    _log->write_status(
        "Will create octree in box with anchor [", _partbox.get_anchor().x(),
        " m, ", _partbox.get_anchor().y(), " m, ", _partbox.get_anchor().z(),
        " m] and sides [", _partbox.get_sides().x(), " m, ",
        _partbox.get_sides().y(), " m, ", _partbox.get_sides().z(), " m]...");
    _log->write_info(
        "In raw units, this corresponds to a box with anchor [",
        rawunitsbox.get_anchor().x(), ", ", rawunitsbox.get_anchor().y(), ", ",
        rawunitsbox.get_anchor().z(), "], and sides [",
        rawunitsbox.get_sides().x(), ", ", rawunitsbox.get_sides().y(), ", ",
        rawunitsbox.get_sides().z(), "].");
  }

  if (!write_stats) {
    _stats_numbin = 0;
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
SPHNGSnapshotDensityFunction::SPHNGSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : SPHNGSnapshotDensityFunction(
          params.get_value< std::string >("DensityFunction:filename"),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "DensityFunction:initial temperature", "8000. K"),
          params.get_value< bool >("DensityFunction:write statistics", false),
          params.get_value< unsigned int >(
              "DensityFunction:statistics number of bins", 200),
          params.get_physical_value< QUANTITY_LENGTH >(
              "DensityFunction:statistics minimum distance", "1.e-5 m"),
          params.get_physical_value< QUANTITY_LENGTH >(
              "DensityFunction:statistics maximum distance", "1. kpc"),
          params.get_value< std::string >("DensityFunction:statistics filename",
                                          "ngb_statistics.txt"),
          params.get_value< bool >("DensityFunction:use new algorithm", false),
          log) {}

/**
 * @brief Destructor.
 *
 * Clean up the octree.
 */
SPHNGSnapshotDensityFunction::~SPHNGSnapshotDensityFunction() {
  delete _octree;
}

/**
 * @brief This routine constructs the internal Octree that is used for neighbour
 * finding.
 */
void SPHNGSnapshotDensityFunction::initialize() {
  _octree = new Octree(_positions, _partbox, false);
  _octree->set_auxiliaries(_smoothing_lengths, Octree::max< double >);

  if (_stats_numbin > 0) {
    if (_log) {
      _log->write_status("Obtaining particle neighbour statistics...");
    }
    const unsigned int numbin = _stats_numbin;
    const double mindist = _stats_mindist;
    const double maxdist = _stats_maxdist;
    const double logmindist = std::log10(mindist);
    const double logmaxdist = std::log10(maxdist);
    const double logddist = logmaxdist - logmindist;
    std::vector< double > ngbstats(numbin);
    double totnumngb = 0.;
    double numsmall = 0.;
    double numlarge = 0.;
    for (unsigned int i = 0; i < _positions.size(); ++i) {
      if (_log && i % (_positions.size() / 10) == 0) {
        _log->write_info("Got statistics for ", i, " of ", _positions.size(),
                         " particles.");
      }
      std::vector< unsigned int > ngbs = _octree->get_ngbs(_positions[i]);
      const unsigned int numngbs = ngbs.size();
      totnumngb += numngbs;
      for (unsigned int j = 0; j < numngbs; ++j) {
        unsigned int index = ngbs[j];
        if (index != i) {
          double r = (_positions[i] - _positions[index]).norm();
          if (r >= mindist) {
            unsigned int ibin =
                (std::log10(r) - logmindist) / logddist * numbin;
            if (ibin < numbin) {
              ngbstats[ibin] += 1.;
            } else {
              numlarge += 1.;
            }
          } else {
            numsmall += 1.;
          }
        }
      }
    }

    std::ofstream statfile(_stats_filename);
    double numinside = totnumngb - numsmall - numlarge;
    statfile << "# statistics account for " << (numinside / totnumngb) * 100.
             << "% of the particles.\n";
    statfile << "# " << (numsmall / totnumngb) * 100.
             << "% of the particles was closer together, "
             << (numlarge / totnumngb) * 100. << "% was further apart.\n";
    statfile << "#\n# r\tfraction\n";
    for (unsigned int i = 0; i < numbin; ++i) {
      double r = std::pow(10., i * logddist / numbin + logmindist);
      statfile << r << "\t" << ngbstats[i] / totnumngb << "\n";
    }
    if (_log) {
      _log->write_status("Wrote statistics file \"", _stats_filename, "\".");
    }
  }
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

double SPHNGSnapshotDensityFunction::full_integral(double phi, double r0,
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
    B3 = r03 / 4. * (-4. / 3. + (r0 / h) - 0.3 * r0h2 + 1. / 30. * r0h3 -
                     1. / 15. * r0h_3 + 8. / 5. * r0h_2);
    B2 = r03 / 4. * (-4. / 3. + (r0 / h) - 0.3 * r0h2 + 1. / 30. * r0h3 -
                     1. / 15. * r0h_3);
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
    I_5 =
        I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

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
    I_5 =
        I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

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
    I_5 =
        I_3 +
        a * (1. + a2) * (1. + a2) / 16. *
            ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

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
    full_int = r0h3 / M_PI * (1. / 6. * I_2 - 3. / 40. * r0h2 * I_4 +
                              1. / 40. * r0h3 * I_5 + B1 / r03 * I0);
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

double SPHNGSnapshotDensityFunction::mass_contribution(
    const Cell &cell, const CoordinateVector<> particle, const double h) {
  double M, Msum;

  Msum = 0.;
  M = 0.;

  std::vector< Face > face_vector = cell.get_faces();

  // Loop over each face of a cell.
  for (unsigned int i = 0; i < face_vector.size(); i++) {

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
 * @brief Function that gives the density for a given cell -> Maya.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues SPHNGSnapshotDensityFunction::operator()(const Cell &cell) const {

  DensityValues values;

  if (_use_new_algorithm) {

    CoordinateVector<> position = cell.get_cell_midpoint();

    // Find the vertex that is furthest away from the cell midpoint.
    std::vector< Face > face_vector = cell.get_faces();
    double radius = 0.0;
    for (unsigned int i = 0; i < face_vector.size(); i++) {
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
    std::vector< unsigned int > ngbs =
        _octree->get_ngbs_sphere(position, radius);
    const unsigned int numngbs = ngbs.size();

    double density = 0.;

    // Loop over all the neighbouring particles and calculate their mass
    // contributions.
    for (unsigned int i = 0; i < numngbs; i++) {
      const unsigned int index = ngbs[i];
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
    values.set_ionic_fraction(ION_He_n, 1.e-6);

  } else {

    const CoordinateVector<> position = cell.get_cell_midpoint();

    double density = 0.;
    std::vector< unsigned int > ngbs = _octree->get_ngbs(position);
    const unsigned int numngbs = ngbs.size();
    for (unsigned int i = 0; i < numngbs; ++i) {
      const unsigned int index = ngbs[i];
      const double r = (position - _positions[index]).norm();
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
    values.set_ionic_fraction(ION_He_n, 1.e-6);
  }

  return values;
}
