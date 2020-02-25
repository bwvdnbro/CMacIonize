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
 * @file AsciiFileTablePhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution that reads photon sources from an Ascii file.
 *
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#ifndef ASCIIFILETABLEPHOTONSOURCEDISTRIBUTION_HPP
#define ASCIIFILETABLEPHOTONSOURCEDISTRIBUTION_HPP

#include "CoordinateVector.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;

/**
 * @brief PhotonSourceDistribution that reads photon sources from an Ascii file.
 */
class AsciiFileTablePhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Positions of the sources in the ASCII file (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Luminosity weights of the sources in the ASCII file. */
  std::vector< double > _weights;

  /*! @brief Total luminosity of all sources in the ASCII file (in s^-1). */
  double _total_luminosity;

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  AsciiFileTablePhotonSourceDistribution(std::string filename, Log *log);
  AsciiFileTablePhotonSourceDistribution(ParameterFile &params, Log *log);

  virtual photonsourcenumber_t get_number_of_sources() const;
  virtual CoordinateVector<> get_position(photonsourcenumber_t index);
  virtual double get_weight(photonsourcenumber_t index) const;
  virtual double get_total_luminosity() const;
};

#endif // ASCIIFILETABLEPHOTONSOURCEDISTRIBUTION_HPP
