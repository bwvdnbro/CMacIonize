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
 * @file GadgetSnapshotPhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution that reads photon sources from a Gadget2 type
 * 3 snapshot file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef GADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_HPP
#define GADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_HPP

#include "PhotonSourceDistribution.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;

/**
 * @brief PhotonSourceDistribution that reads photon sources from a Gadget2 type
 * 3 snapshot file.
 */
class GadgetSnapshotPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Positions of the sources in the snapshot file (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  GadgetSnapshotPhotonSourceDistribution(std::string filename, Log *log = NULL);
  GadgetSnapshotPhotonSourceDistribution(ParameterFile &params,
                                         Log *log = NULL);

  virtual unsigned int get_number_of_sources();
  virtual CoordinateVector<> get_position(unsigned int index);
  virtual double get_weight(unsigned int index);
};

#endif // GADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_HPP
