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
class UVLuminosityFunction;

/**
 * @brief PhotonSourceDistribution that reads photon sources from a Gadget2 type
 * 3 snapshot file.
 */
class GadgetSnapshotPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Positions of the sources in the snapshot file (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief UV luminosities of the sources in the snapshot file (in s^-1). */
  std::vector< double > _luminosities;

  /*! @brief Total luminosity of all sources in the snapshot file (in s^-1). */
  double _total_luminosity;

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  GadgetSnapshotPhotonSourceDistribution(
      std::string filename, std::string formation_time_name, const Box<> box,
      UVLuminosityFunction *luminosity_function,
      double fallback_unit_length_in_SI = 0.,
      double fallback_unit_time_in_SI = 0.,
      double fallback_unit_mass_in_SI = 0., double cutoff_age = 1.577e14,
      bool use_gas = false, double SFR_unit = 0.,
      bool comoving_integration = false, double hubble_parameter = 0.7,
      Log *log = nullptr);
  GadgetSnapshotPhotonSourceDistribution(ParameterFile &params,
                                         Log *log = nullptr);

  virtual photonsourcenumber_t get_number_of_sources() const;
  virtual CoordinateVector<> get_position(photonsourcenumber_t index);
  virtual double get_weight(photonsourcenumber_t index) const;
  virtual double get_total_luminosity() const;
};

#endif // GADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_HPP
