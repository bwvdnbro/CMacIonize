/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SingleSupernovaPhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution without UV sources that mimicks a single
 * supernova explosion with an adjustable energy output.
 *
 * Used to test the hydrodynamical impact of a single supernova explosion on
 * an ISM with some equation of state.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SINGLESUPERNOVAPHOTONSOURCEDISTRIBUTION_HPP
#define SINGLESUPERNOVAPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"

/**
 * @brief PhotonSourceDistribution without UV sources that mimicks a single
 * supernova explosion with an adjustable energy output.
 */
class SingleSupernovaPhotonSourceDistribution
    : public PhotonSourceDistribution {
private:
  /*! @brief Position of the supernova explosion (in m). */
  const CoordinateVector<> _position;

  /*! @brief Life time of the source before it explodes (in s). */
  const double _lifetime;

  /*! @brief Luminosity of the source before it goes supernova (in s^-1). */
  const double _luminosity;

  /*! @brief Energy of the supernova exposion (in J). */
  const double _energy;

  /*! @brief Flag signalling if the supernova already exploded. */
  bool _has_exploded;

public:
  /**
   * @brief Constructor.
   *
   * @param position Position of the supernova explosion (in m).
   * @param lifetime Life time of the source before it explodes (in s).
   * @param luminosity Ionising luminosity of the source before it goes
   * supernova (in s^-1).
   * @param energy Energy of the supernova explosion (in J).
   * @param log Log to write logging info to.
   */
  inline SingleSupernovaPhotonSourceDistribution(
      const CoordinateVector<> position, const double lifetime,
      const double luminosity, const double energy, Log *log = nullptr)
      : _position(position), _lifetime(lifetime), _luminosity(luminosity),
        _energy(energy), _has_exploded(false) {

    if (log != nullptr) {
      log->write_status(
          "Set up SingleSupernovaPhotonSourceDistribution at position [",
          _position.x(), " m, ", _position.y(), " m, ", _position.z(),
          " m], with life time ", _lifetime, " s, ionising luminosity ",
          _luminosity, " s^-1 and with SN energy ", _energy, " J.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - position: Position of the supernova explosion (default: [0. m, 0. m, 0.
   *    m])
   *  - lifetime: Life time of the source before it explodes (default: 10. Myr)
   *  - luminosity: Ionising luminosity of the source before it goes supernova
   *    (default: 1.e49 s^-1)
   *  - energy: Energy of the supernova explosion (default: 1.e51 erg)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SingleSupernovaPhotonSourceDistribution(ParameterFile &params,
                                          Log *log = nullptr)
      : SingleSupernovaPhotonSourceDistribution(
            params.get_physical_vector< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:position", "[0. m, 0. m, 0. m]"),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:lifetime", "10. Myr"),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "PhotonSourceDistribution:luminosity", "1.e49 s^-1"),
            params.get_physical_value< QUANTITY_ENERGY >(
                "PhotonSourceDistribution:energy", "1.e51 erg"),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~SingleSupernovaPhotonSourceDistribution() {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * The PhotonSourceDistribution will return exactly this number of valid
   * and unique positions by successive application of operator().
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_sources() const {
    return _luminosity > 0. && !_has_exploded;
  }

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {
    return _position;
  }

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the photon source, used to determine how many photons are
   * emitted from this particular source.
   */
  virtual double get_weight(photonsourcenumber_t index) const { return 1.; }

  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */
  virtual double get_total_luminosity() const {
    if (!_has_exploded) {
      return _luminosity;
    } else {
      return 0.;
    }
  }

  /**
   * @brief Update the distribution after the system moved to the given time.
   *
   * @param simulation_time Current simulation time (in s).
   * @return True if the distribution changed, false otherwise.
   */
  virtual bool update(const double simulation_time) {

    if (!_has_exploded && simulation_time >= _lifetime) {
      _has_exploded = true;
      return true;
    } else {
      // make sure the PhotonSource is warned if the number of sources is
      // zero
      return _has_exploded;
    }
  }

  /**
   * @brief Add stellar feedback at the given time.
   *
   * @param grid DensityGrid to operate on.
   * @param time Current simulation time (in s).
   */
  virtual void add_stellar_feedback(DensityGrid &grid, const double time) {

    if (!_has_exploded && time >= _lifetime) {
      DensityGrid::iterator cell = grid.get_cell(_position);
      cell.get_hydro_variables().set_energy_term(_energy);
      _has_exploded = true;
    }
  }
};

#endif // SINGLESUPERNOVAPHOTONSOURCEDISTRIBUTION_HPP
