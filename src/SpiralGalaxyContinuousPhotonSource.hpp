/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SpiralGalaxyContinuousPhotonSource.hpp
 *
 * @brief ContinuousPhotonSource implementation that represents a spiral galaxy.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPIRALGALAXYCONTINUOUSPHOTONSOURCE_HPP
#define SPIRALGALAXYCONTINUOUSPHOTONSOURCE_HPP

#include "Box.hpp"
#include "ContinuousPhotonSource.hpp"
#include "Error.hpp"
#include "PhotonSource.hpp"
#include "RandomGenerator.hpp"
#include "Utilities.hpp"

#include <vector>

/**
 * @brief ContinuousPhotonSource implementation that represents a spiral galaxy.
 */
class SpiralGalaxyContinuousPhotonSource : public ContinuousPhotonSource {
private:
  /*! @brief Box containing the luminosity distribution grid (in m). */
  const Box<> _box;

  /*! @brief Number of cells in each dimension. */
  const CoordinateVector< unsigned int > _ncell;

  /*! @brief Luminosity distribution grid. */
  std::vector< double > _luminosities;

  /*! @brief Conversion factor from kpc to m (in m). */
  const double _kpc;

  /*! @brief Cutoff radius of the bulge (in m). */
  const double _rbulge;

  /*! @brief Scale radius of the bulge (in m). */
  const double _r_jaffe;

  /*! @brief Luminosity of the bulge (no unit, as we are only interested in the
   *  distribution function). */
  const double _lb0;

  /*! @brief Scale length @f$r_*@f$ of the stellar disk (in m). */
  const double _rstars;

  /*! @brief Scale height @f$h_*@f$ of the stellar disk (in m). */
  const double _hstars;

  /*! @brief Luminosity of the disk (no unit, as we are only interested in the
   *  distribution function). */
  const double _ld0;

  /**
   * @brief Get the luminosity value at the given position.
   *
   * @param position Position (in m).
   * @return Luminosity (no unit, as we are only interested in the distribution
   * function).
   */
  inline double get_luminosity(const CoordinateVector<> position) const {
    const double w2 = position.x() * position.x() + position.y() * position.y();
    const double w = std::sqrt(w2);
    const double r2 = w2 + position.z() * position.z();
    const double r = std::sqrt(r2);
    double lbulge = 0.;
    if (r < _rbulge && r > 0.2 * _kpc) {
      const double u = r / _r_jaffe;
      const double up1 = 1. + u;
      lbulge = 0.25 * _lb0 / M_PI / (u * u * up1 * up1);
    }
    double ldisk = 0.;
    if (r < 15. * _kpc) {
      ldisk = _ld0 * std::exp(-std::abs(position.z()) / _hstars) *
              std::exp(-w / _rstars);
    }
    return ldisk + lbulge;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the luminosity distribution grid (in m).
   * @param ncell Number of cells in each dimension.
   * @param rstars Scale length @f$r_*@f$ of the stellar disk (in m).
   * @param hstars Scale height @f$h_*@f$ of the stellar disk (in m).
   * @param B_over_T Ratio of the bulge luminosity and the total luminosity,
   * @f$B/T@f$.
   */
  SpiralGalaxyContinuousPhotonSource(
      const Box<> box, const CoordinateVector< unsigned int > ncell,
      double rstars, double hstars, double B_over_T)
      : _box(box), _ncell(ncell), _kpc(3.086e19), _rbulge(2. * _kpc),
        _r_jaffe(0.4 * _kpc),
        _lb0(B_over_T / (_r_jaffe * _r_jaffe * _r_jaffe *
                         (1. - 1. / (1. + _rbulge / _r_jaffe)))),
        _rstars(rstars), _hstars(hstars),
        _ld0(0.25 * (1. - B_over_T) / (M_PI * _hstars * _rstars * _rstars)) {

    const unsigned int ntot = ncell.x() * ncell.y() * ncell.z() + 1;
    _luminosities.resize(ntot);
    unsigned int index = 0;
    double total_luminosity = 0.;
    for (unsigned int ix = 0; ix < ncell.x(); ++ix) {
      for (unsigned int iy = 0; iy < ncell.y(); ++iy) {
        for (unsigned int iz = 0; iz < ncell.z(); ++iz) {
          const CoordinateVector<> position(
              box.get_anchor().x() +
                  (ix + 0.5) * box.get_sides().x() / ncell.x(),
              box.get_anchor().y() +
                  (iy + 0.5) * box.get_sides().y() / ncell.y(),
              box.get_anchor().z() +
                  (iz + 0.5) * box.get_sides().z() / ncell.z());

          _luminosities[index] = get_luminosity(position);
          total_luminosity += _luminosities[index];
          ++index;
        }
      }
    }
    _luminosities[index] = 0.;
    ++index;

    cmac_assert(index == ntot);

    // normalize distribution function and make cumulative
    _luminosities[0] /= total_luminosity;
    for (unsigned int i = 1; i < ntot; ++i) {
      _luminosities[i] /= total_luminosity;
      _luminosities[i] += _luminosities[i - 1];
    }
  }

  /**
   * @brief Get the position and direction of a random photon.
   *
   * @param random_generator RandomGenerator to use.
   * @return std::pair of CoordinateVector instances, specifying a starting
   * position (in m) and direction for an incoming photon.
   */
  virtual std::pair< CoordinateVector<>, CoordinateVector<> >
  get_random_incoming_direction(RandomGenerator &random_generator) const {

    // find out from which cell we emit
    const double x = random_generator.get_uniform_random_double();
    unsigned int active_index =
        Utilities::locate(x, &_luminosities[0], _luminosities.size());

    // get the index of the cell
    const unsigned int ix = active_index / _ncell.y() / _ncell.z();
    active_index -= ix * _ncell.y() * _ncell.z();
    const unsigned int iy = active_index / _ncell.z();
    active_index -= iy * _ncell.z();
    const unsigned int iz = active_index;

    // draw a uniform random position from that cell
    const CoordinateVector<> position(
        _box.get_anchor().x() +
            (ix + random_generator.get_uniform_random_double()) *
                _box.get_sides().x() / _ncell.x(),
        _box.get_anchor().y() +
            (iy + random_generator.get_uniform_random_double()) *
                _box.get_sides().y() / _ncell.y(),
        _box.get_anchor().z() +
            (iz + random_generator.get_uniform_random_double()) *
                _box.get_sides().z() / _ncell.z());
    // get a random isotropic direction
    const CoordinateVector<> direction =
        PhotonSource::get_random_direction(random_generator);

    return std::make_pair(position, direction);
  }

  /**
   * @brief Get the total surface area of the emitting volume.
   *
   * This is used to calculate the total luminosity of the radiation.
   *
   * @return Total surface area of the emitting volume (in m^2).
   */
  virtual double get_total_surface_area() const { return 1.; }
};

#endif // SPIRALGALAXYCONTINUOUSPHOTONSOURCE_HPP
