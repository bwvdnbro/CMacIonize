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

  /*! @brief Fraction of the total luminosity that comes from the bulge. */
  const double _bulge_to_total;

  /*! @brief Constant used to sample from the bulge. */
  const double _rB_over_rJ_plus_rB;

  /*! @brief Cumulative disc luminosity distribution: x coordinates (in m). */
  std::vector< double > _cumulative_distribution_x;

  /*! @brief Cumulative disc luminosity distribution: y coordinates. */
  std::vector< double > _cumulative_distribution_y;

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
      lbulge = 0.25 * _lb0 / (M_PI * u * u * up1 * up1);
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
   * @param log Log to write logging info to.
   */
  SpiralGalaxyContinuousPhotonSource(
      const Box<> box, const CoordinateVector< unsigned int > ncell,
      double rstars, double hstars, double B_over_T, Log *log = nullptr)
      : _box(box), _ncell(ncell), _kpc(3.086e19), _rbulge(2. * _kpc),
        _r_jaffe(0.4 * _kpc),
        _lb0(B_over_T / (_r_jaffe * _r_jaffe * _r_jaffe *
                         (1. - 1. / (1. + _rbulge / _r_jaffe)))),
        _rstars(rstars), _hstars(hstars),
        _ld0(0.25 * (1. - B_over_T) / (M_PI * _hstars * _rstars * _rstars)),
        _bulge_to_total(B_over_T),
        _rB_over_rJ_plus_rB(_rbulge / (_rbulge + _r_jaffe)) {

    _cumulative_distribution_x.resize(1000, 0.);
    _cumulative_distribution_y.resize(1000, 0.);
    for (unsigned int i = 0; i < 999; ++i) {
      const double w = i * 15. * _kpc * 0.001;
      _cumulative_distribution_x[i] = w;
      const double x = w / _rstars;
      _cumulative_distribution_y[i] = 1. - (1. + x) * std::exp(-x);
    }
    _cumulative_distribution_x[999] = 15. * _kpc;
    _cumulative_distribution_y[999] = 1.;

    const unsigned int ntot = ncell.x() * ncell.y() * ncell.z() + 1;
    _luminosities.resize(ntot);
    unsigned int index = 0;
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
          ++index;
        }
      }
    }
    _luminosities[index] = 0.;
    ++index;

    cmac_assert(index == ntot);

    // make cumulative
    for (unsigned int i = 1; i < ntot; ++i) {
      _luminosities[i] += _luminosities[i - 1];
    }

    // normalize
    for (unsigned int i = 0; i < ntot; ++i) {
      _luminosities[i] /= _luminosities.back();
    }

    if (log) {
      log->write_status("Constructed SpiralGalaxyContinuousPhotonSource in a "
                        "box with anchor [",
                        _box.get_anchor().x(), " m, ", _box.get_anchor().y(),
                        " m, ", _box.get_anchor().z(), " m] and sides [",
                        _box.get_sides().x(), " m, ", _box.get_sides().y(),
                        " m, ", _box.get_sides().z(),
                        " m], with a disc scale length of ", _rstars,
                        " m, and a disc scale height of ", _hstars, " m.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SpiralGalaxyContinuousPhotonSource(ParameterFile &params, Log *log = nullptr)
      : SpiralGalaxyContinuousPhotonSource(
            Box<>(
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_anchor", "[-12. kpc, -12. kpc, -12. kpc]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid:box_sides", "[24. kpc, 24. kpc, 24. kpc]")),
            params.get_value< CoordinateVector< unsigned int > >(
                "continuousphotonsource:ncell",
                CoordinateVector< unsigned int >(201)),
            params.get_physical_value< QUANTITY_LENGTH >(
                "continuousphotonsource:r_stars", "5. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "continuousphotonsource:h_stars", "0.6 kpc"),
            params.get_value< double >("continuousphotonsource:B_over_T", 0.2),
            log) {}

  /**
   * @brief Get the position and direction of a random photon.
   *
   * We randomly sample a position from the combined luminosity distribution.
   * First, we choose whether we need to sample the disc or the bulge. This is
   * easy, as we know the bulge to total ratio. We hence just need to draw a
   * random number, and if this number is smaller that the bulge to total ratio,
   * we sample from the bulge. If it is larger, we sample from the disc.
   *
   * The bulge luminosity is given by
   * @f[
   *   L_B(r) = \frac{L_{B,0}}{\left(\frac{r}{r_J}\right)^2 \left(1 +
   *   \frac{r}{r_J}\right)^2},
   * @f]
   * with @f$r_J@f$ a scale parameter. The constant @f$L_{B,0}@f$ is set so that
   * the total integrated luminosity within some cutoff radius @f$r_B@f$ is
   * unity. The integrated luminosity within a radius @f$r@f$ is given by
   * @f[
   *   L_B(<r) = 4\pi{}L_{B,0}r_J^3\left(1-\frac{1}{1+\frac{r}{r_J}}\right),
   * @f]
   * and hence
   * @f[
   *   L_{B,0} = \frac{1+\frac{r_B}{r_J}}{4\pi{}r_J^2r_B}.
   * @f]
   *
   * If @f$u@f$ is a random uniform number in the range [0,1], then we can find
   * a radius @f$r@f$ distributed according to the bulge luminosity distribution
   * by inverting the cumulative distribution @f$L_B(<r)@f$:
   * @f[
   *   r = \frac{r_J}{1 - \left(\frac{r_B}{r_J+r_B}\right)u} - r_J.
   * @f]
   *
   * The disc luminosity is given by
   * @f[
   *   L_D(w,z) = L_{D,0} e^{-\frac{|z|}{h}} e^{-\frac{w}{w_s}},
   * @f]
   * with @f$h@f$ the disc scale height, @f$w = \sqrt{x^2 + y^2}@f$ the radius
   * in cylindrical coordinates, and @f$w_s@f$ the scale radius of the disc.
   * Again, the constant @f$L_{D,0}@f$ is chosen so that the total luminosity
   * within the disc is unity (where we integrate out to infinity). The
   * total integrated luminosity within a cylindrical radius @f$w@f$ and height
   * @f$z@f$ is given by
   * @f[
   *   L_D(<w,<|z|) = 4\pi{}hw_s^2 L_{D,0} \left(1 - e^{-\frac{|z|}{h}}\right)
   *   \left(1 - \left(1 + \frac{w}{w_s}\right)e^{-\frac{w}{w_s}} \right),
   * @f]
   * which evaluated at infinity gives us
   * @f[
   *   L_{D,0} = \frac{1}{4\pi{}hw_s^2}.
   * @f]
   *
   * It is important to note that the total luminosity distribution is in fact
   * the product of two separate luminosity distributions:
   * @f[
   *   L_D(w,z) = L_{D,w}(w) L_{D,z}(z),
   * @f]
   * with
   * @f[
   *   L_{D,w}(w) = \frac{1}{2\pi{}w_s^2} e^{-\frac{w}{w_s}},
   * @f]
   * and
   * @f[
   *   L_{D,z}(z) = \frac{1}{2h} e^{-\frac{|z|}{h}}.
   * @f]
   * If we want to draw a random position from this combined distribution, we
   * can sample the polar coordinate and the height separately.
   *
   * The height distribution can be easily inverted: if @f$u@f$ is a random
   * uniform number in the range [0,1], then
   * @f[
   *   |z| = -h \ln{u}.
   * @f]
   * To get the sign of @f$z@f$, we can sample a second random number.
   *
   * The polar distribution is somewhat more complicated, as the cumulative
   * polar distribution function
   * @f[
   *   L_{D,w}(<w) = 1 - \left(1 + \frac{w}{w_s}\right)e^{-\frac{w}{w_s}}
   * @f]
   * cannot be inverted analytically. We hence pretabulate the cumulative
   * distribution function, and obtain a random cylindrical radius from that
   * table.
   *
   * @param random_generator RandomGenerator to use.
   * @return std::pair of CoordinateVector instances, specifying a starting
   * position (in m) and direction for an incoming photon.
   */
  virtual std::pair< CoordinateVector<>, CoordinateVector<> >
  get_random_incoming_direction(RandomGenerator &random_generator) const {

    // initialize as a point outside the box
    CoordinateVector<> position = _box.get_anchor() - _box.get_sides();

    // we need to check if the sampled coordinate is within the box
    // if not, we need to repeat the sampling procedure until it is
    while (!_box.inside(position)) {
      // choose if we sample from bulge or disc
      const double x_bulge = random_generator.get_uniform_random_double();
      if (x_bulge <= _bulge_to_total) {
        // sample from bulge
        const double u = random_generator.get_uniform_random_double();
        const double r = _r_jaffe / (1. - _rB_over_rJ_plus_rB * u) - _r_jaffe;
        const double phi =
            2. * M_PI * random_generator.get_uniform_random_double();
        const double cost =
            2. * random_generator.get_uniform_random_double() - 1.;
        const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
        position[0] = r * sint * std::cos(phi);
        position[1] = r * sint * std::sin(phi);
        position[2] = r * cost;
      } else {
        // sample from disc
        const double u1 =
            2. * random_generator.get_uniform_random_double() - 1.;
        double z;
        if (u1 > 0.) {
          z = -_hstars * std::log(u1);
        } else {
          z = _hstars * std::log(-u1);
        }
        const double phi =
            2. * M_PI * random_generator.get_uniform_random_double();
        const double u2 = random_generator.get_uniform_random_double();
        unsigned int i = Utilities::locate(u2, &_cumulative_distribution_y[0],
                                           _cumulative_distribution_y.size());
        const double w =
            (u2 - _cumulative_distribution_y[i]) /
            (_cumulative_distribution_y[i + 1] -
             _cumulative_distribution_y[i]) *
            (_cumulative_distribution_x[i + 1] - _cumulative_distribution_x[i]);
        position[0] = w * std::cos(phi);
        position[1] = w * std::sin(phi);
        position[2] = z;
      }
    }
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
