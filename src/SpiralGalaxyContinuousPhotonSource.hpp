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
 *
 * The luminosity profile is a combination of the bulge model of Jaffe, W. 1983,
 * MNRAS, 202, 995 (http://adsabs.harvard.edu/abs/1983MNRAS.202..995J), and a
 * double exponential profile.
 *
 * We manually computed and inverted the cumulative distributions for both
 * models to obtain a random position generating function that can be used in
 * parallel.
 */
class SpiralGalaxyContinuousPhotonSource : public ContinuousPhotonSource {
private:
  /*! @brief Box containing the luminosity distribution grid (in m). */
  const Box<> _box;

  /*! @brief Conversion factor from kpc to m (in m). */
  const double _kpc;

  /*! @brief Inner cutoff radius of the bulge (in m). */
  const double _rC;

  /*! @brief Outer cutoff radius of the bulge (in m). */
  const double _rB;

  /*! @brief Scale radius of the bulge (in m). */
  const double _rJ;

  /*! @brief Scale length @f$r_*@f$ of the stellar disk (in m). */
  const double _r_stars;

  /*! @brief Scale height @f$h_*@f$ of the stellar disk (in m). */
  const double _h_stars;

  /*! @brief Fraction of the total luminosity that comes from the bulge. */
  double _bulge_to_total_ratio;

  /*! @brief Constant used to sample from the bulge. */
  const double _rB_over_rJ_plus_rB;

  /*! @brief Constant used to sample from the bulge. */
  const double _rC_over_rJ_plus_rC;

  /*! @brief Cumulative disc luminosity distribution: x coordinates (in m). */
  std::vector< double > _cumulative_disc_luminosity_distribution_x;

  /*! @brief Cumulative disc luminosity distribution: y coordinates. */
  std::vector< double > _cumulative_disc_luminosity_distribution_y;

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the luminosity distribution grid (in m).
   * @param rstars Scale length @f$r_*@f$ of the stellar disk (in m).
   * @param hstars Scale height @f$h_*@f$ of the stellar disk (in m).
   * @param B_over_T Ratio of the bulge luminosity and the total luminosity,
   * @f$B/T@f$.
   * @param log Log to write logging info to.
   */
  SpiralGalaxyContinuousPhotonSource(const Box<> box, double rstars,
                                     double hstars, double B_over_T,
                                     Log *log = nullptr)
      : _box(box), _kpc(3.086e19), _rC(0.2 * _kpc), _rB(2. * _kpc),
        _rJ(0.4 * _kpc), _r_stars(rstars), _h_stars(hstars),
        _bulge_to_total_ratio(B_over_T), _rB_over_rJ_plus_rB(_rB / (_rB + _rJ)),
        _rC_over_rJ_plus_rC(_rC / (_rC + _rJ)) {

    // the total bulge luminosity is 4\pi L_{B,0} r_j^3 (\frac{r_b}{r_j+r_b})
    // however, we cut off the bulge in the inner r_c = 0.2 kpc, removing a part
    // 4 \pi L_{B,0} r_j^3(\frac{r_c}{r_j+r_c})
    // this means that the part of the bulge we sample has a lower luminosity,
    // and the sampled bulge to total needs to be adjusted
    _bulge_to_total_ratio *= (1. - _rC_over_rJ_plus_rC / _rB_over_rJ_plus_rB);

    // we assume the box is centered on (0., 0., 0.)
    const double rmax = 1.2 * _box.get_anchor().norm();
    const unsigned int nbin = 1000;
    _cumulative_disc_luminosity_distribution_x.resize(nbin + 1, 0.);
    _cumulative_disc_luminosity_distribution_y.resize(nbin + 1, 0.);
    for (unsigned int i = 0; i < nbin; ++i) {
      const double w = i * rmax / nbin;
      _cumulative_disc_luminosity_distribution_x[i] = w;
      const double x = w / _r_stars;
      _cumulative_disc_luminosity_distribution_y[i] =
          1. - (1. + x) * std::exp(-x);
    }
    // we make sure the maximum of our distribution is outside the box, so that
    // positions sampled in that bin always end up being thrown away...
    _cumulative_disc_luminosity_distribution_x[nbin] = rmax;
    _cumulative_disc_luminosity_distribution_y[nbin] = 1.;

    if (log) {
      log->write_status("Constructed SpiralGalaxyContinuousPhotonSource with a "
                        "disc scale length of ",
                        _r_stars, " m, a disc scale height of ", _h_stars,
                        " m, and a bulge to total ratio of ", B_over_T, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - scale length stars: Scale length of the stellar disc (default: 5. kpc)
   *  - scale height stars: Scale height of the stellar disc (default: 0.6 kpc)
   *  - bulge over total ratio: Ratio of the bulge luminosity and the total
   *    luminosity (default: 0.2)
   *
   * @param simulation_box Simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SpiralGalaxyContinuousPhotonSource(const Box<> &simulation_box,
                                     ParameterFile &params, Log *log = nullptr)
      : SpiralGalaxyContinuousPhotonSource(
            simulation_box,
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:scale length stars", "5. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ContinuousPhotonSource:scale height stars", "0.6 kpc"),
            params.get_value< double >(
                "ContinuousPhotonSource:bulge over total ratio", 0.2),
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
   * If @f$u@f$ is a random uniform number in the range @f$[0,1]@f$, then we can
   * find a radius @f$r@f$ distributed according to the bulge luminosity
   * distribution by inverting the cumulative distribution @f$L_B(<r)@f$:
   * @f[
   *   r = \frac{r_J}{\frac{1}{\left(\frac{r_B}{r_J+r_B}\right)u}-1}.
   * @f]
   *
   * For the specific case here, we also use an inner cutoff radius @f$r_C@f$.
   * The integrated luminosity within some radius @f$r \geq{} r_C@f$ is
   * @f[
   *   L_B(r_C \leq{} r) = 4\pi{}L_{B,0}r_J^3\left(\frac{r}{r+r_J} -
   *   \frac{r_C}{r_C + r_J} \right).
   * @f]
   * Still requiring that the total luminosity within @f$[r_C, r_B]@f$ is unity,
   * we get a new value for the constant @f$L_{B,0}@f$:
   * @f[
   *   L_{B,0} = \frac{1}{4\pi{}r_J^3}\frac{1}{\frac{r_B}{r_J+r_B} -
   *   \frac{r_C}{r_J+r_C}} = \frac{1}{4\pi{}r_J^3}\frac{1}{r_{BJ} - r_{CJ}},
   * @f]
   * where we introduced @f$r_{BJ} = \frac{r_B}{r_J+r_B}@f$ and @f$r_{CJ} =
   * \frac{r_C}{r_J+r_C}@f$.
   *
   * Plugging this into the expression for the integrated luminosity, we get
   * @f[
   *   L_B(r_C \leq{} r) = \frac{1}{r_{BJ}-r_{CJ}} \left(\frac{r}{r+r_J} -
   *   r_{CJ} \right).
   * @f]
   * If @f$u@f$ is a random uniform number in the range @f$[0,1]@f$, then
   * @f[
   *   r = \frac{r_J}{\frac{1}{r_{BJ}u + (1-u)r_{CJ}}-1}
   * @f]
   * is distributed according to the truncated bulge luminosity function.
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
   * uniform number in the range @f$[0,1]@f$, then
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
      if (x_bulge <= _bulge_to_total_ratio) {
        // sample from bulge
        const double u = random_generator.get_uniform_random_double();
        const double A =
            u * _rB_over_rJ_plus_rB + (1. - u) * _rC_over_rJ_plus_rC;
        const double r = _rJ / (1. / A - 1.);
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
          z = -_h_stars * std::log(u1);
        } else {
          z = _h_stars * std::log(-u1);
        }
        const double phi =
            2. * M_PI * random_generator.get_uniform_random_double();
        const double u2 = random_generator.get_uniform_random_double();
        unsigned int i = Utilities::locate(
            u2, &_cumulative_disc_luminosity_distribution_y[0],
            _cumulative_disc_luminosity_distribution_y.size());
        const double w =
            _cumulative_disc_luminosity_distribution_x[i] +
            (u2 - _cumulative_disc_luminosity_distribution_y[i]) /
                (_cumulative_disc_luminosity_distribution_y[i + 1] -
                 _cumulative_disc_luminosity_distribution_y[i]) *
                (_cumulative_disc_luminosity_distribution_x[i + 1] -
                 _cumulative_disc_luminosity_distribution_x[i]);
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
