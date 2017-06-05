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
 * @file ExactGeometricTests.hpp
 *
 * @brief Exact geometric tests to test the orientation of a tetrahedron, and
 * to test if a point is inside the circumsphere of a tetrahedron.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EXACTGEOMETRICTESTS_HPP
#define EXACTGEOMETRICTESTS_HPP

#include <CoordinateVector.hpp>

#include <boost/multiprecision/integer.hpp>

/**
 * @brief Exact geometric tests to test the orientation of a tetrahedron, and
 * to test if a point is inside the circumsphere of a tetrahedron.
 */
class ExactGeometricTests {
private:
  /*! @brief Extended precision integer used for the exact tetrahedron
   *  orientation test. */
  typedef boost::multiprecision::int256_t int_orient3d;

  /*! @brief Extended precision integer used for the exact in sphere test. */
  typedef boost::multiprecision::number< boost::multiprecision::cpp_int_backend<
      278, 278, boost::multiprecision::signed_magnitude,
      boost::multiprecision::unchecked, void > >
      int_insphere;

public:
  /**
   * @brief Test the orientation of the tetrahedron that has the four given
   * points as vertices.
   *
   * The test returns a positive result if the fourth vertex is below the plane
   * through the three other vertices, with above the direction from which the
   * three points are ordered counterclockwise.
   *
   * E.g. if the four points are (0, 0, 0), (0, 0, 1), (0, 1, 0), and (1, 0, 0),
   * then this function returns 1.
   *
   * If the four points are exactly coplanar, then this function returns 0.
   *
   * @param a First vertex.
   * @param b Second vertex.
   * @param c Third vertex.
   * @param d Fourth vertex.
   * @return -1, 0, or 1, depending on the orientation of the tetrahedron.
   */
  inline static char orient3d(const CoordinateVector< unsigned long > &a,
                              const CoordinateVector< unsigned long > &b,
                              const CoordinateVector< unsigned long > &c,
                              const CoordinateVector< unsigned long > &d) {
    // the input coordinates should be the 53 bit mantissas of double precision
    // floating point values in the range [1, 2[
    const int_orient3d axp = a.x();
    const int_orient3d ayp = a.y();
    const int_orient3d azp = a.z();

    const int_orient3d bxp = b.x();
    const int_orient3d byp = b.y();
    const int_orient3d bzp = b.z();

    const int_orient3d cxp = c.x();
    const int_orient3d cyp = c.y();
    const int_orient3d czp = c.z();

    const int_orient3d dxp = d.x();
    const int_orient3d dyp = d.y();
    const int_orient3d dzp = d.z();

    // since all values above are positive, their differences can be at most 53
    // bits long as well
    const int_orient3d adx = axp - dxp;
    const int_orient3d ady = ayp - dyp;
    const int_orient3d adz = azp - dzp;

    const int_orient3d bdx = bxp - dxp;
    const int_orient3d bdy = byp - dyp;
    const int_orient3d bdz = bzp - dzp;

    const int_orient3d cdx = cxp - dxp;
    const int_orient3d cdy = cyp - dyp;
    const int_orient3d cdz = czp - dzp;

    // multiplication doubles the number of significant bits, so the numbers
    // below are at most 106 bits
    const int_orient3d bdxcdy = bdx * cdy;
    const int_orient3d cdxbdy = cdx * bdy;

    const int_orient3d cdxady = cdx * ady;
    const int_orient3d adxcdy = adx * cdy;

    const int_orient3d adxbdy = adx * bdy;
    const int_orient3d bdxady = bdx * ady;

    // since we do not know the sign of the terms between brackets, their
    // differences could gain a bit; the factors in between brackets are at most
    // 107 bits. Multiplication with a 53 bit factor means at most 160 bits for
    // the individual terms, and addition of three such terms can gain at most
    // 2 bits.
    // We hence need at most 162 bits of precision.
    const int_orient3d result = adz * (bdxcdy - cdxbdy) +
                                bdz * (cdxady - adxcdy) +
                                cdz * (adxbdy - bdxady);

    if (result > 0) {
      return 1;
    } else if (result < 0) {
      return -1;
    } else {
      return 0;
    }
  }

  /**
   * @brief Check if the fifth given point is inside the circumsphere of the
   * tetrahedron formed by the other four given points.
   *
   * It is assumed that the first four points are the vertices of a positively
   * oriented tetrahedron, as defined by a positive return value of orient3d().
   *
   * If the fifth point is exactly on the circumsphere of the tetrahedron, this
   * functions returns 0.
   *
   * @param a First vertex of the tetrahedron.
   * @param b Second vertex of the tetrahedron.
   * @param c Third vertex of the tetrahedron.
   * @param d Fourth vertex of the tetrahedron.
   * @param e Test point.
   * @return -1, 0, or 1, depending on the outcome of the geometric test.
   */
  inline static char insphere(const CoordinateVector< unsigned long > &a,
                              const CoordinateVector< unsigned long > &b,
                              const CoordinateVector< unsigned long > &c,
                              const CoordinateVector< unsigned long > &d,
                              const CoordinateVector< unsigned long > &e) {
    // the input coordinates should be the 53 bit mantissas of double precision
    // floating point values in the range [1, 2[
    const int_insphere axp = a.x();
    const int_insphere ayp = a.y();
    const int_insphere azp = a.z();

    const int_insphere bxp = b.x();
    const int_insphere byp = b.y();
    const int_insphere bzp = b.z();

    const int_insphere cxp = c.x();
    const int_insphere cyp = c.y();
    const int_insphere czp = c.z();

    const int_insphere dxp = d.x();
    const int_insphere dyp = d.y();
    const int_insphere dzp = d.z();

    const int_insphere exp = e.x();
    const int_insphere eyp = e.y();
    const int_insphere ezp = e.z();

    // since all values above are positive, the differences below can have at
    // most 53 significant bits
    const int_insphere aex = axp - exp;
    const int_insphere aey = ayp - eyp;
    const int_insphere aez = azp - ezp;

    const int_insphere bex = bxp - exp;
    const int_insphere bey = byp - eyp;
    const int_insphere bez = bzp - ezp;

    const int_insphere cex = cxp - exp;
    const int_insphere cey = cyp - eyp;
    const int_insphere cez = czp - ezp;

    const int_insphere dex = dxp - exp;
    const int_insphere dey = dyp - eyp;
    const int_insphere dez = dzp - ezp;

    // the terms in the expressions below have at most 106 significant bits
    // since we do not know the sign of the terms, the differences could gain a
    // bit, so the total expression has at most 107 significant bits
    const int_insphere ab = aex * bey - bex * aey;
    const int_insphere bc = bex * cey - cex * bey;
    const int_insphere cd = cex * dey - dex * cey;
    const int_insphere da = dex * aey - aex * dey;
    const int_insphere ac = aex * cey - cex * aey;
    const int_insphere bd = bex * dey - dex * bey;

    // the terms below have at most 160 bits; we can gain 2 bits through the
    // addition and get a total maximum of 162 significant bits
    const int_insphere abc = aez * bc - bez * ac + cez * ab;
    const int_insphere bcd = bez * cd - cez * bd + dez * bc;
    const int_insphere cda = cez * da + dez * ac + aez * cd;
    const int_insphere dab = dez * ab + aez * bd + bez * da;

    // the terms below have at most 106 significant bits, and the total
    // expression at most 108.
    const int_insphere aenrm2 = aex * aex + aey * aey + aez * aez;
    const int_insphere benrm2 = bex * bex + bey * bey + bez * bez;
    const int_insphere cenrm2 = cex * cex + cey * cey + cez * cez;
    const int_insphere denrm2 = dex * dex + dey * dey + dez * dez;

    // every term has at most 270 significant bits; the addition can gain 2
    // the total number of significant bits is capped at 272
    const int_insphere result =
        denrm2 * abc - cenrm2 * dab + benrm2 * cda - aenrm2 * bcd;

    if (result > 0) {
      return 1;
    } else if (result < 0) {
      return -1;
    } else {
      return 0;
    }
  }

  /**
   * @brief Get the 53 bit mantissa of the given double precision floating point
   * value.
   *
   * @param value Double precision floating point value.
   * @return Integer representation of the 53 bit mantissa.
   */
  inline static unsigned long get_mantissa(double value) {
    return value * (1l << DBL_MANT_DIG);
  }
};

#endif // EXACTGEOMETRICTESTS_HPP
