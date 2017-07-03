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

#include "Configuration.hpp"
#include "CoordinateVector.hpp"

#ifdef HAVE_MULTIPRECISION
#include <boost/multiprecision/integer.hpp>
#endif

/**
 * @brief Exact geometric tests to test the orientation of a tetrahedron, and
 * to test if a point is inside the circumsphere of a tetrahedron.
 */
class ExactGeometricTests {
#ifdef HAVE_MULTIPRECISION
private:
  /*! @brief Extended precision integer used for the exact tetrahedron
   *  orientation test. */
  typedef boost::multiprecision::int256_t int_orient3d;

  /*! @brief Extended precision integer used for the exact in sphere test. */
  typedef boost::multiprecision::number< boost::multiprecision::cpp_int_backend<
      278, 278, boost::multiprecision::signed_magnitude,
      boost::multiprecision::unchecked, void > >
      int_insphere;
#endif

public:
  /**
   * @brief Auxiliary typedef used to extract the mantissa from a double
   * precision floating point value.
   *
   * The variables in the union occupy the same memory, which allows us to
   * access the bytes used to store the double precision floating point value.
   */
  typedef union {
    /*! @brief Double precision floating point value. */
    double dvalue;
    /**
     * @brief Anonymous struct containing the 3 parts of a general double
     * precision floating point value.
     *
     * A general double precision floating point value has a sign @f$s@f$, a
     * mantissa @f$m@f$, and an exponent @f$e@f$, so that the value @f$v@f$ of
     * the double is given by
     * \f[
     * v = s \times 1.m \times 2^{e-1023}
     * \f]
     * (see
     * https://en.wikipedia.org/wiki/Double-precision_floating-point_format).
     *
     * The IEEE 754 standard specifies that @f$s@f$, @f$m@f$ and @f$e@f$ have
     * respectively 1 bit, 52 bit and 11 bit precision. In memory (low to high
     * bits), they are ordered as follows:
     * \f[
     * m e s
     * \f]
     */
    struct {
      /*! @brief Mantissa @f$m@f$. */
      unsigned long mantissa : 52;
      /*! @brief Exponent @f$e@f$. */
      unsigned long exponent : 11;
      /*! @brief Sign @f$s@f$. */
      unsigned long sign : 1;
    } parts;
  } binary_double;

  /**
   * @brief Get the 52 bit mantissa of the given double precision floating point
   * value.
   *
   * @param value Double precision floating point value.
   * @return 52 bit mantissa of that same value.
   */
  inline static unsigned long get_mantissa(double value) {
    binary_double dvalue;
    dvalue.dvalue = value;
    return dvalue.parts.mantissa;
  }

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
  inline static char orient3d_exact(const CoordinateVector< double > &a,
                                    const CoordinateVector< double > &b,
                                    const CoordinateVector< double > &c,
                                    const CoordinateVector< double > &d) {
#ifdef HAVE_MULTIPRECISION
    const int_orient3d axp = get_mantissa(a.x());
    const int_orient3d ayp = get_mantissa(a.y());
    const int_orient3d azp = get_mantissa(a.z());

    const int_orient3d bxp = get_mantissa(b.x());
    const int_orient3d byp = get_mantissa(b.y());
    const int_orient3d bzp = get_mantissa(b.z());

    const int_orient3d cxp = get_mantissa(c.x());
    const int_orient3d cyp = get_mantissa(c.y());
    const int_orient3d czp = get_mantissa(c.z());

    const int_orient3d dxp = get_mantissa(d.x());
    const int_orient3d dyp = get_mantissa(d.y());
    const int_orient3d dzp = get_mantissa(d.z());

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
#else
    cmac_error("Cannot use exact geometric tests, as Boost Multiprecision was "
               "not found on the system!");
#endif
  }

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
   * This version uses double precision floating point arithmetics and
   * determines an error bound to check if round off error could change the
   * outcome of the test. If this is the case, an exact test is used, using the
   * given integer representations of the vertex coordinates.
   *
   * @param ar First vertex (real coordinates, in the range [1,2[).
   * @param br Second vertex (real coordinates, in the range [1,2[).
   * @param cr Third vertex (real coordinates, in the range [1,2[).
   * @param dr Fourth vertex (real coordinates, in the range [1,2[).
   * @return -1, 0, or 1, depending on the orientation of the tetrahedron.
   */
  inline static char orient3d_adaptive(const CoordinateVector<> &ar,
                                       const CoordinateVector<> &br,
                                       const CoordinateVector<> &cr,
                                       const CoordinateVector<> &dr) {

    const CoordinateVector<> ad = ar - dr;
    const CoordinateVector<> bd = br - dr;
    const CoordinateVector<> cd = cr - dr;

    const double bdxcdy = bd.x() * cd.y();
    const double cdxbdy = cd.x() * bd.y();

    const double cdxady = cd.x() * ad.y();
    const double adxcdy = ad.x() * cd.y();

    const double adxbdy = ad.x() * bd.y();
    const double bdxady = bd.x() * ad.y();

    const double errbound =
        1.e-10 * ((std::abs(bdxcdy) + std::abs(cdxbdy)) * std::abs(ad.z()) +
                  (std::abs(cdxady) + std::abs(adxcdy)) * std::abs(bd.z()) +
                  (std::abs(adxbdy) + std::abs(bdxady)) * std::abs(cd.z()));

    const double result = ad.z() * (bdxcdy - cdxbdy) +
                          bd.z() * (cdxady - adxcdy) +
                          cd.z() * (adxbdy - bdxady);

    if (result < -errbound) {
      return -1;
    } else if (result > errbound) {
      return 1;
    } else {
      return orient3d_exact(ar, br, cr, dr);
    }
  }

  /**
   * @brief Check if the fifth given point is inside the circumsphere of the
   * tetrahedron formed by the other four given points.
   *
   * It is assumed that the first four points are the vertices of a positively
   * oriented tetrahedron, as defined by a negative return value of orient3d().
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
  inline static char insphere_exact(const CoordinateVector<> &a,
                                    const CoordinateVector<> &b,
                                    const CoordinateVector<> &c,
                                    const CoordinateVector<> &d,
                                    const CoordinateVector<> &e) {
#ifdef HAVE_MULTIPRECISION
    const int_insphere axp = get_mantissa(a.x());
    const int_insphere ayp = get_mantissa(a.y());
    const int_insphere azp = get_mantissa(a.z());

    const int_insphere bxp = get_mantissa(b.x());
    const int_insphere byp = get_mantissa(b.y());
    const int_insphere bzp = get_mantissa(b.z());

    const int_insphere cxp = get_mantissa(c.x());
    const int_insphere cyp = get_mantissa(c.y());
    const int_insphere czp = get_mantissa(c.z());

    const int_insphere dxp = get_mantissa(d.x());
    const int_insphere dyp = get_mantissa(d.y());
    const int_insphere dzp = get_mantissa(d.z());

    const int_insphere exp = get_mantissa(e.x());
    const int_insphere eyp = get_mantissa(e.y());
    const int_insphere ezp = get_mantissa(e.z());

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
#else
    cmac_error("Cannot use exact geometric tests, as Boost Multiprecision was "
               "not found on the system!");
#endif
  }

  /**
   * @brief Check if the fifth given point is inside the circumsphere of the
   * tetrahedron formed by the other four given points.
   *
   * It is assumed that the first four points are the vertices of a positively
   * oriented tetrahedron, as defined by a negative return value of orient3d().
   *
   * If the fifth point is exactly on the circumsphere of the tetrahedron, this
   * functions returns 0.
   *
   * This version uses double precision floating point arithmetics and
   * determines an error bound to check if round off error could change the
   * outcome of the test. If this is the case, an exact test is used, using the
   * given integer representations of the vertex coordinates.
   *
   * @param ar First vertex (real coordinates, in the range [1,2[).
   * @param br Second vertex (real coordinates, in the range [1,2[).
   * @param cr Third vertex (real coordinates, in the range [1,2[).
   * @param dr Fourth vertex (real coordinates, in the range [1,2[).
   * @param er Fifth vertex (real coordinates, in the range [1,2[).
   * @return -1, 0, or 1, depending on the outcome of the geometric test.
   */
  inline static char insphere_adaptive(const CoordinateVector<> &ar,
                                       const CoordinateVector<> &br,
                                       const CoordinateVector<> &cr,
                                       const CoordinateVector<> &dr,
                                       const CoordinateVector<> &er) {

    const CoordinateVector<> ae = ar - er;
    const CoordinateVector<> be = br - er;
    const CoordinateVector<> ce = cr - er;
    const CoordinateVector<> de = dr - er;

    const double aexbey = ae.x() * be.y();
    const double bexaey = be.x() * ae.y();
    const double ab = aexbey - bexaey;
    const double bexcey = be.x() * ce.y();
    const double cexbey = ce.x() * be.y();
    const double bc = bexcey - cexbey;
    const double cexdey = ce.x() * de.y();
    const double dexcey = de.x() * ce.y();
    const double cd = cexdey - dexcey;
    const double dexaey = de.x() * ae.y();
    const double aexdey = ae.x() * de.y();
    const double da = dexaey - aexdey;
    const double aexcey = ae.x() * ce.y();
    const double cexaey = ce.x() * ae.y();
    const double ac = aexcey - cexaey;
    const double bexdey = be.x() * de.y();
    const double dexbey = de.x() * be.y();
    const double bd = bexdey - dexbey;

    const double abc = ae.z() * bc - be.z() * ac + ce.z() * ab;
    const double bcd = be.z() * cd - ce.z() * bd + de.z() * bc;
    const double cda = ce.z() * da + de.z() * ac + ae.z() * cd;
    const double dab = de.z() * ab + ae.z() * bd + be.z() * da;

    const double aenrm2 = ae.norm2();
    const double benrm2 = be.norm2();
    const double cenrm2 = ce.norm2();
    const double denrm2 = de.norm2();

    const double aezplus = std::abs(ae.z());
    const double bezplus = std::abs(be.z());
    const double cezplus = std::abs(ce.z());
    const double dezplus = std::abs(de.z());
    const double aexbeyplus = std::abs(aexbey);
    const double bexaeyplus = std::abs(bexaey);
    const double bexceyplus = std::abs(bexcey);
    const double cexbeyplus = std::abs(cexbey);
    const double cexdeyplus = std::abs(cexdey);
    const double dexceyplus = std::abs(dexcey);
    const double dexaeyplus = std::abs(dexaey);
    const double aexdeyplus = std::abs(aexdey);
    const double aexceyplus = std::abs(aexcey);
    const double cexaeyplus = std::abs(cexaey);
    const double bexdeyplus = std::abs(bexdey);
    const double dexbeyplus = std::abs(dexbey);
    const double errbound = 1.e-10 * (((cexdeyplus + dexceyplus) * bezplus +
                                       (dexbeyplus + bexdeyplus) * cezplus +
                                       (bexceyplus + cexbeyplus) * dezplus) *
                                          aenrm2 +
                                      ((dexaeyplus + aexdeyplus) * cezplus +
                                       (aexceyplus + cexaeyplus) * dezplus +
                                       (cexdeyplus + dexceyplus) * aezplus) *
                                          benrm2 +
                                      ((aexbeyplus + bexaeyplus) * dezplus +
                                       (bexdeyplus + dexbeyplus) * aezplus +
                                       (dexaeyplus + aexdeyplus) * bezplus) *
                                          cenrm2 +
                                      ((bexceyplus + cexbeyplus) * aezplus +
                                       (cexaeyplus + aexceyplus) * bezplus +
                                       (aexbeyplus + bexaeyplus) * cezplus) *
                                          denrm2);

    const double result =
        (denrm2 * abc - cenrm2 * dab) + (benrm2 * cda - aenrm2 * bcd);

    if (result < -errbound) {
      return -1;
    } else if (result > errbound) {
      return 1;
    } else {
      return insphere_exact(ar, br, cr, dr, er);
    }
  }
};

#endif // EXACTGEOMETRICTESTS_HPP
