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
 * @file LambertW.hpp
 *
 * @brief Own implementation of the Lambert \f$W\f$ functions (branches
 * \f$W_0\f$ and \f$W_{-1}\f$).
 *
 * We only provide the real part of the two real functions Lambert \f$W_0\f$ and
 * Lambert \f$W_{-1}\f$ in their real range, \f$[-\frac{1}{e}, 0[\f$.
 *
 * This file was based on the Lambert W page on Wikipedia:
 * https://en.wikipedia.org/wiki/Lambert_W_function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LAMBERTW_HPP
#define LAMBERTW_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>

/**
 * @brief Own implementation of the Lambert \f$W\f$ functions (branches
 * \f$W_0\f$ and \f$W_{-1}\f$).
 */
class LambertW {
private:
  /**
   * @brief Get an initial guess for the value of the given branch of the
   * Lambert \f$W\f$ function for the given input value \f$r\f$.
   *
   * For the Lambert \f$W_0\f$ function, we use the fifth order Taylor expansion
   * as given on Wikipedia. For the Lambert \f$W_{-1}\f$ function, we make use
   * of the fact that the function is almost symmetric to Lambert \f$W_0\f$
   * w.r.t. \f$W = -1\f$ close to \f$r = -\frac{1}{e}\f$, so that \f$-2-W_0\f$
   * should be a close initial guess that makes sure our Newton method ends up
   * on the right branch.
   *
   * @param r Input value \f$r\f$.
   * @param branch Desired branch of the function (only 0 and -1 are supported,
   * default: 0).
   * @return Initial guess for the value of \f$W(r)\f$.
   */
  inline static double initial_guess(const double r, const int branch = 0) {
    // precompute r powers for maximal speed
    const double r2 = r * r;
    const double r3 = r2 * r;
    const double r4 = r2 * r2;
    const double r5 = r4 * r;
    const double w = r - r2 + 1.5 * r3 - (8. / 3.) * r4 + (125. / 24.) * r5;
    if (branch == 0) {
      return w;
    } else if (branch == -1) {
      // mirror the W0 initial guess w.r.t. W = -1
      return -2. - w;
    } else {
      // we only support the 0 and -1 branches
      std::cerr << "Unsupported Lambert W branch: " << branch << std::endl;
      std::abort();
    }
  }

  /**
   * @brief Perform a single Newton step to improve on a given guess value for
   * the Lambert \f$W\f$ function for the given input value \f$r\f$.
   *
   * Based on the formula given on Wikipedia.
   *
   * @param w Guess for the value of \f$W(r)\f$. This guess can be the result of
   * a call to initial_guess, or can be a result of a previous call to this
   * function.
   * @param r Input value \f$r\f$.
   * @return Better guess for \f$W(r)\f$.
   */
  inline static double newton_step(const double w, const double r) {
    // precompute terms for maximal speed
    const double expw = std::exp(w);
    const double wexpw = w * expw;
    return w - (wexpw - r) / (expw + wexpw);
  }

public:
  /**
   * @brief Get the value of the given branch of the Lambert \f$W\f$ function
   * for the given input value \f$r\f$.
   *
   * We use a Newton method to iteratively solve the defining equation
   * \f[
   *   r = W(r) {\rm{}e}^{W(r)}.
   * \f]
   *
   * @param r Input value \f$r\f$ (only \f$-\frac{1}{e} \leq{} r < 0\f$ are
   * supported).
   * @param branch Desired branch of the function (only 0 and -1 are supported,
   * default: 0).
   * @param tolerance Required relative accuracy for the result
   * (default: 1.e-10).
   * @return Value of \f$W(r)\f$ up to the required relative accuracy.
   */
  inline static double lambert_w(const double r, const int branch = 0,
                                 const double tolerance = 1.e-10) {
    if (r >= 0. || r < -1. / M_E) {
      std::cerr << "Input value for Lambert W outside supported range: " << r
                << " (supported range: " << (-1. / M_E) << " --> " << 0. << ")!"
                << std::endl;
      std::abort();
    }
    double w0 = initial_guess(r, branch);
    double w1 = newton_step(w0, r);
    while (std::abs(w0 - w1) > std::abs(w0 + w1) * tolerance) {
      w0 = w1;
      w1 = newton_step(w0, r);
    }
    return w1;
  }
};

#endif // LAMBERTW_HPP
