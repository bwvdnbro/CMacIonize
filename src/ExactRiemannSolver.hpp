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
 * @file ExactRiemannSolver.hpp
 *
 * @brief Exact Riemann solver.
 *
 * This Riemann solver is based on the ExactRiemannSolver class in the public
 * simulation code Shadowfax (Vandenbroucke & De Rijcke, 2016).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EXACTRIEMANNSOLVER_HPP
#define EXACTRIEMANNSOLVER_HPP

#include "Error.hpp"
#include "RiemannSolver.hpp"

#include <algorithm>
#include <cmath>

/**
 * @brief Exact Riemann solver.
 */
class ExactRiemannSolver : public RiemannSolver {
private:
  /*! @brief Adiabatic index @f$\gamma{}@f$. */
  const double _gamma;

  /*! @brief @f$\frac{\gamma+1}{2\gamma}@f$ */
  const double _gp1d2g;

  /*! @brief @f$\frac{\gamma-1}{2\gamma}@f$ */
  const double _gm1d2g;

  /*! @brief @f$\frac{\gamma-1}{\gamma+1}@f$ */
  const double _gm1dgp1;

  /*! @brief @f$\frac{2}{\gamma+1}@f$ */
  const double _tdgp1;

  /*! @brief @f$\frac{2}{\gamma-1}@f$ */
  const double _tdgm1;

  /*! @brief @f$\frac{\gamma-1}{2}@f$ */
  const double _gm1d2;

  /*! @brief @f$\frac{2\gamma}{\gamma-1}@f$ */
  const double _tgdgm1;

  /*! @brief @f$\frac{1}{\gamma}@f$ */
  const double _ginv;

  /*! @brief @f$\frac{1}{\gamma-1}@f$ */
  const double _gm1inv;

  /**
   * @brief Get the soundspeed corresponding to the given density and pressure.
   *
   * @param rho Density value.
   * @param P Pressure value.
   * @return Soundspeed.
   */
  inline double get_soundspeed(double rho, double P) const {
    const double result = std::sqrt(_gamma * P / rho);
    cmac_assert(result == result);
    return result;
  }

  /**
   * @brief Riemann fL or fR function.
   *
   * @param P Pressure of the left or right state.
   * @param A @f$\frac{2}{\gamma+1}\frac{1}{\rho}@f$.
   * @param B @f$\frac{\gamma-1}{\gamma+1}P@f$.
   * @param Pinv @f$\frac{1}{P}@f$.
   * @param afac @f$\frac{2}{\gamma-1}a@f$.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the fL or fR function.
   */
  inline double fb(double P, double A, double B, double Pinv, double afac,
                   double Pstar) const {
    if (Pstar > P) {
      const double result = (Pstar - P) * std::sqrt(A / (Pstar + B));
      cmac_assert(result == result);
      return result;
    } else {
      const double result = afac * (std::pow(Pstar * Pinv, _gm1d2g) - 1.);
      cmac_assert(result == result);
      return result;
    }
  }

  /**
   * @brief Riemann f function.
   *
   * @param PL Pressure of the left state.
   * @param AL @f$\frac{2}{\gamma+1}\frac{1}{\rho_L}@f$.
   * @param BL @f$\frac{\gamma-1}{\gamma+1}P_L@f$.
   * @param PLinv @f$\frac{1}{P_L}@f$.
   * @param aLfac @f$\frac{2}{\gamma-1}a_L@f$.
   * @param PR Pressure of the right state.
   * @param AR @f$\frac{2}{\gamma+1}\frac{1}{\rho_R}@f$.
   * @param BR @f$\frac{\gamma-1}{\gamma+1}P_R@f$.
   * @param PRinv @f$\frac{1}{P_R}@f$.
   * @param aRfac @f$\frac{2}{\gamma-1}a_R@f$.
   * @param udiff @f$u_R - u_L@f$.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the Riemann f function.
   */
  inline double f(double PL, double AL, double BL, double PLinv, double aLfac,
                  double PR, double AR, double BR, double PRinv, double aRfac,
                  double udiff, double Pstar) const {
    return fb(PL, AL, BL, PLinv, aLfac, Pstar) +
           fb(PR, AR, BR, PRinv, aRfac, Pstar) + udiff;
  }

  /**
   * @brief Derivative of the Riemann fL or fR function.
   *
   * @param P Pressure of the left or right state.
   * @param A @f$\frac{2}{\gamma+1}\frac{1}{\rho}@f$.
   * @param B @f$\frac{\gamma-1}{\gamma+1}P@f$.
   * @param Pinv @f$\frac{1}{P}@f$.
   * @param rhoainv @f$\frac{1}{\rho a}@f$.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the derivative of the Riemann fL or fR function.
   */
  inline double fprimeb(double P, double A, double B, double Pinv,
                        double rhoainv, double Pstar) const {
    if (Pstar > P) {
      const double C = 1. / (Pstar + B);
      const double result = (1. - 0.5 * (Pstar - P) * C) * std::sqrt(A * C);
      cmac_assert(result == result);
      return result;
    } else {
      const double result = std::pow(Pstar * Pinv, -_gp1d2g) * rhoainv;
      cmac_assert(result == result);
      return result;
    }
  }

  /**
   * @brief Derivative of the Riemann f function.
   *
   * @param PL Pressure of the left state.
   * @param AL @f$\frac{2}{\gamma+1}\frac{1}{\rho_L}@f$.
   * @param BL @f$\frac{\gamma-1}{\gamma+1}P_L@f$.
   * @param PLinv @f$\frac{1}{P_L}@f$.
   * @param rhoLaLinv @f$\frac{1}{\rho_L a_L}@f$.
   * @param PR Pressure of the right state.
   * @param AR @f$\frac{2}{\gamma+1}\frac{1}{\rho_R}@f$.
   * @param BR @f$\frac{\gamma-1}{\gamma+1}P_R@f$.
   * @param PRinv @f$\frac{1}{P_R}@f$.
   * @param rhoRaRinv @f$\frac{1}{\rho_R a_R}@f$.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the derivative of the Riemann f function.
   */
  inline double fprime(double PL, double AL, double BL, double PLinv,
                       double rhoLaLinv, double PR, double AR, double BR,
                       double PRinv, double rhoRaRinv, double Pstar) const {
    return fprimeb(PL, AL, BL, PLinv, rhoLaLinv, Pstar) +
           fprimeb(PR, AR, BR, PRinv, rhoRaRinv, Pstar);
  }

  /**
   * @brief Riemann gL or gR function.
   *
   * @param A @f$\frac{2}{\gamma+1}\frac{1}{\rho}@f$.
   * @param B @f$\frac{\gamma-1}{\gamma+1}P@f$.
   * @param Pstar (Temporary) pressure in the middle state.
   * @return Value of the gL or gR function.
   */
  inline double gb(double A, double B, double Pstar) const {
    const double result = std::sqrt(A / (Pstar + B));
    cmac_assert(result == result);
    return result;
  }

  /**
   * @brief Get an initial guess for the pressure in the middle state.
   *
   * @param PL Left state pressure.
   * @param aL Left state soundspeed.
   * @param AL @f$\frac{2}{\gamma+1}\frac{1}{\rho_L}@f$.
   * @param BL @f$\frac{\gamma-1}{\gamma+1}P_L@f$.
   * @param PR Right state pressure.
   * @param aR Right state soundspeed.
   * @param AR @f$\frac{2}{\gamma+1}\frac{1}{\rho_R}@f$.
   * @param BR @f$\frac{\gamma-1}{\gamma+1}P_R@f$.
   * @param udiff @f$u_R - u_L@f$.
   * @return Initial guess for the pressure in the middle state.
   */
  inline double guess_P(double PL, double aL, double AL, double BL, double PR,
                        double aR, double AR, double BR, double udiff) const {
    double Pguess;
    const double Pmin = std::min(PL, PR);
    const double Pmax = std::max(PL, PR);
    const double qmax = Pmax / Pmin;
    cmac_assert(qmax == qmax);
    const double PLpPR = PL + PR;
    const double smallP = 5.e-9 * PLpPR;
    double Ppv = 0.5 * PLpPR - 0.125 * udiff * PLpPR * (aL + aR);
    Ppv = std::max(smallP, Ppv);
    if (qmax <= 2. && Pmin <= Ppv && Ppv <= Pmax) {
      Pguess = Ppv;
    } else {
      if (Ppv < Pmin) {
        // two rarefactions
        Pguess =
            std::pow((aL + aR - _gm1d2 * udiff) / (aL * std::pow(PL, -_gm1d2g) +
                                                   aR * std::pow(PR, -_gm1d2g)),
                     _tgdgm1);
        cmac_assert(Pguess == Pguess);
      } else {
        // two shocks
        const double gL = gb(AL, BL, Ppv);
        const double gR = gb(AR, BR, Ppv);
        Pguess = (gL * PL + gR * PR - udiff) / (gL + gR);
        cmac_assert(Pguess == Pguess);
      }
    }
    // Toro: "Not that approximate solutions may predict, incorrectly, a
    // negative value for pressure (...). Thus in order to avoid negative guess
    // values we introduce the small positive constant _tolerance"
    // (tolerance is 1.e-8 in this case)
    Pguess = std::max(smallP, Pguess);
    return Pguess;
  }

  /**
   * @brief Find the pressure of the middle state by using Brent's method.
   *
   * @param PL Pressure of the left state.
   * @param AL @f$\frac{2}{\gamma+1}\frac{1}{\rho_L}@f$.
   * @param BL @f$\frac{\gamma-1}{\gamma+1}P_L@f$.
   * @param PLinv @f$\frac{1}{P_L}@f$.
   * @param aLfac @f$\frac{2}{\gamma-1}a_L@f$.
   * @param PR Pressure of the right state.
   * @param AR @f$\frac{2}{\gamma+1}\frac{1}{\rho_R}@f$.
   * @param BR @f$\frac{\gamma-1}{\gamma+1}P_R@f$.
   * @param PRinv @f$\frac{1}{P_R}@f$.
   * @param aRfac @f$\frac{2}{\gamma-1}a_R@f$.
   * @param udiff @f$u_R - u_L@f$.
   * @param Plow Lower bound guess for the pressure of the middle state.
   * @param Phigh Higher bound guess for the pressure of the middle state.
   * @param fPlow Value of the pressure function for the lower bound guess.
   * @param fPhigh Value of the pressure function for the upper bound guess.
   * @return Pressure of the middle state, with a 1.e-8 relative error
   * precision.
   */
  inline double solve_brent(double PL, double AL, double BL, double PLinv,
                            double aLfac, double PR, double AR, double BR,
                            double PRinv, double aRfac, double udiff,
                            double Plow, double Phigh, double fPlow,
                            double fPhigh) const {
    double a = Plow;
    double b = Phigh;
    double c = 0.;
    double d = 1e230;

    double fa = fPlow;
    double fb = fPhigh;
    double fc = 0.;

    double s = 0.;
    double fs = 0.;

    if (fa * fb > 0.) {
      cmac_error("Equal sign function values provided to solve_brent (%g %g)!",
                 fa, fb);
    }

    // if |f(a)| < |f(b)| then swap (a,b) end if
    if (std::abs(fa) < std::abs(fb)) {
      double tmp = a;
      a = b;
      b = tmp;
      tmp = fa;
      fa = fb;
      fb = tmp;
    }

    c = a;
    fc = fa;
    bool mflag = true;

    uint_fast32_t itcount = 0;
    // we bail after 1e4 iterations; whatever the value then should be good
    // enough
    while (itcount < 1e4 && !(fb == 0.) &&
           (std::abs(a - b) > 5.e-9 * (a + b))) {
      if ((fa != fc) && (fb != fc)) {
        // Inverse quadratic interpolation
        const double famfbinv = 1. / (fa - fb);
        const double famfcinv = 1. / (fa - fc);
        const double fbmfcinv = 1. / (fb - fc);
        s = a * fb * fc * famfbinv * famfcinv -
            b * fa * fc * famfbinv * fbmfcinv +
            c * fa * fb * famfcinv * fbmfcinv;
      } else {
        // Secant Rule
        s = b - fb * (b - a) / (fb - fa);
      }

      const double tmp2 = 0.25 * (3. * a + b);
      if (!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b))) ||
          (mflag && (std::abs(s - b) >= 0.5 * std::abs(b - c))) ||
          (!mflag && (std::abs(s - b) >= 0.5 * std::abs(c - d))) ||
          (mflag && (std::abs(b - c) < 5.e-9 * (b + c))) ||
          (!mflag && (std::abs(c - d) < 5.e-9 * (c + d)))) {
        s = 0.5 * (a + b);
        mflag = true;
      } else {
        mflag = false;
      }
      fs = f(PL, AL, BL, PLinv, aLfac, PR, AR, BR, PRinv, aRfac, udiff, s);
      d = c;
      c = b;
      fc = fb;
      if (fa * fs < 0.) {
        b = s;
        fb = fs;
      } else {
        a = s;
        fa = fs;
      }

      // if |f(a)| < |f(b)| then swap (a,b) end if
      if (std::abs(fa) < std::abs(fb)) {
        double tmp = a;
        a = b;
        b = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
      }

      ++itcount;
    }
    return b;
  }

  /**
   * @brief Sample the Riemann problem solution for a position in the right
   * shock wave regime.
   *
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param PRinv @f$\frac{1}{P_R}@f$.
   * @param ustar Velocity of the middle state.
   * @param Pstar Pressure of the middle state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   */
  inline void sample_right_shock_wave(double rhoR, double uR, double PR,
                                      double aR, double PRinv, double ustar,
                                      double Pstar, double &rhosol,
                                      double &usol, double &Psol,
                                      double dxdt = 0.) const {
    // variable used twice below
    const double PdPR = Pstar * PRinv;
    // get the shock speed
    const double SR = uR + aR * std::sqrt(_gp1d2g * PdPR + _gm1d2g);
    if (SR > dxdt) {
      /// middle state (shock) regime
      rhosol = rhoR * (PdPR + _gm1dgp1) / (_gm1dgp1 * PdPR + 1.);
      usol = ustar;
      Psol = Pstar;
    } else {
      /// right state regime
      rhosol = rhoR;
      usol = uR;
      Psol = PR;
    }
  }

  /**
   * @brief Sample the Riemann problem solution for a position in the right
   * rarefaction wave regime.
   *
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param PRinv @f$\frac{1}{P_R}@f$.
   * @param ustar Velocity of the middle state.
   * @param Pstar Pressure of the middle state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   */
  inline void sample_right_rarefaction_wave(double rhoR, double uR, double PR,
                                            double aR, double PRinv,
                                            double ustar, double Pstar,
                                            double &rhosol, double &usol,
                                            double &Psol,
                                            double dxdt = 0.) const {
    // get the velocity of the head of the rarefaction wave
    const double SHR = uR + aR;
    if (SHR > dxdt) {
      /// rarefaction wave regime
      // variable used twice below
      const double PdPR = Pstar * PRinv;
      // get the velocity of the tail of the rarefaction wave
      const double STR = ustar + aR * std::pow(PdPR, _gm1d2g);
      if (STR > dxdt) {
        /// middle state regime
        rhosol = rhoR * std::pow(PdPR, _ginv);
        usol = ustar;
        Psol = Pstar;
      } else {
        /// rarefaction fan regime
        // variable used twice below
        const double base = _tdgp1 - _gm1dgp1 * (uR - dxdt) / aR;
        rhosol = rhoR * std::pow(base, _tdgm1);
        usol = _tdgp1 * (-aR + _gm1d2 * uR + dxdt);
        Psol = PR * std::pow(base, _tgdgm1);
      }
    } else {
      /// right state regime
      rhosol = rhoR;
      usol = uR;
      Psol = PR;
    }
  }

  /**
   * @brief Sample the Riemann problem solution in the right state regime.
   *
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param PRinv @f$\frac{1}{P_R}@f$.
   * @param ustar Velocity of the middle state.
   * @param Pstar Pressure of the middle state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   */
  inline void sample_right_state(double rhoR, double uR, double PR, double aR,
                                 double PRinv, double ustar, double Pstar,
                                 double &rhosol, double &usol, double &Psol,
                                 double dxdt = 0.) const {
    if (Pstar > PR) {
      /// shock wave
      sample_right_shock_wave(rhoR, uR, PR, aR, PRinv, ustar, Pstar, rhosol,
                              usol, Psol, dxdt);
    } else {
      /// rarefaction wave
      sample_right_rarefaction_wave(rhoR, uR, PR, aR, PRinv, ustar, Pstar,
                                    rhosol, usol, Psol, dxdt);
    }
  }

  /**
   * @brief Sample the Riemann problem solution for a position in the left shock
   *  wave regime.
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param PLinv @f$\frac{1}{P_L}@f$.
   * @param ustar Velocity of the middle state.
   * @param Pstar Pressure of the middle state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   */
  inline void sample_left_shock_wave(double rhoL, double uL, double PL,
                                     double aL, double PLinv, double ustar,
                                     double Pstar, double &rhosol, double &usol,
                                     double &Psol, double dxdt = 0.) const {
    // variable used twice below
    const double PdPL = Pstar * PLinv;
    // get the shock speed
    const double SL = uL - aL * std::sqrt(_gp1d2g * PdPL + _gm1d2g);
    if (SL < dxdt) {
      /// middle state (shock) regime
      rhosol = rhoL * (PdPL + _gm1dgp1) / (_gm1dgp1 * PdPL + 1.);
      usol = ustar;
      Psol = Pstar;
    } else {
      /// left state regime
      rhosol = rhoL;
      usol = uL;
      Psol = PL;
    }
  }

  /**
   * @brief Sample the Riemann problem solution for a position in the left
   * rarefaction wave regime.
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param PLinv @f$\frac{1}{P_L}@f$.
   * @param ustar Velocity of the middle state.
   * @param Pstar Pressure of the middle state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   */
  inline void sample_left_rarefaction_wave(double rhoL, double uL, double PL,
                                           double aL, double PLinv,
                                           double ustar, double Pstar,
                                           double &rhosol, double &usol,
                                           double &Psol,
                                           double dxdt = 0.) const {
    // get the velocity of the head of the rarefaction wave
    const double SHL = uL - aL;
    if (SHL < dxdt) {
      /// rarefaction wave regime
      // variable used twice below
      const double PdPL = Pstar * PLinv;
      // get the velocity of the tail of the rarefaction wave
      const double STL = ustar - aL * std::pow(PdPL, _gm1d2g);
      if (STL > dxdt) {
        /// rarefaction fan regime
        // variable used twice below
        const double base = _tdgp1 + _gm1dgp1 * (uL - dxdt) / aL;
        rhosol = rhoL * std::pow(base, _tdgm1);
        usol = _tdgp1 * (aL + _gm1d2 * uL + dxdt);
        Psol = PL * std::pow(base, _tgdgm1);
      } else {
        /// middle state regime
        rhosol = rhoL * std::pow(PdPL, _ginv);
        usol = ustar;
        Psol = Pstar;
      }
    } else {
      /// left state regime
      rhosol = rhoL;
      usol = uL;
      Psol = PL;
    }
  }

  /**
   * @brief Sample the Riemann problem solution in the left state regime.
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param PLinv @f$\frac{1}{P_L}@f$.
   * @param ustar Velocity of the middle state.
   * @param Pstar Pressure of the middle state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   */
  inline void sample_left_state(double rhoL, double uL, double PL, double aL,
                                double PLinv, double ustar, double Pstar,
                                double &rhosol, double &usol, double &Psol,
                                double dxdt = 0.) const {
    if (Pstar > PL) {
      /// shock wave
      sample_left_shock_wave(rhoL, uL, PL, aL, PLinv, ustar, Pstar, rhosol,
                             usol, Psol, dxdt);
    } else {
      /// rarefaction wave
      sample_left_rarefaction_wave(rhoL, uL, PL, aL, PLinv, ustar, Pstar,
                                   rhosol, usol, Psol, dxdt);
    }
  }

  /**
   * @brief Sample the vacuum Riemann problem if the right state is a vacuum.
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   * @return Flag indicating wether the left state (-1), the right state (1), or
   * a vacuum state (0) was sampled.
   */
  inline int_fast32_t sample_right_vacuum(double rhoL, double uL, double PL,
                                          double aL, double &rhosol,
                                          double &usol, double &Psol,
                                          double dxdt = 0.) const {
    if (uL - aL < dxdt) {
      /// vacuum regime
      // get the vacuum rarefaction wave speed
      const double SL = uL + _tdgm1 * aL;
      if (SL > dxdt) {
        /// rarefaction wave regime
        // variable used twice below
        const double base = _tdgp1 + _gm1dgp1 * (uL - dxdt) / aL;
        rhosol = rhoL * std::pow(base, _tdgm1);
        usol = _tdgp1 * (aL + _gm1d2 * uL + dxdt);
        Psol = PL * std::pow(base, _tgdgm1);
        return -1;
      } else {
        /// vacuum
        rhosol = 0.;
        usol = 0.;
        Psol = 0.;
        return 0;
      }
    } else {
      /// left state regime
      rhosol = rhoL;
      usol = uL;
      Psol = PL;
      return -1;
    }
  }

  /**
   * @brief Sample the vacuum Riemann problem if the left state is a vacuum.
   *
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   * @return Flag indicating wether the left state (-1), the right state (1), or
   * a vacuum state (0) was sampled.
   */
  inline int_fast32_t sample_left_vacuum(double rhoR, double uR, double PR,
                                         double aR, double &rhosol,
                                         double &usol, double &Psol,
                                         double dxdt = 0.) const {
    if (dxdt < uR + aR) {
      /// vacuum regime
      // get the vacuum rarefaction wave speed
      const double SR = uR - _tdgm1 * aR;
      if (SR < dxdt) {
        /// rarefaction wave regime
        // variable used twice below
        const double base = _tdgp1 - _gm1dgp1 * (uR - dxdt) / aR;
        rhosol = rhoR * std::pow(base, _tdgm1);
        usol = _tdgp1 * (-aR + _tdgm1 * uR + dxdt);
        Psol = PR * std::pow(base, _tgdgm1);
        return 1;
      } else {
        /// vacuum
        rhosol = 0.;
        usol = 0.;
        Psol = 0.;
        return 0;
      }
    } else {
      /// right state regime
      rhosol = rhoR;
      usol = uR;
      Psol = PR;
      return 1;
    }
  }

  /**
   * @brief Sample the vacuum Riemann problem in the case vacuum is generated in
   * between the left and right state.
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   * @return Flag indicating wether the left state (-1), the right state (1), or
   * a vacuum state (0) was sampled.
   */
  inline int_fast32_t
  sample_vacuum_generation(double rhoL, double uL, double PL, double aL,
                           double rhoR, double uR, double PR, double aR,
                           double &rhosol, double &usol, double &Psol,
                           double dxdt) const {
    // get the speeds of the left and right rarefaction waves
    const double SR = uR - _tdgm1 * aR;
    const double SL = uL + _tdgm1 * aL;
    if (SR > dxdt && SL < dxdt) {
      /// vacuum
      rhosol = 0.;
      usol = 0.;
      Psol = 0.;
      return 0;
    } else {
      if (SL < dxdt) {
        /// right state
        if (dxdt < uR + aR) {
          /// right rarefaction wave regime
          // variable used twice below
          const double base = _tdgp1 - _gm1dgp1 * (uR - dxdt) / aR;
          rhosol = rhoR * std::pow(base, _tdgm1);
          usol = _tdgp1 * (-aR + _tdgm1 * uR + dxdt);
          Psol = PR * std::pow(base, _tgdgm1);
        } else {
          /// right state regime
          rhosol = rhoR;
          usol = uR;
          Psol = PR;
        }
        return 1;
      } else {
        /// left state
        if (dxdt > uL - aL) {
          /// left rarefaction wave regime
          // variable used twice below
          const double base = _tdgp1 + _gm1dgp1 * (uL - dxdt) / aL;
          rhosol = rhoL * std::pow(base, _tdgm1);
          usol = _tdgp1 * (aL + _tdgm1 * uL + dxdt);
          Psol = PL * std::pow(base, _tgdgm1);
        } else {
          /// left state regime
          rhosol = rhoL;
          usol = uL;
          Psol = PL;
        }
        return -1;
      }
    }
  }

  /**
   * @brief Vacuum Riemann solver.
   *
   * This solver is called when one or both states have a zero density, or when
   * the vacuum generation condition is satisfied (meaning vacuum is generated
   * in the middle state, although strictly speaking there is no "middle"
   * state if vacuum is involved).
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   * @return Flag indicating wether the left state (-1), the right state (1), or
   * a vacuum state (0) was sampled.
   */
  inline int_fast32_t solve_vacuum(double rhoL, double uL, double PL, double aL,
                                   double rhoR, double uR, double PR, double aR,
                                   double &rhosol, double &usol, double &Psol,
                                   double dxdt = 0.) const {
    // if both states are vacuum, the solution is also vacuum
    if (rhoL == 0. && rhoR == 0.) {
      rhosol = 0.;
      usol = 0.;
      Psol = 0.;
      return 0;
    }

    if (rhoR == 0.) {
      /// vacuum right state
      return sample_right_vacuum(rhoL, uL, PL, aL, rhosol, usol, Psol, dxdt);
    } else if (rhoL == 0.) {
      /// vacuum left state
      return sample_left_vacuum(rhoR, uR, PR, aR, rhosol, usol, Psol, dxdt);
    } else {
      /// vacuum "middle" state
      return sample_vacuum_generation(rhoL, uL, PL, aL, rhoR, uR, PR, aR,
                                      rhosol, usol, Psol, dxdt);
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * Note that we limit the adiabatic index to a value of 1.00000001. If an
   * adiabatic index of 1 is requested (isothermal gas), we handle the equation
   * of state elsewhere and assume this small value to make sure all derived
   * quantities are well defined.
   *
   * @param gamma Adiabatic index @f$\gamma{}@f$.
   */
  ExactRiemannSolver(const double gamma)
      : _gamma(std::max(gamma, 1.00000001)),
        _gp1d2g(0.5 * (_gamma + 1.) / _gamma),
        _gm1d2g(0.5 * (_gamma - 1.) / _gamma),
        _gm1dgp1((_gamma - 1) / (_gamma + 1.)), _tdgp1(2. / (_gamma + 1.)),
        _tdgm1(2. / (_gamma - 1.)), _gm1d2(0.5 * (_gamma - 1.)),
        _tgdgm1(2. * _gamma / (_gamma - 1.)), _ginv(1. / _gamma),
        _gm1inv(1. / (_gamma - 1.)) {

    if (gamma < 1.) {
      cmac_error("The adiabatic index needs to be 1 or larger!")
    }
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~ExactRiemannSolver() {}

  /**
   * @brief Solve the Riemann problem with the given left and right state.
   *
   * @param rhoL Left state density.
   * @param uL Left state velocity.
   * @param PL Left state pressure.
   * @param rhoR Right state density.
   * @param uR Right state velocity.
   * @param PR Right state pressure.
   * @param rhosol Density solution.
   * @param usol Velocity solution.
   * @param Psol Pressure solution.
   * @param dxdt Point in velocity space where we want to sample the solution.
   * @return Flag signaling whether the left state (-1), the right state (1), or
   * a vacuum state (0) was sampled.
   */
  inline int_fast32_t solve(const double rhoL, const double uL, const double PL,
                            const double rhoR, const double uR, const double PR,
                            double &rhosol, double &usol, double &Psol,
                            const double dxdt = 0.) const {

    // handle vacuum
    if (rhoL == 0. || rhoR == 0.) {
      double aL, aR;
      if (rhoL == 0.) {
        aL = 0.;
      } else {
        aL = get_soundspeed(rhoL, PL);
      }
      if (rhoR == 0.) {
        aR = 0.;
      } else {
        aR = get_soundspeed(rhoR, PR);
      }
      return solve_vacuum(rhoL, uL, PL, aL, rhoR, uR, PR, aR, rhosol, usol,
                          Psol, dxdt);
    }

    // get the soundspeeds
    const double aL = get_soundspeed(rhoL, PL);
    cmac_assert(aL == aL);
    const double aR = get_soundspeed(rhoR, PR);
    cmac_assert(aR == aR);
    const double aLfac = _tdgm1 * aL;
    const double aRfac = _tdgm1 * aR;
    const double udiff = uR - uL;

    // handle vacuum generation
    if (aLfac + aRfac <= udiff) {
      return solve_vacuum(rhoL, uL, PL, aL, rhoR, uR, PR, aR, rhosol, usol,
                          Psol, dxdt);
    }

    // precompute some variables
    const double AL = _tdgp1 / rhoL;
    cmac_assert(AL == AL);
    const double BL = _gm1dgp1 * PL;
    const double PLinv = 1. / PL;
    cmac_assert(PLinv == PLinv);
    const double rhoLaLinv = 1. / (rhoL * aL);
    cmac_assert(rhoLaLinv == rhoLaLinv);
    const double AR = _tdgp1 / rhoR;
    cmac_assert(AR == AR);
    const double BR = _gm1dgp1 * PR;
    const double PRinv = 1. / PR;
    cmac_assert(PRinv == PRinv);
    const double rhoRaRinv = 1. / (rhoR * aR);
    cmac_assert(rhoRaRinv == rhoRaRinv);

    // find the pressure and velocity in the middle state
    // since this is an exact Riemann solver, this is an iterative process,
    // whereby we basically find the root of a function (the Riemann f function
    // defined above)
    // we start by using a Newton-Raphson method, since we do not have an
    // interval in which the function changes sign
    // however, as soon as we have such an interval, we switch to a much more
    // robust root finding method (Brent's method). We do this because the
    // Newton-Raphson method in some cases can overshoot and return a negative
    // pressure, for which the Riemann f function is not defined. Brent's method
    // will never stroll outside of the initial interval in which the function
    // changes sign.
    double Pstar = 0.;
    double Pguess = guess_P(PL, aL, AL, BL, PR, aR, AR, BR, udiff);
    // we only store this variable to store the sign of the function for
    // pressure zero
    // we need to find a larger pressure for which this sign changes to have an
    // interval where we can use Brent's method
    double fPstar =
        f(PL, AL, BL, PLinv, aLfac, PR, AR, BR, PRinv, aRfac, udiff, Pstar);
    double fPguess =
        f(PL, AL, BL, PLinv, aLfac, PR, AR, BR, PRinv, aRfac, udiff, Pguess);
    if (fPstar * fPguess >= 0.) {
      // Newton-Raphson until convergence or until usable interval is
      // found to use Brent's method
      while (std::abs(Pstar - Pguess) > 5.e-9 * (Pstar + Pguess) &&
             fPguess < 0.) {
        Pstar = Pguess;
        fPstar = fPguess;
        Pguess -= fPguess / fprime(PL, AL, BL, PLinv, rhoLaLinv, PR, AR, BR,
                                   PRinv, rhoRaRinv, Pguess);
        fPguess = f(PL, AL, BL, PLinv, aLfac, PR, AR, BR, PRinv, aRfac, udiff,
                    Pguess);
      }
    }

    // As soon as there is a suitable interval: use Brent's method
    if (std::abs(Pstar - Pguess) > 5.e-9 * (Pstar + Pguess) && fPguess > 0.) {
      Pstar = solve_brent(PL, AL, BL, PLinv, aLfac, PR, AR, BR, PRinv, aRfac,
                          udiff, Pstar, Pguess, fPstar, fPguess);
    } else {
      Pstar = Pguess;
    }

    // the middle state velocity is fixed once the middle state pressure is
    // known
    const double ustar =
        0.5 * ((uL + uR) + (fb(PR, AR, BR, PRinv, aRfac, Pstar) -
                            fb(PL, AL, BL, PLinv, aLfac, Pstar)));

    // we now have solved the Riemann problem: we have the left, middle and
    // right state, and this completely fixes the solution
    // we just need to sample the solution for x/t = 0.
    if (ustar < dxdt) {
      // right state
      sample_right_state(rhoR, uR, PR, aR, PRinv, ustar, Pstar, rhosol, usol,
                         Psol, dxdt);
      return 1;
    } else {
      // left state
      sample_left_state(rhoL, uL, PL, aL, PLinv, ustar, Pstar, rhosol, usol,
                        Psol, dxdt);
      return -1;
    }
  }

  /**
   * @brief Solve the Riemann problem with the given left and right state and
   * get the resulting flux accross an interface.
   *
   * @param rhoL Left state density.
   * @param uL Left state velocity.
   * @param PL Left state pressure.
   * @param rhoR Right state density.
   * @param uR Right state velocity.
   * @param PR Right state pressure.
   * @param mflux Mass flux solution.
   * @param pflux Momentum flux solution.
   * @param Eflux Energy flux solution.
   * @param normal Surface normal of the interface.
   * @param vface Velocity of the interface, used to boost the fluxes.
   */
  virtual void solve_for_flux(const double rhoL, const CoordinateVector<> uL,
                              const double PL, const double rhoR,
                              const CoordinateVector<> uR, const double PR,
                              double &mflux, CoordinateVector<> &pflux,
                              double &Eflux, const CoordinateVector<> normal,
                              const CoordinateVector<> vface = 0.) const {

    // boost the velocities to the interface frame (and use new variables,
    // as we still want to use the old value of uL for other neighbours)
    const CoordinateVector<> uLface = uL - vface;
    const CoordinateVector<> uRface = uR - vface;

    // project the velocities onto the surface normal
    const double vL = CoordinateVector<>::dot_product(uLface, normal);
    const double vR = CoordinateVector<>::dot_product(uRface, normal);

    // solve the Riemann problem
    double rhosol, vsol, Psol;
    const int flag = solve(rhoL, vL, PL, rhoR, vR, PR, rhosol, vsol, Psol);

    // if the solution was vacuum, there is no flux
    if (flag != 0) {
      // deproject the velocity
      CoordinateVector<> usol;
      if (flag == -1) {
        vsol -= vL;
        usol = uLface + vsol * normal;
      } else {
        vsol -= vR;
        usol = uRface + vsol * normal;
      }

      // rho*e = rho*u + 0.5*rho*v^2 = P/(gamma-1.) + 0.5*rho*v^2
      double rhoesol;
      if (_gamma > 1.) {
        rhoesol = 0.5 * rhosol * usol.norm2() + Psol * _gm1inv;
      } else {
        // this flux will be ignored, but we make sure it has a sensible value
        rhoesol = 0.5 * rhosol * usol.norm2();
      }
      vsol = CoordinateVector<>::dot_product(usol, normal);

      // get the fluxes
      mflux = rhosol * vsol;
      pflux = rhosol * vsol * usol + Psol * normal;
      Eflux = (rhoesol + Psol) * vsol;

      // de-boost fluxes to fixed reference frame
      const double vface2 = vface.norm2();
      Eflux +=
          CoordinateVector<>::dot_product(vface, pflux) + 0.5 * vface2 * mflux;
      pflux += mflux * vface;
    }
  }
};

#endif // EXACTRIEMANNSOLVER_HPP
