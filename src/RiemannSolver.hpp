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
 * @file RiemannSolver.hpp
 *
 * @brief Exact Riemann solver.
 *
 * This Riemann solver is based on the ExactRiemannSolver class in the public
 * simulation code Shadowfax (Vandenbroucke & De Rijcke, 2016).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RIEMANNSOLVER_HPP
#define RIEMANNSOLVER_HPP

#include "Error.hpp"

#include <algorithm>
#include <cmath>

/**
 * @brief Exact Riemann solver.
 */
class RiemannSolver {
private:
  /*! @brief Adiabatic index \f$\gamma{}\f$. */
  double _gamma;

  /*! @brief \f$\frac{\gamma+1}{2\gamma}\f$ */
  double _gp1d2g;

  /*! @brief \f$\frac{\gamma-1}{2\gamma}\f$ */
  double _gm1d2g;

  /*! @brief \f$\frac{\gamma-1}{\gamma+1}\f$ */
  double _gm1dgp1;

  /*! @brief \f$\frac{2}{\gamma+1}\f$ */
  double _tdgp1;

  /*! @brief \f$\frac{2}{\gamma-1}\f$ */
  double _tdgm1;

  /*! @brief \f$\frac{\gamma-1}{2}\f$ */
  double _gm1d2;

  /*! @brief \f$\frac{2\gamma}{\gamma-1}\f$ */
  double _tgdgm1;

  /*! @brief \f$\frac{1}{\gamma}\f$ */
  double _ginv;

  /**
   * @brief Get the soundspeed corresponding to the given density and pressure.
   *
   * @param rho Density value.
   * @param P Pressure value.
   * @return Soundspeed.
   */
  inline double get_soundspeed(double rho, double P) {
    return std::sqrt(_gamma * P / rho);
  }

  /**
   * @brief Riemann fL or fR function.
   *
   * @param rho Density of the left or right state.
   * @param P Pressure of the left or right state.
   * @param a Soundspeed of the left or right state.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the fL or fR function.
   */
  inline double fb(double rho, double P, double a, double Pstar) {
    double fval = 0.;
    if (Pstar > P) {
      double A = _tdgp1 / rho;
      double B = _gm1dgp1 * P;
      fval = (Pstar - P) * std::sqrt(A / (Pstar + B));
    } else {
      fval = _tdgm1 * a * (std::pow(Pstar / P, _gm1d2g) - 1.);
    }
    return fval;
  }

  /**
   * @brief Riemann f function.
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the Riemann f function.
   */
  inline double f(double rhoL, double uL, double PL, double aL, double rhoR,
                  double uR, double PR, double aR, double Pstar) {
    return fb(rhoL, PL, aL, Pstar) + fb(rhoR, PR, aR, Pstar) + (uR - uL);
  }

  /**
   * @brief Derivative of the Riemann fL or fR function.
   *
   * @param rho Density of the left or right state.
   * @param P Pressure of the left or right state.
   * @param a Soundspeed of the left or right state.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the derivative of the Riemann fL or fR function.
   */
  inline double fprimeb(double rho, double P, double a, double Pstar) {
    double fval = 0.;
    if (Pstar > P) {
      double A = _tdgp1 / rho;
      double B = _gm1dgp1 * P;
      fval =
          (1. - 0.5 * (Pstar - P) / (B + Pstar)) * std::sqrt(A / (Pstar + B));
    } else {
      fval = 1. / (rho * a) * std::pow(Pstar / P, -_gp1d2g);
    }
    return fval;
  }

  /**
   * @brief Derivative of the Riemann f function.
   *
   * @param rhoL Density of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param rhoR Density of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param Pstar (Temporary) pressure of the middle state.
   * @return Value of the derivative of the Riemann f function.
   */
  inline double fprime(double rhoL, double PL, double aL, double rhoR,
                       double PR, double aR, double Pstar) {
    return fprimeb(rhoL, PL, aL, Pstar) + fprimeb(rhoR, PR, aR, Pstar);
  }

  /**
   * @brief Riemann gL or gR function.
   *
   * @param rho Density of the left or right state.
   * @param P Pressure of the left or right state.
   * @param Pstar (Temporary) pressure in the middle state.
   * @return Value of the gL or gR function.
   */
  inline double gb(double rho, double P, double Pstar) {
    double A = _tdgp1 / rho;
    double B = _gm1dgp1 * P;
    return std::sqrt(A / (Pstar + B));
  }

  /**
   * @brief Get an initial guess for the pressure in the middle state.
   *
   * @param rhoL Left state density.
   * @param uL Left state velocity.
   * @param PL Left state pressure.
   * @param aL Left state soundspeed.
   * @param rhoR Right state density.
   * @param uR Right state velocity.
   * @param PR Right state pressure.
   * @param aR Right state soundspeed.
   * @return Initial guess for the pressure in the middle state.
   */
  inline double guess_P(double rhoL, double uL, double PL, double aL,
                        double rhoR, double uR, double PR, double aR) {
    double Pguess;
    double Pmin = std::min(PL, PR);
    double Pmax = std::max(PL, PR);
    double qmax = Pmax / Pmin;
    double Ppv = 0.5 * (PL + PR) - 0.125 * (uR - uL) * (PL + PR) * (aL + aR);
    Ppv = std::max(5.e-9 * (PL + PR), Ppv);
    if (qmax <= 2. && Pmin <= Ppv && Ppv <= Pmax) {
      Pguess = Ppv;
    } else {
      if (Ppv < Pmin) {
        // two rarefactions
        Pguess = std::pow(
            (aL + aR - _gm1d2 * (uR - uL)) /
                (aL / std::pow(PL, _gm1d2g) + aR / std::pow(PR, _gm1d2g)),
            _tgdgm1);
      } else {
        // two shocks
        double gL = gb(rhoL, PL, Ppv);
        double gR = gb(rhoR, PR, Ppv);
        Pguess = (gL * PL + gR * PR - uR + uL) / (gL + gR);
      }
    }
    // Toro: "Not that approximate solutions may predict, incorrectly, a
    // negative value for pressure (...). Thus in order to avoid negative guess
    // values we introduce the small positive constant _tolerance"
    // (tolerance is 1.e-8 in this case)
    Pguess = std::max(5.e-9 * (PL + PR), Pguess);
    return Pguess;
  }

  /**
   * @brief Find the pressure of the middle state by using Brent's method.
   *
   * @param rhoL Density of the left state.
   * @param uL Velocity of the left state.
   * @param PL Pressure of the left state.
   * @param aL Soundspeed of the left state.
   * @param rhoR Density of the right state.
   * @param uR Velocity of the right state.
   * @param PR Pressure of the right state.
   * @param aR Soundspeed of the right state.
   * @param Plow Lower bound guess for the pressure of the middle state.
   * @param Phigh Higher bound guess for the pressure of the middle state.
   * @return Pressure of the middle state, with a 1.e-8 relative error
   * precision.
   */
  inline double solve_brent(double rhoL, double uL, double PL, double aL,
                            double rhoR, double uR, double PR, double aR,
                            double Plow, double Phigh) {
    double a = Plow;
    double b = Phigh;
    double c = 0.;
    double d = 1e230;

    double fa = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, a);
    double fb = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, b);
    double fc = 0.;

    double s = 0.;
    double fs = 0.;

    if (fa * fb > 0.) {
      cmac_error("Equal sign function values provided to solve_brent!");
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
    int i = 0;

    while (!(fb == 0.) && (std::abs(a - b) > 5.e-9 * (a + b))) {
      if ((fa != fc) && (fb != fc)) {
        // Inverse quadratic interpolation
        s = a * fb * fc / (fa - fb) / (fa - fc) +
            b * fa * fc / (fb - fa) / (fb - fc) +
            c * fa * fb / (fc - fa) / (fc - fb);
      } else {
        // Secant Rule
        s = b - fb * (b - a) / (fb - fa);
      }

      double tmp2 = 0.25 * (3. * a + b);
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
      fs = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, s);
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
      i++;
    }
    return b;
  }

  //  inline void sample_left_state(){
  //    if(Pstar > PL){
  //      // shock wave
  //    } else {
  //      // rarefaction wave
  //    }
  //  }

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Adiabatic index \f$\gamma{}\f$.
   */
  RiemannSolver(double gamma) : _gamma(gamma) {
    // related quantities:
    _gp1d2g = 0.5 * (_gamma + 1.) / _gamma; // gamma plus 1 divided by 2 gamma
    _gm1d2g = 0.5 * (_gamma - 1.) / _gamma; // gamma minus 1 divided by 2 gamma
    _gm1dgp1 =
        (_gamma - 1.) / (_gamma + 1.); // gamma minus 1 divided by gamma plus 1
    _tdgp1 = 2. / (_gamma + 1.);       // two divided by gamma plus 1
    _tdgm1 = 2. / (_gamma - 1.);       // two divided by gamma minus 1
    _gm1d2 = 0.5 * (_gamma - 1.);      // gamma minus 1 divided by 2
    _tgdgm1 =
        2. * _gamma / (_gamma - 1.); // two times gamma divided by gamma minus 1
    _ginv = 1. / _gamma;             // gamma inverse
  }

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
   */
  inline void solve(double rhoL, double uL, double PL, double rhoR, double uR,
                    double PR, double &rhosol, double &usol, double &Psol) {
    double aL = get_soundspeed(rhoL, PL);
    double aR = get_soundspeed(rhoR, PR);

    if (rhoL == 0. || rhoR == 0.) {
      cmac_error("Vacuum Riemann solver is not implemented yet!");
    }

    if (2. * aL / (_gamma - 1.) + 2. * aR / (_gamma - 1.) <= uR - uL) {
      cmac_error("Vacuum Riemann solver is not implemented yet!");
    } else {
      double Pstar = 0.;
      double Pguess = guess_P(rhoL, uL, PL, aL, rhoR, uR, PR, aR);
      double fPstar = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, Pstar);
      double fPguess = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, Pguess);
      if (fPstar * fPguess >= 0.) {
        // Newton-Raphson until convergence or until usable interval is
        // found to use Brent's method
        while (std::abs(Pstar - Pguess) > 5.e-9 * (Pstar + Pguess) &&
               fPguess < 0.) {
          Pstar = Pguess;
          Pguess =
              Pguess - fPguess / fprime(rhoL, PL, aL, rhoR, PR, aR, Pguess);
          fPguess = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, Pguess);
        }
      }

      // As soon as there is a suitable interval: use Brent's method
      if (std::abs(Pstar - Pguess) > 5.e-9 * (Pstar + Pguess) && fPguess > 0.) {
        Pstar = solve_brent(rhoL, uL, PL, aL, rhoR, uR, PR, aR, Pstar, Pguess);
      } else {
        Pstar = Pguess;
      }

      double ustar = 0.5 * (uL + uR) +
                     0.5 * (fb(rhoR, PR, aR, Pstar) - fb(rhoL, PL, aL, Pstar));
      (void)ustar;

      // sample the solution
      //      if(ustar < 0.){
      //        // left state
      //        sample_left_state();
      //      } else {
      //        // right state
      //        sample_right_state();
      //      }
      //        double vhalf;
      //        if(u < 0) {
      //            solution = WR;
      //            double pdpR = p / WR.p();
      //            if(p > WR.p()) {
      //                // shockwave
      //                double SR = vR + aR * sqrt(_gp1d2g * pdpR + _gm1d2g);
      //                if(SR > 0) {
      //                    solution.set_rho(WR.rho() * (pdpR + _gm1dgp1) /
      //                                     (_gm1dgp1 * pdpR + 1.));
      //                    solution.set_p(p);
      //                    vhalf = u - vR;
      //                } else {
      //                    // solution = WR
      //                    vhalf = 0.;
      //                }
      //            } else {
      //                // rarefaction wave
      //                double SHR = vR + aR;
      //                if(SHR > 0) {
      //                    double STR = u + aR * pow(pdpR, _gm1d2g);
      //                    if(STR <= 0) {
      //                        solution.set_rho(
      //                                WR.rho() *
      //                                pow(_tdgp1 - _gm1dgp1 / aR * vR,
      //                                _tdgm1));
      //                        vhalf = _tdgp1 * (-aR + _gm1d2 * vR) - vR;
      //                        solution.set_p(WR.p() * pow(_tdgp1 - _gm1dgp1 /
      //                        aR * vR,
      //                                                    _tgdgm1));
      //                    } else {
      //                        solution.set_rho(WR.rho() * pow(pdpR, _ginv));
      //                        solution.set_p(p);
      //                        vhalf = u - vR;
      //                    }
      //                } else {
      //                    // solution = WR
      //                    vhalf = 0.;
      //                }
      //            }
      //        } else {
      //            solution = WL;
      //            double pdpL = p / WL.p();
      //            if(p > WL.p()) {
      //                // shockwave
      //                double SL = vL - aL * sqrt(_gp1d2g * pdpL + _gm1d2g);
      //                if(SL < 0) {
      //                    solution.set_rho(WL.rho() * (pdpL + _gm1dgp1) /
      //                                     (_gm1dgp1 * pdpL + 1.));
      //                    solution.set_p(p);
      //                    vhalf = u - vL;
      //                } else {
      //                    // solution = WL
      //                    vhalf = 0.;
      //                }
      //            } else {
      //                // rarefaction wave
      //                double SHL = vL - aL;
      //                if(SHL < 0) {
      //                    double STL = u - aL * pow(pdpL, _gm1d2g);
      //                    if(STL > 0) {
      //                        solution.set_rho(
      //                                WL.rho() *
      //                                pow(_tdgp1 + _gm1dgp1 / aL * vL,
      //                                _tdgm1));
      //                        vhalf = _tdgp1 * (aL + _gm1d2 * vL) - vL;
      //                        solution.set_p(WL.p() * pow(_tdgp1 + _gm1dgp1 /
      //                        aL * vL,
      //                                                    _tgdgm1));
      //                    } else {
      //                        solution.set_rho(WL.rho() * pow(pdpL, _ginv));
      //                        vhalf = u - vL;
      //                        solution.set_p(p);
      //                    }
      //                } else {
      //                    // solution = WL
      //                    vhalf = 0.;
      //                }
      //            }
      //        }

      //        solution[1] += vhalf * n[0];
      //        solution[2] += vhalf * n[1];
      //#if ndim_ == 3
      //        solution[3] += vhalf * n[2];
      //#endif
    }
  }
};

#endif // RIEMANNSOLVER_HPP
