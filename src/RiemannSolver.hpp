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
  inline double guess_p(double rhoL, double uL, double PL, double aL,
                        double rhoR, double uR, double PR, double aR) {
    double pguess;
    double pmin = std::min(PL, PR);
    double pmax = std::max(PL, PR);
    double qmax = pmax / pmin;
    double ppv = 0.5 * (PL + PR) - 0.125 * (uR - uL) * (PL + PR) * (aL + aR);
    ppv = std::max(5.e-9 * (PL + PR), ppv);
    if (qmax <= 2. && pmin <= ppv && ppv <= pmax) {
      pguess = ppv;
    } else {
      if (ppv < pmin) {
        // two rarefactions
        pguess = std::pow(
            (aL + aR - _gm1d2 * (uR - uL)) /
                (aL / std::pow(PL, _gm1d2g) + aR / std::pow(PR, _gm1d2g)),
            _tgdgm1);
      } else {
        // two shocks
        double gL = gb(rhoL, PL, ppv);
        double gR = gb(rhoR, PR, ppv);
        pguess = (gL * PL + gR * PR - uR + uL) / (gL + gR);
      }
    }
    // Toro: "Not that approximate solutions may predict, incorrectly, a
    // negative value for pressure (...). Thus in order to avoid negative guess
    // values we introduce the small positive constant _tolerance"
    // (tolerance is 1.e-8 in this case)
    pguess = std::max(5.e-9 * (PL + PR), pguess);
    return pguess;
  }

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
   * @param rhostar Density solution.
   * @param ustar Velocity solution.
   * @param Pstar Pressure solution.
   */
  inline void solve(double rhoL, double uL, double PL, double rhoR, double uR,
                    double PR, double &rhostar, double &ustar, double &Pstar) {
    double aL = get_soundspeed(rhoL, PL);
    double aR = get_soundspeed(rhoR, PR);

    if (rhoL == 0. || rhoR == 0.) {
      cmac_error("Vacuum Riemann solver is not implemented yet!");
    }

    if (2. * aL / (_gamma - 1.) + 2. * aR / (_gamma - 1.) <= uR - uL) {
      cmac_error("Vacuum Riemann solver is not implemented yet!");
    } else {
      double p = 0.;
      double pguess = guess_p(rhoL, uL, PL, aL, rhoR, uR, PR, aR);
      double fp = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, p);
      double fpguess = f(rhoL, uL, PL, aL, rhoR, uR, PR, aR, pguess);
      (void)fp;
      (void)fpguess;
      //        if(fp < _cutoff) {
      //            while(fabs(p - pguess) > _tolerance * 0.5 * (p + pguess)) {
      //                p = pguess;
      //                pguess = pguess - fpguess / fprime(pguess, WL, WR, aL,
      //                aR);
      //                if(pguess < 0.) {
      //                    pguess = 0.;
      //                }
      //                fpguess = f(pguess, WL, WR, vL, vR, aL, aR);
      //            }
      //            p = pguess;
      //        } else {
      //            if(fp * fpguess >= 0.) {
      //                // Newton-Raphson until convergence or until usable
      //                interval is
      //                // found to use Brent's method
      //                while(fabs(p - pguess) > _tolerance * 0.5 * (p + pguess)
      //                &&
      //                      fpguess < 0.) {
      //                    p = pguess;
      //                    pguess = pguess - fpguess / fprime(pguess, WL, WR,
      //                    aL, aR);
      //                    fpguess = f(pguess, WL, WR, vL, vR, aL, aR);
      //                }
      //            }
      //            // As soon as there is a usable interval: use Brent's method
      //            if(fabs(p - pguess) > _tolerance * 0.5 * (p + pguess) &&
      //               fpguess > 0.) {
      //                p = 0.;
      //                _counterBrent++;
      //                p = BrentsMethodSolve(p, pguess, _tolerance, WL, WR, vL,
      //                vR, aL,
      //                                      aR);
      //            } else {
      //                p = pguess;
      //            }
      //        }

      //        double u = 0.5 * (vL + vR) + 0.5 * (fb(p, WR, aR) - fb(p, WL,
      //        aL));

      //        // set mach number
      //        if(p > WL.p()) {
      //            // left shock
      //            mach = std::max(mach, sqrt(_gp1d2g * p / WL.p() + _gm1d2g));
      //        }
      //        if(p > WR.p()) {
      //            // right shock
      //            mach = std::max(mach, sqrt(_gp1d2g * p / WR.p() + _gm1d2g));
      //        }

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
