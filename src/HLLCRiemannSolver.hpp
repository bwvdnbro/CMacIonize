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
 * @file HLLCRiemannSolver.hpp
 *
 * @brief HLLC Riemann solver.
 *
 * Based on the HLLC Riemann solver in Shadowfax (Vandenbroucke & De Rijcke,
 * 2016).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HLLCRIEMANNSOLVER_HPP
#define HLLCRIEMANNSOLVER_HPP

#include "RiemannSolver.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>

/**
 * @brief HLLC Riemann solver.
 */
class HLLCRiemannSolver : public RiemannSolver {
private:
  /*! @brief Polytropic index of the gas, @f$\gamma{}@f$. */
  const double _gamma;

  /*! @brief @f$\frac{2}{\gamma{}-1}@f$. */
  const double _tdgm1;

  /*! @brief @f$\frac{\gamma{}+1}{2\gamma{}}@f$. */
  const double _gp1d2g;

  /*! @brief @f$\frac{1}{\gamma{}-1}@f$. */
  const double _odgm1;

  /*! @brief @f$\frac{\gamma{}-1}{2}@f$. */
  const double _gm1d2;

  /*! @brief @f$\frac{\gamma{}-1}{\gamma{}+1}@f$. */
  const double _gm1dgp1;

  /*! @brief @f$\frac{2\gamma{}}{\gamma{}-1}@f$. */
  const double _tgdgm1;

  /*! @brief @f$\frac{2}{\gamma{}+1}@f$. */
  const double _tdgp1;

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
   * @return Flag indicating wether the right state (1) or a vacuum state (0)
   * was sampled.
   */
  inline int_fast32_t sample_right_vacuum(double rhoL, double uL, double PL,
                                          double aL, double &rhosol,
                                          double &usol, double &Psol) const {
    if (uL < aL) {
      /// vacuum regime
      // get the vacuum rarefaction wave speed
      const double SL = uL + _tdgm1 * aL;
      if (SL > 0.) {
        /// rarefaction wave regime
        // variable used twice below
        const double base = _tdgp1 + _gm1dgp1 * uL / aL;
        rhosol = rhoL * std::pow(base, _tdgm1);
        usol = _tdgp1 * (aL + _gm1d2 * uL);
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
   * @return Flag indicating wether the left state (-1) or a vacuum state (0)
   * was sampled.
   */
  inline int_fast32_t sample_left_vacuum(double rhoR, double uR, double PR,
                                         double aR, double &rhosol,
                                         double &usol, double &Psol) const {
    if (-aR < uR) {
      /// vacuum regime
      // get the vacuum rarefaction wave speed
      const double SR = uR - _tdgm1 * aR;
      if (SR < 0.) {
        /// rarefaction wave regime
        // variable used twice below
        const double base = _tdgp1 - _gm1dgp1 * uR / aR;
        rhosol = rhoR * std::pow(base, _tdgm1);
        usol = _tdgp1 * (-aR + _tdgm1 * uR);
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
   * @return Flag indicating wether the left state (-1), the right state (1), or
   * a vacuum state (0) was sampled.
   */
  inline int_fast32_t
  sample_vacuum_generation(double rhoL, double uL, double PL, double aL,
                           double rhoR, double uR, double PR, double aR,
                           double &rhosol, double &usol, double &Psol) const {

    // get the speeds of the left and right rarefaction waves
    const double SR = uR - _tdgm1 * aR;
    const double SL = uL + _tdgm1 * aL;
    if (SR > 0. && SL < 0.) {
      /// vacuum
      rhosol = 0.;
      usol = 0.;
      Psol = 0.;
      return 0;
    } else {
      if (SL < 0.) {
        /// right state
        if (-aR < uR) {
          /// right rarefaction wave regime
          // variable used twice below
          const double base = _tdgp1 - _gm1dgp1 * uR / aR;
          rhosol = rhoR * std::pow(base, _tdgm1);
          usol = _tdgp1 * (-aR + _tdgm1 * uR);
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
        if (aL > uL) {
          /// left rarefaction wave regime
          // variable used twice below
          const double base = _tdgp1 + _gm1dgp1 * uL / aL;
          rhosol = rhoL * std::pow(base, _tdgm1);
          usol = _tdgp1 * (aL + _tdgm1 * uL);
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
   * @param rhoL Left state density.
   * @param uL Left state velocity.
   * @param PL Left state pressure.
   * @param uLface Left state velocity deboosted to the interface frame.
   * @param vL Left state velocity along the interface normal.
   * @param aL Left state sound speed.
   * @param rhoR Right state density.
   * @param uR Right state velocity.
   * @param PR Right state pressure.
   * @param uRface Right state velocity deboosted to the interface frame.
   * @param vR Right state velocity along the interface normal.
   * @param aR Right state sound speed.
   * @param mflux Mass flux solution.
   * @param pflux Momentum flux solution.
   * @param Eflux Energy flux solution.
   * @param normal Surface normal of the interface.
   * @param vface Velocity of the interface, used to boost the fluxes.
   */
  inline void
  solve_vacuum_flux(const double rhoL, const CoordinateVector<> uL,
                    const double PL, const CoordinateVector<> uLface,
                    const double vL, const double aL, const double rhoR,
                    const CoordinateVector<> uR, const double PR,
                    const CoordinateVector<> uRface, const double vR,
                    const double aR, double &mflux, CoordinateVector<> &pflux,
                    double &Eflux, const CoordinateVector<> normal,
                    const CoordinateVector<> vface = 0.) const {

    // solve the Riemann problem
    double rhosol, vsol, Psol;
    int_fast32_t flag;
    if (rhoR == 0.) {
      /// vacuum right state
      flag = sample_right_vacuum(rhoL, vL, PL, aL, rhosol, vsol, Psol);
    } else if (rhoL == 0.) {
      /// vacuum left state
      flag = sample_left_vacuum(rhoR, vR, PR, aR, rhosol, vsol, Psol);
    } else {
      /// vacuum "middle" state
      flag = sample_vacuum_generation(rhoL, vL, PL, aL, rhoR, vR, PR, aR,
                                      rhosol, vsol, Psol);
    }

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
        rhoesol = 0.5 * rhosol * usol.norm2() + Psol * _odgm1;
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

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Polytropic index of the gas.
   */
  inline HLLCRiemannSolver(const double gamma)
      : _gamma(gamma), _tdgm1(2. / (gamma - 1.)),
        _gp1d2g(0.5 * (gamma + 1.) / gamma), _odgm1(1. / (gamma - 1.)),
        _gm1d2(0.5 * (gamma - 1.)), _gm1dgp1((gamma - 1.) / (gamma + 1.)),
        _tgdgm1(2. * gamma / (gamma - 1.)), _tdgp1(2. / (gamma + 1.)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~HLLCRiemannSolver() {}

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

    if (rhoL == 0. && rhoR == 0.) {
      // pure vacuum: all fluxes are 0
      return;
    }

    // deboost the velocities to the frame of the interface
    const CoordinateVector<> uLface = uL - vface;
    const CoordinateVector<> uRface = uR - vface;

    // get the velocities along the surface normal of the interface
    const double vL = CoordinateVector<>::dot_product(uLface, normal);
    const double vR = CoordinateVector<>::dot_product(uRface, normal);

    const double rhoLinv = 1. / rhoL;
    const double rhoRinv = 1. / rhoR;
    const double aL = std::sqrt(_gamma * PL * rhoLinv);
    const double aR = std::sqrt(_gamma * PR * rhoRinv);

    const double vdiff = vR - vL;
    const double abar = aL + aR;

    // Handle vacuum: vacuum does not require iteration and is always exact
    if (rhoL == 0. || rhoR == 0. || _tdgm1 * abar <= vdiff) {
      solve_vacuum_flux(rhoL, uL, PL, uLface, vL, aL, rhoR, uR, PR, uRface, vR,
                        aR, mflux, pflux, Eflux, normal, vface);
      return;
    }

    // STEP 1: pressure estimate
    const double rhobar = rhoL + rhoR;
    const double pPVRS = 0.5 * ((PL + PR) - 0.25 * vdiff * rhobar * abar);
    const double pstar = std::max(0., pPVRS);

    // STEP 2: wave speed estimates
    // all these speeds are along the interface normal, since uL and uR are
    double qL = 1.;
    if (pstar > PL) {
      qL = std::sqrt(1. + _gp1d2g * (pstar / PL - 1.));
    }
    double qR = 1.;
    if (pstar > PR) {
      qR = std::sqrt(1. + _gp1d2g * (pstar / PR - 1.));
    }

    // we only use the relative speeds below, so these expressions differ from
    // Toro (2009)
    const double SLmvL = -aL * qL;
    const double SRmvR = aR * qR;
    const double Sstar = (PR - PL + rhoL * vL * SLmvL - rhoR * vR * SRmvR) /
                         (rhoL * SLmvL - rhoR * SRmvR);

    if (Sstar >= 0.) {
      const double rhoLvL = rhoL * vL;
      const double vL2 = uLface.norm2();
      const double eL = PL * _odgm1 * rhoLinv + 0.5 * vL2;
      const double SL = SLmvL + vL;

      mflux = rhoLvL;
      pflux = rhoLvL * uLface + PL * normal;
      Eflux = rhoLvL * eL + PL * vL;

      if (SL < 0.) {
        const double starfac = SLmvL / (SL - Sstar) - 1.;
        const double SLrhoL = SL * rhoL;
        const double SstarmvL = Sstar - vL;
        const double SLrhoLstarfac = SLrhoL * starfac;
        const double SLrhoLSstarmvL = SLrhoL * SstarmvL;

        mflux += SLrhoLstarfac;
        pflux += SLrhoLstarfac * uLface + SLrhoLSstarmvL * normal;
        Eflux +=
            SLrhoLstarfac * eL + SLrhoLSstarmvL * (Sstar + PL / (rhoL * SLmvL));
      }
    } else {
      const double rhoRvR = rhoR * vR;
      const double vR2 = uRface.norm2();
      const double eR = PR * _odgm1 * rhoRinv + 0.5 * vR2;
      const double SR = SRmvR + vR;

      mflux = rhoRvR;
      pflux = rhoRvR * uRface + PR * normal;
      Eflux = rhoRvR * eR + PR * vR;

      if (SR > 0.) {
        const double starfac = SRmvR / (SR - Sstar) - 1.;
        const double SRrhoR = SR * rhoR;
        const double SstarmvR = Sstar - vR;
        const double SRrhoRstarfac = SRrhoR * starfac;
        const double SRrhoRSstarmvR = SRrhoR * SstarmvR;

        mflux += SRrhoRstarfac;
        pflux += SRrhoRstarfac * uRface + SRrhoRSstarmvR * normal;
        Eflux +=
            SRrhoRstarfac * eR + SRrhoRSstarmvR * (Sstar + PR / (rhoR * SRmvR));
      }
    }

    const double vface2 = vface.norm2();
    Eflux +=
        CoordinateVector<>::dot_product(vface, pflux) + 0.5 * vface2 * mflux;
    pflux += mflux * vface;
  }
};

#endif // HLLCRIEMANNSOLVER_HPP
