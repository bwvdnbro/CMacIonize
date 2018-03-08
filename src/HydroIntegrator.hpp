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
 * @file HydroIntegrator.hpp
 *
 * @brief Class that performs the hydrodynamical integration.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROINTEGRATOR_HPP
#define HYDROINTEGRATOR_HPP

#include "BondiProfile.hpp"
#include "DensityGrid.hpp"
#include "DensityGridTraversalJobMarket.hpp"
#include "GradientCalculator.hpp"
#include "HydroBoundaryConditions.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "RiemannSolver.hpp"
#include "SimulationBox.hpp"
#include "Timer.hpp"

/*! @brief Uncomment this to switch off second order integration. */
//#define NO_SECOND_ORDER

#include <cfloat>

/*! @brief Stop the serial time timer and start the parallel time timer. */
#define hydro_start_parallel_timing_block()                                    \
  serial_timer.stop();                                                         \
  parallel_timer.start();

/*! @brief Stop the parallel time timer and start the serial time timer. */
#define hydro_stop_parallel_timing_block()                                     \
  parallel_timer.stop();                                                       \
  serial_timer.start();

/**
 * @brief Class that performs the hydrodynamical integration.
 */
class HydroIntegrator {
private:
  /*! @brief Adiabatic index of the gas. */
  const double _gamma;

  /*! @brief Adiabatic index minus one. */
  const double _gm1;

  /*! @brief Inverse of adiabatic index minus one. */
  const double _gm1_inv;

  /*! @brief Flag indicating whether we use radiative heating or not. */
  const bool _do_radiative_heating;

  /*! @brief Flag indicating whether we want radiative cooling or not. */
  const bool _do_radiative_cooling;

  /*! @brief Assumed temperature for neutral gas (in K). */
  const double _neutral_temperature;

  /*! @brief Assumed temperature for ionised gas (in K). */
  const double _ionised_temperature;

  /*! @brief Exact Riemann solver used to solve the Riemann problem. */
  const RiemannSolver _solver;

  /*! @brief Boundary conditions to apply to each boundary. */
  const HydroBoundaryConditionType _boundaries[6];

  /*! @brief Bondi profile used for Bondi boundary conditions. */
  const BondiProfile *_bondi_profile;

  /**
   * @brief Get the HydroBoundaryConditionType corresponding to the given type
   * string.
   *
   * @param type std::string representation of a boundary condition type.
   * @return Corresponding HydroBoundaryConditionType.
   */
  static HydroBoundaryConditionType get_boundary_type(std::string type) {
    if (type == "periodic") {
      return HYDRO_BOUNDARY_PERIODIC;
    } else if (type == "reflective") {
      return HYDRO_BOUNDARY_REFLECTIVE;
    } else if (type == "inflow") {
      return HYDRO_BOUNDARY_INFLOW;
    } else if (type == "bondi") {
      return HYDRO_BOUNDARY_BONDI;
    } else {
      cmac_error("Unknown boundary condition type: %s!", type.c_str());
      return HYDRO_BOUNDARY_INVALID;
    }
  }

  /**
   * @brief Functor that does the flux computation for a single cell.
   */
  class HydroFluxComputation {
  private:
    /*! @brief Reference to the underlying HydroIntegrator. */
    const HydroIntegrator &_hydro_integrator;

    /*! @brief Reference to the grid. */
    const DensityGrid &_grid;

    /*! @brief Iterator to the end of the grid. */
    const DensityGrid::iterator &_grid_end;

    /*! @brief Integration time step (in s). */
    const double _timestep;

    /**
     * @brief Per face slope limiter for a single quantity.
     *
     * Based on the slope limiter described in one of the appendices of Hopkins'
     * GIZMO paper.
     *
     * @param phimid0 Reconstructed value of the quantity at the interface.
     * @param phiL Value at the left of the interface.
     * @param phiR Value at the right of the interface.
     * @param dnrm_over_r Ratio of the distance between the left cell midpoint
     * to the face midpoint and the distances between left and right cell
     * midpoint.
     * @return Limited value of the quantity at the interface.
     */
    inline static double limit(const double phimid0, const double phiL,
                               const double phiR, const double dnrm_over_r) {

      const static double psi1 = 0.5;
      const static double psi2 = 0.25;

      const double delta1 = psi1 * std::abs(phiL - phiR);
      const double delta2 = psi2 * std::abs(phiL - phiR);

      const double phimin = std::min(phiL, phiR);
      const double phimax = std::max(phiL, phiR);

      const double phibar = phiL + dnrm_over_r * (phiR - phiL);

      // if sign(phimax+delta1) == sign(phimax)
      double phiplus;
      if ((phimax + delta1) * phimax > 0.) {
        phiplus = phimax + delta1;
      } else {
        const double absphimax = std::abs(phimax);
        phiplus = phimax * absphimax / (absphimax + delta1);
      }

      // if sign(phimin-delta1) == sign(phimin)
      double phiminus;
      if ((phimin - delta1) * phimin > 0.) {
        phiminus = phimin - delta1;
      } else {
        const double absphimin = std::abs(phimin);
        phiminus = phimin * absphimin / (absphimin + delta1);
      }

      double phimid;
      if (phiL == phiR) {
        phimid = phiL;
      } else {
        if (phiL < phiR) {
          phimid = std::max(phiminus, std::min(phibar + delta2, phimid0));
        } else {
          phimid = std::min(phiplus, std::max(phibar - delta2, phimid0));
        }
      }
      return phimid;
    }

  public:
    /**
     * @brief Constructor.
     *
     * @param hydro_integrator Reference to the underlying HydroIntegrator.
     * @param grid Reference to the DensityGrid.
     * @param grid_end Iterator to the end of the grid.
     * @param timestep Integration time step (in s).
     */
    inline HydroFluxComputation(const HydroIntegrator &hydro_integrator,
                                const DensityGrid &grid,
                                const DensityGrid::iterator &grid_end,
                                double timestep)
        : _hydro_integrator(hydro_integrator), _grid(grid), _grid_end(grid_end),
          _timestep(timestep) {}

    /**
     * @brief Do the flux computation for a single cell of the grid.
     *
     * @param cell DensityGrid::iterator pointing to a grid cell.
     */
    inline void operator()(DensityGrid::iterator &cell) {

      const CoordinateVector<> posL = cell.get_cell_midpoint();
      const double rhoL = cell.get_hydro_variables().get_primitives_density();
      const CoordinateVector<> uL =
          cell.get_hydro_variables().get_primitives_velocity();
      const double PL = cell.get_hydro_variables().get_primitives_pressure();

      const CoordinateVector<> gradrhoL =
          cell.get_hydro_variables().primitive_gradients(0);
      const CoordinateVector< CoordinateVector<> > graduL(
          cell.get_hydro_variables().primitive_gradients(1),
          cell.get_hydro_variables().primitive_gradients(2),
          cell.get_hydro_variables().primitive_gradients(3));
      const CoordinateVector<> gradPL =
          cell.get_hydro_variables().primitive_gradients(4);
      auto ngbs = cell.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        DensityGrid::iterator ngb = std::get< 0 >(*ngbit);
        // the midpoint is only used if we use a second order scheme
        const CoordinateVector<> midpoint = std::get< 1 >(*ngbit);
        const CoordinateVector<> normal = std::get< 2 >(*ngbit);
        const double surface_area = std::get< 3 >(*ngbit);
        const CoordinateVector<> posR = posL + std::get< 4 >(*ngbit);

        // get the right state
        double rhoR;
        CoordinateVector<> uR;
        double PR;
        CoordinateVector<> gradrhoR;
        CoordinateVector< CoordinateVector<> > graduR;
        CoordinateVector<> gradPR;
        CoordinateVector<> vframe;
        if (ngb != _grid_end) {
          rhoR = ngb.get_hydro_variables().get_primitives_density();
          uR = ngb.get_hydro_variables().get_primitives_velocity();
          PR = ngb.get_hydro_variables().get_primitives_pressure();
          gradrhoR = ngb.get_hydro_variables().primitive_gradients(0);
          graduR[0] = ngb.get_hydro_variables().primitive_gradients(1);
          graduR[1] = ngb.get_hydro_variables().primitive_gradients(2);
          graduR[2] = ngb.get_hydro_variables().primitive_gradients(3);
          gradPR = ngb.get_hydro_variables().primitive_gradients(4);
          vframe = _grid.get_interface_velocity(cell, ngb, midpoint);
        } else if (_hydro_integrator._boundaries[0] == HYDRO_BOUNDARY_BONDI) {
          const double r = posR.norm();
          double vR;
          _hydro_integrator._bondi_profile->get_hydrodynamic_variables(r, rhoR,
                                                                       vR, PR);
          uR[0] = vR * posR.x() / r;
          uR[1] = vR * posR.y() / r;
          uR[2] = vR * posR.z() / r;
          // we assume the gradients are just zero
        } else {
          // apply boundary conditions
          rhoR = rhoL;
          uR = uL;
          PR = PL;
          gradrhoR = gradrhoL;
          graduR = graduL;
          gradPR = gradPL;
          for (uint_fast8_t i = 0; i < 3; ++i) {
            if ((normal[i] < 0. &&
                 _hydro_integrator._boundaries[2 * i] ==
                     HYDRO_BOUNDARY_REFLECTIVE) ||
                (normal[i] > 0. &&
                 _hydro_integrator._boundaries[2 * i + 1] ==
                     HYDRO_BOUNDARY_REFLECTIVE)) {
              uR[i] = -uR[i];
              gradrhoR[i] = -gradrhoR[i];
              // we only invert the gradient components not orthogonal to the
              // face; the component orthogonal to the face has the same
              // gradient
              graduR[(i + 1) % 3][i] = -graduR[(i + 1) % 3][i];
              graduR[(i + 2) % 3][i] = -graduR[(i + 2) % 3][i];
              gradPR[i] = -gradPR[i];
            }
          }
        }

        // do the second order spatial gradient extrapolation
        const CoordinateVector<> dL = midpoint - posL;
        double rhoL_prime =
            rhoL + CoordinateVector<>::dot_product(gradrhoL, dL);
        CoordinateVector<> uL_prime(
            uL[0] + CoordinateVector<>::dot_product(graduL[0], dL),
            uL[1] + CoordinateVector<>::dot_product(graduL[1], dL),
            uL[2] + CoordinateVector<>::dot_product(graduL[2], dL));
        double PL_prime = PL + CoordinateVector<>::dot_product(gradPL, dL);

        const CoordinateVector<> dR = midpoint - posR;
        double rhoR_prime =
            rhoR + CoordinateVector<>::dot_product(gradrhoR, dR);
        CoordinateVector<> uR_prime(
            uR[0] + CoordinateVector<>::dot_product(graduR[0], dR),
            uR[1] + CoordinateVector<>::dot_product(graduR[1], dR),
            uR[2] + CoordinateVector<>::dot_product(graduR[2], dR));
        double PR_prime = PR + CoordinateVector<>::dot_product(gradPR, dR);

        // apply the per face slope limiter
        const double rinv = 1. / (posL - posR).norm();
        const double dL_over_r = dL.norm() * rinv;
        const double dR_over_r = dR.norm() * rinv;

        rhoL_prime = limit(rhoL_prime, rhoL, rhoR, dL_over_r);
        uL_prime[0] = limit(uL_prime[0], uL[0], uR[0], dL_over_r);
        uL_prime[1] = limit(uL_prime[1], uL[1], uR[1], dL_over_r);
        uL_prime[2] = limit(uL_prime[2], uL[2], uR[2], dL_over_r);
        PL_prime = limit(PL_prime, PL, PR, dL_over_r);

        rhoR_prime = limit(rhoR_prime, rhoR, rhoL, dR_over_r);
        uR_prime[0] = limit(uR_prime[0], uR[0], uL[0], dR_over_r);
        uR_prime[1] = limit(uR_prime[1], uR[1], uL[1], dR_over_r);
        uR_prime[2] = limit(uR_prime[2], uR[2], uL[2], dR_over_r);
        PR_prime = limit(PR_prime, PR, PL, dR_over_r);

        // boost the velocities to the interface frame (and use new variables,
        // as we still want to use the old value of uL for other neighbours)
        const CoordinateVector<> uLframe = uL_prime - vframe;
        const CoordinateVector<> uRframe = uR_prime - vframe;

        // project the velocities onto the surface normal
        const double vL = CoordinateVector<>::dot_product(uLframe, normal);
        const double vR = CoordinateVector<>::dot_product(uRframe, normal);

        // solve the Riemann problem
        double rhosol, vsol, Psol;
        const int flag = _hydro_integrator._solver.solve(
            rhoL_prime, vL, PL_prime, rhoR_prime, vR, PR_prime, rhosol, vsol,
            Psol);

        // if the solution was vacuum, there is no flux
        if (flag != 0) {
          // deproject the velocity
          CoordinateVector<> usol;
          if (flag == -1) {
            vsol -= vL;
            usol = uLframe + vsol * normal;
          } else {
            vsol -= vR;
            usol = uRframe + vsol * normal;
          }

          // rho*e = rho*u + 0.5*rho*v^2 = P/(gamma-1.) + 0.5*rho*v^2
          double rhoesol =
              0.5 * rhosol * usol.norm2() + Psol * _hydro_integrator._gm1_inv;
          vsol = CoordinateVector<>::dot_product(usol, normal);

          // get the fluxes
          const double mflux = rhosol * vsol * surface_area * _timestep;
          CoordinateVector<> pflux = rhosol * vsol * usol + Psol * normal;
          pflux *= surface_area * _timestep;
          double eflux = (rhoesol + Psol) * vsol * surface_area * _timestep;

          // de-boost fluxes to fixed reference frame
          const double vframe2 = vframe.norm2();
          eflux += CoordinateVector<>::dot_product(vframe, pflux) +
                   0.5 * vframe2 * mflux;
          pflux += mflux * vframe;

          cell.get_hydro_variables().delta_conserved(0) += mflux;
          cell.get_hydro_variables().delta_conserved(1) += pflux.x();
          cell.get_hydro_variables().delta_conserved(2) += pflux.y();
          cell.get_hydro_variables().delta_conserved(3) += pflux.z();
          cell.get_hydro_variables().delta_conserved(4) += eflux;
        }
      }
    }
  };

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Adiabatic index of the gas.
   * @param do_radiative_heating Flag indicating whether to use radiative
   * heating or not.
   * @param do_radiative_cooling Flag indicating whether to use radiative
   * cooling or not.
   * @param neutral_temperature Assumed temperature for neutral gas (in K).
   * @param ionised_temperature Assumed temperature for ionised gas (in K).
   * @param boundary_xlow Type of boundary for the lower x boundary.
   * @param boundary_xhigh Type of boundary for the upper x boundary.
   * @param boundary_ylow Type of boundary for the lower y boundary.
   * @param boundary_yhigh Type of boundary for the upper y boundary.
   * @param boundary_zlow Type of boundary for the lower z boundary.
   * @param boundary_zhigh Type of boundary for the upper z boundary.
   * @param box_periodicity Periodicity flags for the grid box (used to check
   * the validity of the boundary condition types).
   * @param bondi_profile BondiProfile object used for Bondi inflow boundary
   * conditions.
   */
  inline HydroIntegrator(double gamma, bool do_radiative_heating,
                         bool do_radiative_cooling,
                         double neutral_temperature = 100.,
                         double ionised_temperature = 1.e4,
                         std::string boundary_xlow = "reflective",
                         std::string boundary_xhigh = "reflective",
                         std::string boundary_ylow = "reflective",
                         std::string boundary_yhigh = "reflective",
                         std::string boundary_zlow = "reflective",
                         std::string boundary_zhigh = "reflective",
                         CoordinateVector< bool > box_periodicity =
                             CoordinateVector< bool >(false),
                         const BondiProfile *bondi_profile = nullptr)
      : _gamma(gamma), _gm1(_gamma - 1.), _gm1_inv(1. / _gm1),
        _do_radiative_heating(do_radiative_heating),
        _do_radiative_cooling(do_radiative_cooling),
        _neutral_temperature(neutral_temperature),
        _ionised_temperature(ionised_temperature), _solver(gamma),
        _boundaries{get_boundary_type(boundary_xlow),
                    get_boundary_type(boundary_xhigh),
                    get_boundary_type(boundary_ylow),
                    get_boundary_type(boundary_yhigh),
                    get_boundary_type(boundary_zlow),
                    get_boundary_type(boundary_zhigh)},
        _bondi_profile(bondi_profile) {

    if (_boundaries[0] == HYDRO_BOUNDARY_PERIODIC) {
      if (_boundaries[1] != HYDRO_BOUNDARY_PERIODIC) {
        cmac_error("Periodic boundaries in x only work if both x boundaries "
                   "are periodic!");
      }
      if (!box_periodicity[0]) {
        cmac_error("Periodic boundaries in x only work if the grid box is also "
                   "periodic in x!");
      }
    }
    if (_boundaries[2] == HYDRO_BOUNDARY_PERIODIC) {
      if (_boundaries[3] != HYDRO_BOUNDARY_PERIODIC) {
        cmac_error("Periodic boundaries in y only work if both y boundaries "
                   "are periodic!");
      }
      if (!box_periodicity[1]) {
        cmac_error("Periodic boundaries in y only work if the grid box is also "
                   "periodic in y!");
      }
    }
    if (_boundaries[4] == HYDRO_BOUNDARY_PERIODIC) {
      if (_boundaries[5] != HYDRO_BOUNDARY_PERIODIC) {
        cmac_error("Periodic boundaries in z only work if both z boundaries "
                   "are periodic!");
      }
      if (!box_periodicity[2]) {
        cmac_error("Periodic boundaries in z only work if the grid box is also "
                   "periodic in z!");
      }
    }

    if (_boundaries[0] == HYDRO_BOUNDARY_BONDI && _bondi_profile == nullptr) {
      cmac_error(
          "Bondi inflow boundaries only work if a Bondi profile is given.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - polytropic index: Polytropic index @f$\gamma{}@f$ of the gas (default:
   *    5. / 3.)
   *  - radiative heating: Is radiative heating enabled (default: true)?
   *  - radiative cooling: Is radiative cooling enabled (default: false)?
   *  - neutral temperature: Assumed temperature for neutral gas
   *    (default: 100. K)
   *  - ionised temperature: Assumed temperature for ionised gas
   *    (default: 1.e4 K)
   *  - boundary x low: Boundary condition type for the lower x boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary x high: Boundary condition type for the upper x boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary y low: Boundary condition type for the lower y boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary y high: Boundary condition type for the upper y boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary z low: Boundary condition type for the lower z boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary z high: Boundary condition type for the upper z boundary
   *    (periodic/reflective/inflow, default: reflective)
   *
   * @param simulation_box SimulationBox.
   * @param params ParameterFile to read from.
   */
  inline HydroIntegrator(const SimulationBox &simulation_box,
                         ParameterFile &params)
      : HydroIntegrator(
            params.get_value< double >("HydroIntegrator:polytropic index",
                                       5. / 3.),
            params.get_value< bool >("HydroIntegrator:radiative heating", true),
            params.get_value< bool >("HydroIntegrator:radiative cooling",
                                     false),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "HydroIntegrator:neutral temperature", "100. K"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "HydroIntegrator:ionised temperature", "1.e4 K"),
            params.get_value< std::string >("HydroIntegrator:boundary x low",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary x high",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary y low",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary y high",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary z low",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary z high",
                                            "reflective"),
            simulation_box.get_periodicity(), new BondiProfile(params)) {}

  /**
   * @brief Destructor.
   *
   * Clean up the Bondi profile.
   */
  inline ~HydroIntegrator() {
    if (_bondi_profile != nullptr) {
      delete _bondi_profile;
    }
  }

  /**
   * @brief Initialize the hydro variables for the given DensityGrid.
   *
   * @param grid DensityGrid to operate on.
   */
  inline void initialize_hydro_variables(DensityGrid &grid) const {

    const double hydrogen_mass =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    const double boltzmann_k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double volume = it.get_volume();
      const double number_density =
          it.get_ionization_variables().get_number_density();
      const double temperature =
          it.get_ionization_variables().get_temperature();

      const double density = number_density * hydrogen_mass;
      const CoordinateVector<> velocity =
          it.get_hydro_variables().get_primitives_velocity();
      // we assume a completely neutral or completely ionized gas
      double pressure = density * boltzmann_k * temperature / hydrogen_mass;
      if (temperature >= _ionised_temperature) {
        // ionized gas has a lower mean molecular mass
        pressure *= 2.;
      }

      // set the density and pressure (the velocity has been set by
      // DensityGrid::initialize)
      it.get_hydro_variables().set_primitives_density(density);
      it.get_hydro_variables().set_primitives_pressure(pressure);

      const double mass = density * volume;
      const CoordinateVector<> momentum = mass * velocity;
      const double ekin = CoordinateVector<>::dot_product(velocity, momentum);
      // E = V*(rho*u + 0.5*rho*v^2) = V*(P/(gamma-1) + 0.5*rho*v^2)
      const double total_energy = volume * pressure / _gm1 + 0.5 * ekin;

      // set conserved variables
      it.get_hydro_variables().set_conserved_mass(mass);
      it.get_hydro_variables().set_conserved_momentum(momentum);
      it.get_hydro_variables().set_conserved_total_energy(total_energy);
    }

    grid.set_grid_velocity(_gamma);
  }

  /**
   * @brief Get the maximal system time step that will lead to a stable
   * integration.
   *
   * @param grid DensityGrid on which to operate.
   * @return Maximal system time step that yields a stable integration (in s).
   */
  inline double get_maximal_timestep(DensityGrid &grid) const {

    double dtmin = DBL_MAX;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double rho = it.get_hydro_variables().get_primitives_density();
      const double P = it.get_hydro_variables().get_primitives_pressure();
      const double cs = std::sqrt(_gamma * P / rho);
      const double v =
          it.get_hydro_variables().get_primitives_velocity().norm();
      const double V = it.get_volume();
      const double R = std::cbrt(0.75 * V / M_PI);
      const double dt = 0.2 * R / (cs + v);
      dtmin = std::min(dt, dtmin);
    }
    return dtmin;
  }

  /**
   * @brief Do a single hydrodynamical time step.
   *
   * @param grid DensityGrid on which to operate.
   * @param timestep Time step over which to evolve the system (in s).
   * @param serial_timer Timer that times the time spent in serial parts of the
   * algorithm.
   * @param parallel_timer Timer that times the time spent in parallel parts of
   * the algorithm.
   */
  inline void do_hydro_step(DensityGrid &grid, double timestep,
                            Timer &serial_timer, Timer &parallel_timer) const {

    const DensityGrid::iterator grid_end = grid.end();
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());

#ifndef NO_SECOND_ORDER
    // if second order scheme: compute gradients for primitive variables
    GradientCalculator::GradientComputation gradient_computation(_boundaries,
                                                                 grid_end);
    WorkDistributor<
        DensityGridTraversalJobMarket<
            GradientCalculator::GradientComputation >,
        DensityGridTraversalJob< GradientCalculator::GradientComputation > >
        gradient_workers;
    DensityGridTraversalJobMarket< GradientCalculator::GradientComputation >
        gradient_jobs(grid, gradient_computation, block);
    hydro_start_parallel_timing_block();
    gradient_workers.do_in_parallel(gradient_jobs);
    hydro_stop_parallel_timing_block();
#endif

    // do the second order time prediction step
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double halfdt = 0.5 * timestep;

      // get primitive variables
      const double rho = it.get_hydro_variables().get_primitives_density();
      const CoordinateVector<> u =
          it.get_hydro_variables().get_primitives_velocity();
      const double P = it.get_hydro_variables().get_primitives_pressure();

      // get primitive gradients
      const CoordinateVector<> drho =
          it.get_hydro_variables().primitive_gradients(0);
      const CoordinateVector<> dux =
          it.get_hydro_variables().primitive_gradients(1);
      const CoordinateVector<> duy =
          it.get_hydro_variables().primitive_gradients(2);
      const CoordinateVector<> duz =
          it.get_hydro_variables().primitive_gradients(3);
      const CoordinateVector<> dP =
          it.get_hydro_variables().primitive_gradients(4);

      // compute updated variables
      const double divv = dux.x() + duy.y() + duz.z();
      const double rho_inv = 1. / rho;
      const double rho_new =
          rho -
          halfdt * (rho * divv + CoordinateVector<>::dot_product(u, drho));
      double ux_new = u.x() - halfdt * (u.x() * divv + rho_inv * dP.x());
      double uy_new = u.y() - halfdt * (u.y() * divv + rho_inv * dP.y());
      double uz_new = u.z() - halfdt * (u.z() * divv + rho_inv * dP.z());
      const double P_new =
          P -
          halfdt * (_gamma * P * divv + CoordinateVector<>::dot_product(u, dP));

      // add gravitational contribution
      const CoordinateVector<> a =
          it.get_hydro_variables().get_gravitational_acceleration();
      ux_new += halfdt * a.x();
      uy_new += halfdt * a.y();
      uz_new += halfdt * a.z();

      // update variables
      it.get_hydro_variables().primitives(0) = rho_new;
      it.get_hydro_variables().primitives(1) = ux_new;
      it.get_hydro_variables().primitives(2) = uy_new;
      it.get_hydro_variables().primitives(3) = uz_new;
      it.get_hydro_variables().primitives(4) = P_new;
    }

    // do the flux computation (in parallel)
    HydroFluxComputation hydro_flux_computation(*this, grid, grid_end,
                                                timestep);

    WorkDistributor< DensityGridTraversalJobMarket< HydroFluxComputation >,
                     DensityGridTraversalJob< HydroFluxComputation > >
        workers;
    DensityGridTraversalJobMarket< HydroFluxComputation > jobs(
        grid, hydro_flux_computation, block);
    hydro_start_parallel_timing_block();
    workers.do_in_parallel(jobs);
    hydro_stop_parallel_timing_block();

    // do radiation (if enabled)
    if (_do_radiative_heating || _do_radiative_cooling) {
      const double boltzmann_k =
          PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
      const double mH = PhysicalConstants::get_physical_constant(
          PHYSICALCONSTANT_PROTON_MASS);

      for (auto it = grid.begin(); it != grid.end(); ++it) {
        const IonizationVariables &ionization_variables =
            it.get_ionization_variables();

        const double xH = ionization_variables.get_ionic_fraction(ION_H_n);
        const double mpart = xH * mH + 0.5 * (1. - xH) * mH;
        const double Tgas =
            _ionised_temperature * (1. - xH) + _neutral_temperature * xH;
        const double ugas = boltzmann_k * Tgas / _gm1 / mpart;
        const double uold = it.get_hydro_variables().get_primitives_pressure() /
                            _gm1 /
                            it.get_hydro_variables().get_primitives_density();
        const double du = ugas - uold;
        const double dE = it.get_hydro_variables().get_conserved_mass() * du;
        if (_do_radiative_heating && dE > 0.) {
          it.get_hydro_variables().delta_conserved(4) -= dE;
        }
        if (_do_radiative_cooling && dE < 0.) {
          it.get_hydro_variables().delta_conserved(4) -= dE;
        }
      }
    }

    // update conserved variables
    for (auto it = grid.begin(); it != grid.end(); ++it) {

      it.get_hydro_variables().conserved(0) -=
          it.get_hydro_variables().delta_conserved(0);
      it.get_hydro_variables().conserved(1) -=
          it.get_hydro_variables().delta_conserved(1);
      it.get_hydro_variables().conserved(2) -=
          it.get_hydro_variables().delta_conserved(2);
      it.get_hydro_variables().conserved(3) -=
          it.get_hydro_variables().delta_conserved(3);
      it.get_hydro_variables().conserved(4) -=
          it.get_hydro_variables().delta_conserved(4);

      // add gravity
      const CoordinateVector<> a =
          it.get_hydro_variables().get_gravitational_acceleration();
      const double m = it.get_hydro_variables().get_conserved_mass();
      const CoordinateVector<> p =
          it.get_hydro_variables().get_conserved_momentum();
      it.get_hydro_variables().conserved(1) += timestep * m * a.x();
      it.get_hydro_variables().conserved(2) += timestep * m * a.y();
      it.get_hydro_variables().conserved(3) += timestep * m * a.z();
      it.get_hydro_variables().conserved(4) +=
          timestep * CoordinateVector<>::dot_product(p, a);

      // reset time differences
      it.get_hydro_variables().delta_conserved(0) = 0.;
      it.get_hydro_variables().delta_conserved(1) = 0.;
      it.get_hydro_variables().delta_conserved(2) = 0.;
      it.get_hydro_variables().delta_conserved(3) = 0.;
      it.get_hydro_variables().delta_conserved(4) = 0.;
    }

    grid.evolve(timestep);

    const double hydrogen_mass =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    const double boltzmann_k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    // convert conserved variables to primitive variables
    // also set the number density and temperature to the correct value
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double volume = it.get_volume();
      const double mass = it.get_hydro_variables().get_conserved_mass();
      const CoordinateVector<> momentum =
          it.get_hydro_variables().get_conserved_momentum();
      const double total_energy =
          it.get_hydro_variables().get_conserved_total_energy();

      double density, pressure;
      CoordinateVector<> velocity;
      if (mass <= 0.) {
        if (mass < 0.) {
          cmac_error("Negative mass for cell!");
        }
        // vacuum
        density = 0.;
        velocity = CoordinateVector<>(0.);
        pressure = 0.;
      } else {
        density = mass / volume;
        velocity = momentum / mass;
        // E = V*(rho*u + 0.5*rho*v^2) = (V*P/(gamma-1) + 0.5*m*v^2)
        // P = (E - 0.5*m*v^2)*(gamma-1)/V
        pressure = _gm1 *
                   (total_energy -
                    0.5 * CoordinateVector<>::dot_product(velocity, momentum)) /
                   volume;
      }

      cmac_assert(density >= 0.);
      cmac_assert(pressure >= 0.);

      it.get_hydro_variables().set_primitives_density(density);
      it.get_hydro_variables().set_primitives_velocity(velocity);
      it.get_hydro_variables().set_primitives_pressure(pressure);

      IonizationVariables &ionization_variables = it.get_ionization_variables();

      ionization_variables.set_number_density(density / hydrogen_mass);
      const double mean_molecular_mass =
          ionization_variables.get_ionic_fraction(ION_H_n) * hydrogen_mass +
          0.5 * (1. - ionization_variables.get_ionic_fraction(ION_H_n)) *
              hydrogen_mass;
      ionization_variables.set_temperature(mean_molecular_mass * pressure /
                                           boltzmann_k / density);

      cmac_assert(ionization_variables.get_number_density() >= 0.);
      cmac_assert(ionization_variables.get_temperature() >= 0.);
    }

    grid.set_grid_velocity(_gamma);
  }
};

#endif // HYDROINTEGRATOR_HPP
