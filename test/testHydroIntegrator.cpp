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
 * @file testHydroIntegrator.cpp
 *
 * @brief Unit test for the HydroIntegrator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "DensityFunction.hpp"
#include "ExactRiemannSolver.hpp"
#include "HydroIntegrator.hpp"
#include "RandomGenerator.hpp"
#include "VoronoiDensityGrid.hpp"
#include "VoronoiGeneratorDistribution.hpp"
#include <fstream>

/**
 * @brief DensityFunction implementation that sets up a basic Sod shock problem.
 */
class SodShockDensityFunction : public DensityFunction {
public:
  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {
    const CoordinateVector<> position = cell.get_cell_midpoint();
    const double hydrogen_mass =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    const double boltzmann_k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    const double density_unit = 1. / hydrogen_mass;
    const double temperature_unit = hydrogen_mass / boltzmann_k;
    DensityValues values;
    if (position.x() < 0.5) {
      values.set_number_density(density_unit);
      values.set_temperature(temperature_unit);
    } else {
      values.set_number_density(0.125 * density_unit);
      values.set_temperature(0.8 * temperature_unit);
    }
    return values;
  }
};

/**
 * @brief 1D Voronoi grid generator distribution.
 */
class OneDVoronoiGeneratorDistribution : public VoronoiGeneratorDistribution {
private:
  /*! @brief Index of the last returned generator position. */
  generatornumber_t _last_index;

public:
  /**
   * @brief Get the number of positions returned by this distribution.
   *
   * @return 100.
   */
  virtual generatornumber_t get_number_of_positions() const { return 100; }

  /**
   * @brief Get a generator position.
   *
   * @return Generator position (in m).
   */
  virtual CoordinateVector<> get_position() {
    CoordinateVector<> pos((_last_index + 0.5) * 0.01, 0.5, 0.5);
    ++_last_index;
    cmac_assert(_last_index <= 100);
    return pos;
  }
};

/**
 * @brief Unit test for the HydroIntegrator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// Flux symmetry check
  {
    const uint_fast32_t nrho = 10;
    const uint_fast32_t nu = 10;
    const uint_fast32_t nP = 10;

    const double rhoL = 1.;
    const CoordinateVector<> uL(1.);
    const double PL = 1.;

    const CoordinateVector<> posL(0.1, 0.1, 0.1);
    const CoordinateVector<> posR(0.2, 0.2, 0.2);
    const CoordinateVector<> midpoint(0.15, 0.15, 0.15);

    const CoordinateVector<> dL = midpoint - posL;
    const CoordinateVector<> dR = midpoint - posR;
    const double rinv = 1. / ((posL - posR).norm());
    const double dL_over_r = dL.norm() * rinv;
    const double dR_over_r = dR.norm() * rinv;

    CoordinateVector<> normal(1., 1., 1.);
    normal /= normal.norm();
    const CoordinateVector<> vframe;
    const double surface_area = 0.1;
    const double timestep = 0.1;

    const HLLCRiemannSolver solver(5. / 3.);

    RandomGenerator random_generator(42);

    double mfluxL, mfluxR, EfluxL, EfluxR;
    CoordinateVector<> pfluxL, pfluxR;
    for (uint_fast32_t irho = 0; irho < nrho; ++irho) {
      for (uint_fast32_t iu = 0; iu < nu; ++iu) {
        for (uint_fast32_t iP = 0; iP < nP; ++iP) {
          const double rhoR = std::pow(10., -10. + 0.2 * (irho + 0.5));
          const CoordinateVector<> uR(std::pow(10., -10. + 0.4 * (iu + 0.5)) -
                                      std::pow(10., 10. - 0.4 * (iu + 0.5)));
          const double PR = std::pow(10., -10. + 0.2 * (iP + 0.5));

          const CoordinateVector<> gradrhoL(
              2. * random_generator.get_uniform_random_double() - 1.,
              2. * random_generator.get_uniform_random_double() - 1.,
              2. * random_generator.get_uniform_random_double() - 1.);
          const CoordinateVector< CoordinateVector<> > graduL(
              CoordinateVector<>(
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.),
              CoordinateVector<>(
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.),
              CoordinateVector<>(
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.));
          const CoordinateVector<> gradPL(
              2. * random_generator.get_uniform_random_double() - 1.,
              2. * random_generator.get_uniform_random_double() - 1.,
              2. * random_generator.get_uniform_random_double() - 1.);

          const CoordinateVector<> gradrhoR =
              rhoR *
              CoordinateVector<>(
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.,
                  2. * random_generator.get_uniform_random_double() - 1.);
          const CoordinateVector< CoordinateVector<> > graduR(
              uR.x() *
                  CoordinateVector<>(
                      2. * random_generator.get_uniform_random_double() - 1.,
                      2. * random_generator.get_uniform_random_double() - 1.,
                      2. * random_generator.get_uniform_random_double() - 1.),
              uR.y() *
                  CoordinateVector<>(
                      2. * random_generator.get_uniform_random_double() - 1.,
                      2. * random_generator.get_uniform_random_double() - 1.,
                      2. * random_generator.get_uniform_random_double() - 1.),
              uR.z() *
                  CoordinateVector<>(
                      2. * random_generator.get_uniform_random_double() - 1.,
                      2. * random_generator.get_uniform_random_double() - 1.,
                      2. * random_generator.get_uniform_random_double() - 1.));
          const CoordinateVector<> gradPR =
              PR * CoordinateVector<>(
                       2. * random_generator.get_uniform_random_double() - 1.,
                       2. * random_generator.get_uniform_random_double() - 1.,
                       2. * random_generator.get_uniform_random_double() - 1.);

          HydroIntegrator::HydroFluxComputation::compute_fluxes(
              rhoL, uL, PL, rhoR, uR, PR, gradrhoL, graduL, gradPL, gradrhoR,
              graduR, gradPR, dL, dR, dL_over_r, dR_over_r, solver, normal,
              vframe, surface_area, timestep, true, mfluxL, pfluxL, EfluxL);
          HydroIntegrator::HydroFluxComputation::compute_fluxes(
              rhoR, uR, PR, rhoL, uL, PL, gradrhoR, graduR, gradPR, gradrhoL,
              graduL, gradPL, dR, dL, dR_over_r, dL_over_r, solver,
              -1. * normal, vframe, surface_area, timestep, true, mfluxR,
              pfluxR, EfluxR);

          assert_condition_message(mfluxL == -mfluxR, "%.17g %.17g", mfluxL,
                                   mfluxR);
          assert_condition_message(pfluxL == -1. * pfluxR,
                                   "[%.17g %.17g %.17g] [%.17g %.17g %.17g]",
                                   pfluxL.x(), pfluxL.y(), pfluxL.z(),
                                   pfluxR.x(), pfluxR.y(), pfluxR.z());
          assert_condition_message(EfluxL == -EfluxR, "%.17g %.17g", EfluxL,
                                   EfluxR);

          const double mLfluxlimit =
              random_generator.get_uniform_random_double();
          const double pL2fluxlimit =
              random_generator.get_uniform_random_double();
          const double ELfluxlimit =
              random_generator.get_uniform_random_double();
          const double mRfluxlimit =
              random_generator.get_uniform_random_double();
          const double pR2fluxlimit =
              random_generator.get_uniform_random_double();
          const double ERfluxlimit =
              random_generator.get_uniform_random_double();
          const double facL = HydroIntegrator::HydroFluxComputation::limit_flux(
              mfluxL, pfluxL.norm2(), EfluxL, mLfluxlimit, pL2fluxlimit,
              ELfluxlimit, mRfluxlimit, pR2fluxlimit, ERfluxlimit, true, true,
              false);
          const double facR = HydroIntegrator::HydroFluxComputation::limit_flux(
              mfluxR, pfluxR.norm2(), EfluxR, mRfluxlimit, pR2fluxlimit,
              ERfluxlimit, mLfluxlimit, pL2fluxlimit, ELfluxlimit, true, true,
              false);
          assert_condition(facL == facR);
        }
      }
    }
  }
  /// Flux limiter symmetry check
  {
    RandomGenerator random_generator(42);
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      const double mLfluxlimit = random_generator.get_uniform_random_double();
      const double pL2fluxlimit = random_generator.get_uniform_random_double();
      const double ELfluxlimit = random_generator.get_uniform_random_double();
      const double mRfluxlimit = random_generator.get_uniform_random_double();
      const double pR2fluxlimit = random_generator.get_uniform_random_double();
      const double ERfluxlimit = random_generator.get_uniform_random_double();
      const double mflux =
          2. * random_generator.get_uniform_random_double() - 1.;
      const double pflux2 = random_generator.get_uniform_random_double();
      const double Eflux =
          2. * random_generator.get_uniform_random_double() - 1.;
      const double facL = HydroIntegrator::HydroFluxComputation::limit_flux(
          mflux, pflux2, Eflux, mLfluxlimit, pL2fluxlimit, ELfluxlimit,
          mRfluxlimit, pR2fluxlimit, ERfluxlimit, true, true, false);
      const double facR = HydroIntegrator::HydroFluxComputation::limit_flux(
          -mflux, pflux2, -Eflux, mRfluxlimit, pR2fluxlimit, ERfluxlimit,
          mLfluxlimit, pL2fluxlimit, ELfluxlimit, true, true, false);
      assert_condition(facL == facR);
    }
  }

  /// Cartesian grid
  {
    HydroIntegrator integrator(5. / 3., false, false, 0.2, "HLLC", 0.);

    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    CoordinateVector< int_fast32_t > ncell(100, 1, 1);
    SodShockDensityFunction density_function;
    density_function.initialize();
    CoordinateVector< bool > periodic(true, true, true);
    CartesianDensityGrid grid(box, ncell, periodic, true);
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, density_function);

    integrator.initialize_hydro_variables(grid);

    const InternalHydroUnits *hydro_units = integrator.get_internal_units();
    assert_condition(hydro_units != nullptr);

    // write initial snapshot
    {
      std::ofstream snapfile("hydro_snap_0.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        const HydroVariables &hydro_vars = it.get_hydro_variables();

        snapfile << it.get_cell_midpoint().x() << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_DENSITY >(
                        hydro_vars.get_primitives_density())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_VELOCITY >(
                        hydro_vars.get_primitives_velocity().x())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_PRESSURE >(
                        hydro_vars.get_primitives_pressure())
                 << "\n";
        mtot += hydro_units->convert_to_SI_units< QUANTITY_MASS >(
            hydro_vars.get_conserved_mass());
        etot += hydro_units->convert_to_SI_units< QUANTITY_ENERGY >(
            hydro_vars.get_conserved_total_energy());
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }

    Timer serial_timer, parallel_timer;
    for (uint_fast8_t i = 0; i < 100; ++i) {
      integrator.do_hydro_step(grid, 0.001, serial_timer, parallel_timer);
    }

    // write final snapshot
    {
      std::ofstream snapfile("hydro_snap_1.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        const HydroVariables &hydro_vars = it.get_hydro_variables();

        snapfile << it.get_cell_midpoint().x() << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_DENSITY >(
                        hydro_vars.get_primitives_density())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_VELOCITY >(
                        hydro_vars.get_primitives_velocity().x())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_PRESSURE >(
                        hydro_vars.get_primitives_pressure())
                 << "\n";
        mtot += hydro_units->convert_to_SI_units< QUANTITY_MASS >(
            hydro_vars.get_conserved_mass());
        etot += hydro_units->convert_to_SI_units< QUANTITY_ENERGY >(
            hydro_vars.get_conserved_total_energy());
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }

    // output reference solution
    {
      std::ofstream snapfile("hydro_ref_1.txt");
      const ExactRiemannSolver solver(5. / 3.);
      const double t = 0.1;
      for (uint_fast32_t i = 0; i < 1000; ++i) {
        const double x = (i + 0.5) * 0.001;
        double rhosol, usol, Psol;
        solver.solve(1., 0., 1., 0.125, 0., 0.1, rhosol, usol, Psol,
                     (x - 0.5) / t);
        snapfile << x << "\t" << rhosol << "\t" << usol << "\t" << Psol << "\n";
      }
    }
  }

  /// Voronoi grid
  {
    HydroIntegrator integrator(5. / 3., false, false, 0.2, "HLLC", 0.);

    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    SodShockDensityFunction density_function;
    density_function.initialize();
    CoordinateVector< bool > periodic(false, false, false);
    OneDVoronoiGeneratorDistribution *generators =
        new OneDVoronoiGeneratorDistribution();
    VoronoiDensityGrid grid(generators, box, "Old", 0, periodic, true);
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, density_function);

    integrator.initialize_hydro_variables(grid);

    const InternalHydroUnits *hydro_units = integrator.get_internal_units();
    assert_condition(hydro_units != nullptr);

    // write initial snapshot
    {
      std::ofstream snapfile("hydro_voronoi_snap_0.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        const HydroVariables &hydro_vars = it.get_hydro_variables();

        snapfile << it.get_cell_midpoint().x() << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_DENSITY >(
                        hydro_vars.get_primitives_density())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_VELOCITY >(
                        hydro_vars.get_primitives_velocity().x())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_PRESSURE >(
                        hydro_vars.get_primitives_pressure())
                 << "\n";
        mtot += hydro_units->convert_to_SI_units< QUANTITY_MASS >(
            hydro_vars.get_conserved_mass());
        etot += hydro_units->convert_to_SI_units< QUANTITY_ENERGY >(
            hydro_vars.get_conserved_total_energy());
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }

    Timer serial_timer, parallel_timer;
    for (uint_fast8_t i = 0; i < 100; ++i) {
      integrator.do_hydro_step(grid, 0.001, serial_timer, parallel_timer);
    }

    // write final snapshot
    {
      std::ofstream snapfile("hydro_voronoi_snap_1.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        const HydroVariables &hydro_vars = it.get_hydro_variables();

        snapfile << it.get_cell_midpoint().x() << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_DENSITY >(
                        hydro_vars.get_primitives_density())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_VELOCITY >(
                        hydro_vars.get_primitives_velocity().x())
                 << "\t"
                 << hydro_units->convert_to_SI_units< QUANTITY_PRESSURE >(
                        hydro_vars.get_primitives_pressure())
                 << "\n";
        mtot += hydro_units->convert_to_SI_units< QUANTITY_MASS >(
            hydro_vars.get_conserved_mass());
        etot += hydro_units->convert_to_SI_units< QUANTITY_ENERGY >(
            hydro_vars.get_conserved_total_energy());
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }
  }

  return 0;
}
