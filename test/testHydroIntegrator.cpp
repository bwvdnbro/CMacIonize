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
  virtual DensityValues operator()(const Cell &cell) const {
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

  /// Cartesian grid
  {
    HydroIntegrator integrator(5. / 3., false, false, 0.2, "HLLC", 0.);

    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    CoordinateVector< int_fast32_t > ncell(100, 1, 1);
    SodShockDensityFunction density_function;
    density_function.initialize();
    CoordinateVector< bool > periodic(false, true, true);
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
