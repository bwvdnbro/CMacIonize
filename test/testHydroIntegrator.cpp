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
#include "CartesianDensityGrid.hpp"
#include "DensityFunction.hpp"
#include "HydroIntegrator.hpp"
#include "RiemannSolver.hpp"
#include "VoronoiDensityGrid.hpp"
#include "VoronoiGeneratorDistribution.hpp"
#include <fstream>

/**
 * @brief DensityFunction implementation that sets up a basic Sod shock problem.
 */
class SodShockDensityFunction : public DensityFunction {
public:
  /**
   * @brief Get the density field at the given position.
   *
   * @param position Position.
   * @return DensityValues at that position.
   */
  virtual DensityValues operator()(CoordinateVector<> position) const {
    const double hydrogen_mass = 1.6737236e-27;
    const double boltzmann_k = 1.38064852e-23;
    double density_unit = 1. / hydrogen_mass;
    double temperature_unit = hydrogen_mass / boltzmann_k;
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
  unsigned int _last_index;

public:
  /**
   * @brief Get the number of positions returned by this distribution.
   *
   * @return 100.
   */
  virtual unsigned int get_number_of_positions() const { return 100; }

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
    HydroIntegrator integrator(5. / 3., false);

    Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    CoordinateVector< int > ncell(100, 1, 1);
    SodShockDensityFunction density_function;
    CoordinateVector< bool > periodic(false, true, true);
    CartesianDensityGrid grid(box, ncell, density_function, periodic, true);
    std::pair< unsigned long, unsigned long > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block);

    integrator.initialize_hydro_variables(grid);

    // write initial snapshot
    {
      std::ofstream snapfile("hydro_snap_0.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        snapfile << it.get_cell_midpoint().x() << "\t"
                 << it.get_hydro_primitive_density() << "\t"
                 << it.get_hydro_primitive_velocity_x() << "\t"
                 << it.get_hydro_primitive_pressure() << "\n";
        mtot += it.get_hydro_conserved_mass();
        etot += it.get_hydro_conserved_total_energy();
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }

    for (unsigned int i = 0; i < 100; ++i) {
      integrator.do_hydro_step(grid, 0.001);
    }

    // write final snapshot
    {
      std::ofstream snapfile("hydro_snap_1.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        snapfile << it.get_cell_midpoint().x() << "\t"
                 << it.get_hydro_primitive_density() << "\t"
                 << it.get_hydro_primitive_velocity_x() << "\t"
                 << it.get_hydro_primitive_pressure() << "\n";
        mtot += it.get_hydro_conserved_mass();
        etot += it.get_hydro_conserved_total_energy();
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }

    // output reference solution
    {
      std::ofstream snapfile("hydro_ref_1.txt");
      const RiemannSolver solver(5. / 3.);
      const double t = 0.1;
      for (unsigned int i = 0; i < 1000; ++i) {
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
    HydroIntegrator integrator(5. / 3., false);

    Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
    SodShockDensityFunction density_function;
    CoordinateVector< bool > periodic(false, false, false);
    OneDVoronoiGeneratorDistribution *generators =
        new OneDVoronoiGeneratorDistribution();
    VoronoiDensityGrid grid(generators, density_function, box, 0, periodic,
                            true, 0.001, 5. / 3.);
    std::pair< unsigned long, unsigned long > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block);

    integrator.initialize_hydro_variables(grid);

    // write initial snapshot
    {
      std::ofstream snapfile("hydro_voronoi_snap_0.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        snapfile << it.get_cell_midpoint().x() << "\t"
                 << it.get_hydro_primitive_density() << "\t"
                 << it.get_hydro_primitive_velocity_x() << "\t"
                 << it.get_hydro_primitive_pressure() << "\n";
        mtot += it.get_hydro_conserved_mass();
        etot += it.get_hydro_conserved_total_energy();
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }

    for (unsigned int i = 0; i < 100; ++i) {
      integrator.do_hydro_step(grid, 0.001);
    }

    // write final snapshot
    {
      std::ofstream snapfile("hydro_voronoi_snap_1.txt");
      double mtot = 0.;
      double etot = 0.;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        snapfile << it.get_cell_midpoint().x() << "\t"
                 << it.get_hydro_primitive_density() << "\t"
                 << it.get_hydro_primitive_velocity_x() << "\t"
                 << it.get_hydro_primitive_pressure() << "\n";
        mtot += it.get_hydro_conserved_mass();
        etot += it.get_hydro_conserved_total_energy();
      }
      cmac_status("Total mass: %g, total energy: %g", mtot, etot);
    }
  }

  return 0;
}
