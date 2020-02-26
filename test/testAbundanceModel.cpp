/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testAbundanceModel.cpp
 *
 * @brief Unit test for the AbundanceModel class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "AbundanceModel.hpp"
#include "Assert.hpp"
#include "FixedValueAbundanceModel.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"
#include "SolarMetallicityAbundanceModel.hpp"

/**
 * @brief Unit test for the AbundanceModel class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // test FixedValueAbundanceModel
  {
    ParameterFile params;
    RandomGenerator random_generator(42);
    double abundances[NUMBER_OF_ELEMENTNAMES];
    for (int_fast32_t i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
      abundances[i] = random_generator.get_uniform_random_double();
      params.add_value("AbundanceModel:" + get_element_name(i),
                       std::to_string(abundances[i]));
    }

    FixedValueAbundanceModel fixed_value_model(params);
    const Abundances returned_abundances = fixed_value_model.get_abundances();
    for (int_fast32_t i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
      // since we mimic the actual parameter file construction process,
      // our random abundance values are first converted to a string,
      // which can lead to a loss of precision. That is why we do not
      // expect the values to match exactly, but allow for some relative
      // error.
      assert_values_equal_rel(returned_abundances.get_abundance(i),
                              abundances[i], 1.e-5);
    }
  }

  // test SolarMetallicityAbundanceModel
  {
    const double tolerance = 1.e-15;

    SolarMetallicityAbundanceModel solar_metallicity(-3.31);
    const Abundances solar_metallicity_abundances =
        solar_metallicity.get_abundances();
    SolarMetallicityAbundanceModel low_metallicity(-5.);
    const Abundances low_metallicity_abundances =
        low_metallicity.get_abundances();
    SolarMetallicityAbundanceModel high_metallicity(-3.);
    const Abundances high_metallicity_abundances =
        high_metallicity.get_abundances();
#ifdef HAS_HELIUM
    assert_values_equal_rel(
        solar_metallicity_abundances.get_abundance(ELEMENT_He),
        std::pow(10., -1.07), tolerance);
    assert_values_equal_rel(
        low_metallicity_abundances.get_abundance(ELEMENT_He),
        std::pow(10., -1.07), tolerance);
    assert_values_equal_rel(
        high_metallicity_abundances.get_abundance(ELEMENT_He),
        std::pow(10., -1.07), tolerance);
#endif
#ifdef HAS_CARBON
    assert_values_equal_rel(
        solar_metallicity_abundances.get_abundance(ELEMENT_C),
        std::pow(10., -3.57), tolerance);
    assert_values_equal_rel(low_metallicity_abundances.get_abundance(ELEMENT_C),
                            std::pow(10., -5.26), tolerance);
    assert_values_equal_rel(
        high_metallicity_abundances.get_abundance(ELEMENT_C),
        std::pow(10., -3.26), tolerance);
#endif
#ifdef HAS_NITROGEN
    assert_values_equal_rel(
        solar_metallicity_abundances.get_abundance(ELEMENT_N),
        std::pow(10., -4.17), tolerance);
    assert_values_equal_rel(low_metallicity_abundances.get_abundance(ELEMENT_N),
                            std::pow(10., -6.6), tolerance);
    assert_values_equal_rel(
        high_metallicity_abundances.get_abundance(ELEMENT_N),
        std::pow(10., -4.), tolerance);
#endif
#ifdef HAS_OXYGEN
    assert_values_equal_rel(
        solar_metallicity_abundances.get_abundance(ELEMENT_O),
        std::pow(10., -3.31), tolerance);
    assert_values_equal_rel(low_metallicity_abundances.get_abundance(ELEMENT_O),
                            1.e-5, tolerance);
    assert_values_equal_rel(
        high_metallicity_abundances.get_abundance(ELEMENT_O), 1.e-3, tolerance);
#endif
#ifdef HAS_NEON
    assert_values_equal_rel(
        solar_metallicity_abundances.get_abundance(ELEMENT_Ne),
        std::pow(10., -4.07), tolerance);
    assert_values_equal_rel(
        low_metallicity_abundances.get_abundance(ELEMENT_Ne),
        std::pow(10., -5.76), tolerance);
    assert_values_equal_rel(
        high_metallicity_abundances.get_abundance(ELEMENT_Ne),
        std::pow(10., -3.76), tolerance);
#endif
#ifdef HAS_SULPHUR
    assert_values_equal_rel(
        solar_metallicity_abundances.get_abundance(ELEMENT_S),
        std::pow(10., -4.88), tolerance);
    assert_values_equal_rel(low_metallicity_abundances.get_abundance(ELEMENT_S),
                            std::pow(10., -6.57), tolerance);
    assert_values_equal_rel(
        high_metallicity_abundances.get_abundance(ELEMENT_S),
        std::pow(10., -4.57), tolerance);
#endif

    // make sure code compiles if no elements are active
    (void)tolerance;
    (void)solar_metallicity_abundances;
    (void)low_metallicity_abundances;
    (void)high_metallicity_abundances;
  }

  return 0;
}
