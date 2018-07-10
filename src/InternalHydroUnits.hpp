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
 * @file InternalHydroUnits.hpp
 *
 * @brief Internal unit system for hydrodynamics.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef INTERNALHYDROUNITS_HPP
#define INTERNALHYDROUNITS_HPP

#include "UnitConverter.hpp"

/**
 * @brief Internal unit system for hydrodynamics.
 *
 * The idea is to rescale all hydrodynamical quantities so that the average
 * values are 1 initially. This can help to reduce round off error problems in
 * astrophysical setups with very low density values in SI units.
 *
 * The hydrodynamical unit system is completely determined by setting a length,
 * time and mass unit, @f$U_L@f$, @f$U_t@f$ and @f$U_M@f$ respectively. The
 * corresponding units of density (@f$U_{\rho{}}@f$), velocity (@f$U_v@f$),
 * pressure (@f$U_P@f$), momentum (@f$U_p@f$) and energy (@f$U_E@f$) are then
 * given by
 * \f[
 *   U_{\rho{}} = \frac{U_M}{U_L^3},
 * \f]
 * \f[
 *   U_v = \frac{U_L}{U_t},
 * \f]
 * \f[
 *   U_P = \frac{U_M}{U_L U_t^2},
 * \f]
 * \f[
 *   U_p = \frac{U_M U_L}{U_t}
 * \f]
 * and
 * \f[
 *   U_E = \frac{U_M U_L^2}{U_t^2}.
 * \f]
 *
 * Our rescaling requirement can be expressed as
 * \f[
 *   \langle{} \rho{} \rangle{} = U_{\rho{}}
 * \f]
 * and
 * \f[
 *   \langle{} P \rangle{} = U_P,
 * \f]
 * where @f$\langle{} x \rangle{}@f$ represents the average value of quantity
 * @f$x@f$ over all cells. Note that the ratio of the density and pressure is
 * proportional to the sound speed and hence already defines the velocity unit;
 * we cannot impose a similar constraint on the velocity. To constrain the
 * remaining free unit, we will instead rescale the length unit so that the
 * average box dimensions are 1 as well.
 *
 * The entire procedure is now well defined: we first find the average box
 * size @f$L_{\rm{}box}@f$ among the 3 coordinate dimensions, and set the length
 * unit accordingly:
 * \f[
 *   U_L = L_{\rm{}box}.
 * \f]
 * We then find the average density and pressure over all cells, and use the
 * density constraint to find the mass unit:
 * \f[
 *   U_M = \langle{} \rho{} \rangle{} U_L^3.
 * \f]
 * Finally, we use the pressure constraint to set the time unit
 * \f[
 *   U_t = \sqrt{\frac{\langle{} P \rangle{} U_L}{U_M}}.
 * \f]
 *
 * Note that for now, we only perform this rescaling for the hydrodynamical
 * quantities; all other quantities (number densities, temperatures, reaction
 * rates...) are still expressed in SI units. For this reason, we will not
 * actually rescale the length variables. This means we have to add some
 * additional unit conversions in the hydro related methods that use length
 * variables, as well as in the conversion from hydrodynamical quantities to
 * photoionization quantities.
 */
class InternalHydroUnits {
private:
  /*! @brief Conversion factors from SI units to internal units. */
  double _internal_conversion_factors[NUMBER_OF_QUANTITIES];

  /*! @brief Conversion factors from internal units to SI units. */
  double _SI_conversion_factors[NUMBER_OF_QUANTITIES];

  /**
   * @brief Compute conversion factors from and to the internal unit system.
   *
   * @param average_box_size Average box length @f$L_{\rm{}box}@f$ (in m).
   * @param average_density Average density @f$\langle{} \rho{} \rangle{}@f$
   * (in kg m^-3).
   * @param average_pressure Average pressure @f$\langle{} P \rangle{}@f$
   * (in kg m^-1 s^-2).
   * @param to_internal Array with conversion factors to internal units.
   * @param to_SI Array with conversion factors to SI units.
   */
  static inline void compute_conversion_factors(
      const double average_box_size, const double average_density,
      const double average_pressure, double to_internal[NUMBER_OF_QUANTITIES],
      double to_SI[NUMBER_OF_QUANTITIES]) {

    const double length_unit_in_SI = average_box_size;
    const double mass_unit_in_SI = average_density * length_unit_in_SI *
                                   length_unit_in_SI * length_unit_in_SI;
    const double time_unit_in_SI =
        std::sqrt(mass_unit_in_SI / (average_pressure * length_unit_in_SI));

    const Unit length_unit(length_unit_in_SI, 1, 0, 0, 0, 0, 0);
    const Unit time_unit(time_unit_in_SI, 0, 1, 0, 0, 0, 0);
    const Unit mass_unit(mass_unit_in_SI, 0, 0, 1, 0, 0, 0);
    const Unit temperature_unit(1., 0, 0, 0, 1, 0, 0);
    const Unit current_unit(1., 0, 0, 0, 0, 1, 0);
    const Unit angle_unit(1., 0, 0, 0, 0, 0, 1);

    for (int quantity = 0; quantity < NUMBER_OF_QUANTITIES; ++quantity) {
      const Unit SI_unit = UnitConverter::get_SI_unit(quantity);
      const Unit internal_unit =
          SI_unit.get_equivalent(length_unit, time_unit, mass_unit,
                                 temperature_unit, current_unit, angle_unit);
      to_internal[quantity] = 1. / internal_unit;
      to_SI[quantity] = 1. * internal_unit;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param average_box_size Average box length @f$L_{\rm{}box}@f$ (in m).
   * @param average_density Average density @f$\langle{} \rho{} \rangle{}@f$
   * (in kg m^-3).
   * @param average_pressure Average pressure @f$\langle{} P \rangle{}@f$
   * (in kg m^-1 s^-2).
   */
  InternalHydroUnits(const double average_box_size,
                     const double average_density,
                     const double average_pressure) {

    compute_conversion_factors(average_box_size, average_density,
                               average_pressure, _internal_conversion_factors,
                               _SI_conversion_factors);
  }

  /**
   * @brief Reset the internal unit system.
   *
   * @param average_box_size Average box length @f$L_{\rm{}box}@f$ (in m).
   * @param average_density Average density @f$\langle{} \rho{} \rangle{}@f$
   * (in kg m^-3).
   * @param average_pressure Average pressure @f$\langle{} P \rangle{}@f$
   * (in kg m^-1 s^-2).
   */
  inline void reset_units(const double average_box_size,
                          const double average_density,
                          const double average_pressure) {

    compute_conversion_factors(average_box_size, average_density,
                               average_pressure, _internal_conversion_factors,
                               _SI_conversion_factors);
  }

  /**
   * @brief Convert the given variable representing the given template quantity
   * to internal units.
   *
   * @param value Value to convert (in SI units).
   * @return Value in internal units.
   */
  template < Quantity _quantity_ >
  inline double convert_to_internal_units(const double value) const {
    return value * _internal_conversion_factors[_quantity_];
  }

  /**
   * @brief Convert the given variable representing the given template quantity
   * to SI units.
   *
   * @param value Value to convert (in internal units).
   * @return Value in SI units.
   */
  template < Quantity _quantity_ >
  inline double convert_to_SI_units(const double value) const {
    return value * _SI_conversion_factors[_quantity_];
  }

  /**
   * @brief Get the SI value of the internal unit corresponding to the given
   * template quantity.
   *
   * @return SI value of the internal unit corresponding to the quantity.
   */
  template < Quantity _quantity_ > inline double get_unit_SI_value() const {
    return _SI_conversion_factors[_quantity_];
  }

  /**
   * @brief Get the internal value of the SI unit corresponding to the given
   * template quantity.
   *
   * @return Internal value of the SI unit corresponding to the quantity.
   */
  template < Quantity _quantity_ >
  inline double get_unit_internal_value() const {
    return _internal_conversion_factors[_quantity_];
  }
};

#endif // INTERNALHYDROUNITS_HPP
