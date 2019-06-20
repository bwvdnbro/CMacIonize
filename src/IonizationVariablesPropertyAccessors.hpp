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
 * @file IonizationVariablesPropertyAccessors.hpp
 *
 * @brief PropertyAccessors for the IonizationVariables that should be
 * communicated over MPI.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONVARIABLESPROPERTYACCESSORS_HPP
#define IONIZATIONVARIABLESPROPERTYACCESSORS_HPP

#include "DensityGrid.hpp"

/**
 * @brief PropertyAccessor for the mean intensity integrals.
 */
template < IonName _ion_ > class MeanIntensityPropertyAccessor {
public:
  /**
   * @brief Get the mean intensity value for the template IonName.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @return Mean intensity value for that cell (in m^3).
   */
  inline static double get_value(const DensityGrid::iterator &it) {
    return it.get_ionization_variables().get_mean_intensity(_ion_);
  }

  /**
   * @brief Set the mean intensity value for the template IonName.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @param mean_intensity Mean intensity value for that cell (in m^3).
   */
  inline static void set_value(DensityGrid::iterator &it,
                               double mean_intensity) {
    it.get_ionization_variables().set_mean_intensity(_ion_, mean_intensity);
  }
};

/**
 * @brief PropertyAccessor for the ionic fractions.
 */
template < IonName _ion_ > class IonicFractionPropertyAccessor {
public:
  /**
   * @brief Get the ionic fraction for the template IonName.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @return Ionic fraction for that cell.
   */
  inline static double get_value(const DensityGrid::iterator &it) {
    return it.get_ionization_variables().get_ionic_fraction(_ion_);
  }

  /**
   * @brief Set the ionic fraction for the template IonName.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @param ionic_fraction Ionic fraction for that cell.
   */
  inline static void set_value(DensityGrid::iterator &it,
                               double ionic_fraction) {
    it.get_ionization_variables().set_ionic_fraction(_ion_, ionic_fraction);
  }
};

/**
 * @brief PropertyAccessor for the heating.
 */
template < HeatingTermName _heating_term_ > class HeatingPropertyAccessor {
public:
  /**
   * @brief Get the heating for the template HeatingTermName.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @return Heating for that cell (in m^3 s^-1).
   */
  inline static double get_value(const DensityGrid::iterator &it) {
    return it.get_ionization_variables().get_heating(_heating_term_);
  }

  /**
   * @brief Set the heating for the template HeatingTermName.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @param heating Heating for that cell (in m^3 s^-1).
   */
  inline static void set_value(DensityGrid::iterator &it, double heating) {
    it.get_ionization_variables().set_heating(_heating_term_, heating);
  }
};

/**
 * @brief PropertyAccessor for the temperature.
 */
class TemperaturePropertyAccessor {
public:
  /**
   * @brief Get the temperature value.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @return Temperature of that cell (in K).
   */
  inline static double get_value(const DensityGrid::iterator &it) {
    return it.get_ionization_variables().get_temperature();
  }

  /**
   * @brief Set the temperature value.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @param temperature Temperature of that cell (in K).
   */
  inline static void set_value(DensityGrid::iterator &it, double temperature) {
    it.get_ionization_variables().set_temperature(temperature);
  }
};

/**
 * @brief PropertyAccessor for the number density.
 */
class NumberDensityPropertyAccessor {
public:
  /**
   * @brief Get the number density value.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @return Number density of that cell (in m^-3).
   */
  inline static double get_value(const DensityGrid::iterator &it) {
    return it.get_ionization_variables().get_number_density();
  }

  /**
   * @brief Set the number density value.
   *
   * @param it DensityGrid::iterator pointing to a cell.
   * @param number_density Number density of that cell (in m^-3).
   */
  inline static void set_value(DensityGrid::iterator &it,
                               double number_density) {
    it.get_ionization_variables().set_number_density(number_density);
  }
};

#endif // IONIZATIONVARIABLESPROPERTYACCESSORS_HPP
