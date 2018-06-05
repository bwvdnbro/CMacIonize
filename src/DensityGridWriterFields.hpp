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
 * @file DensityGridWriterFields.hpp
 *
 * @brief Variable field selection for DensityGrid output.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDWRITERFIELDS_HPP
#define DENSITYGRIDWRITERFIELDS_HPP

#include "DensityGrid.hpp"
#include "Error.hpp"
#include "ParameterFile.hpp"

#include <cinttypes>
#include <string>

/**
 * @brief Convenient indices for all supported output fields.
 */
enum DensityGridField {
  DENSITYGRIDFIELD_COORDINATES = 0,
  DENSITYGRIDFIELD_NUMBER_DENSITY,
  DENSITYGRIDFIELD_TEMPERATURE,
  DENSITYGRIDFIELD_NEUTRAL_FRACTION,
  DENSITYGRIDFIELD_COOLING,
  DENSITYGRIDFIELD_PHOTOIONIZATION_RATE,
  DENSITYGRIDFIELD_HEATING_RATE,
  DENSITYGRIDFIELD_COSMIC_RAY_FACTOR,
  DENSITYGRIDFIELD_DENSITY,
  DENSITYGRIDFIELD_VELOCITIES,
  DENSITYGRIDFIELD_PRESSURE,
  DENSITYGRIDFIELD_MASS,
  DENSITYGRIDFIELD_TOTAL_ENERGY,
  DENSITYGRIDFIELD_NUMBER
};

/**
 * @brief Convenient indices for all supported output field types.
 */
enum DensityGridFieldType {
  DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE = 0,
  DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE,
  DENSITYGRIDFIELDTYPE_NUMBER
};

/**
 * @brief Functionality to handle DensityGridFields and DensityGridFieldTypes.
 */
class DensityGridWriterFields {
  /// STATIC FUNCTIONS
public:
  /**
   * @brief Get the DensityGridFieldType corresponding to the given
   * DensityGridField.
   *
   * @param field_name DensityGridField.
   * @return Corresponding DensityGridFieldType.
   */
  inline static int_fast32_t get_type(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELD_COORDINATES:
      return DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE;
    case DENSITYGRIDFIELD_NUMBER_DENSITY:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_TEMPERATURE:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_NEUTRAL_FRACTION:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_COOLING:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_PHOTOIONIZATION_RATE:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_HEATING_RATE:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_COSMIC_RAY_FACTOR:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_DENSITY:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_VELOCITIES:
      return DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE;
    case DENSITYGRIDFIELD_PRESSURE:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_MASS:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    case DENSITYGRIDFIELD_TOTAL_ENERGY:
      return DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return DENSITYGRIDFIELDTYPE_NUMBER;
    }
  }

  /**
   * @brief Get the name corresponding to the given DensityGridField.
   *
   * @param field_name DensityGridField.
   * @return Corresponding human readable name.
   */
  inline static std::string get_name(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELD_COORDINATES:
      return "Coordinates";
    case DENSITYGRIDFIELD_NUMBER_DENSITY:
      return "NumberDensity";
    case DENSITYGRIDFIELD_TEMPERATURE:
      return "Temperature";
    case DENSITYGRIDFIELD_NEUTRAL_FRACTION:
      return "NeutralFraction";
    case DENSITYGRIDFIELD_COOLING:
      return "Cooling";
    case DENSITYGRIDFIELD_PHOTOIONIZATION_RATE:
      return "PhotoionizationRate";
    case DENSITYGRIDFIELD_HEATING_RATE:
      return "HeatingRate";
    case DENSITYGRIDFIELD_COSMIC_RAY_FACTOR:
      return "CosmicRayFactor";
    case DENSITYGRIDFIELD_DENSITY:
      return "Density";
    case DENSITYGRIDFIELD_VELOCITIES:
      return "Velocities";
    case DENSITYGRIDFIELD_PRESSURE:
      return "Pressure";
    case DENSITYGRIDFIELD_MASS:
      return "Mass";
    case DENSITYGRIDFIELD_TOTAL_ENERGY:
      return "TotalEnergy";
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return "";
    }
  }

  /**
   * @brief Is the given DensityGridField an ion property?
   *
   * @param field_name DensityGridField.
   * @return True if the given field has a value for each ion.
   */
  inline static bool is_ion_property(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELD_COORDINATES:
      return false;
    case DENSITYGRIDFIELD_NUMBER_DENSITY:
      return false;
    case DENSITYGRIDFIELD_TEMPERATURE:
      return false;
    case DENSITYGRIDFIELD_NEUTRAL_FRACTION:
      return true;
    case DENSITYGRIDFIELD_COOLING:
      return true;
    case DENSITYGRIDFIELD_PHOTOIONIZATION_RATE:
      return true;
    case DENSITYGRIDFIELD_HEATING_RATE:
      return false;
    case DENSITYGRIDFIELD_COSMIC_RAY_FACTOR:
      return false;
    case DENSITYGRIDFIELD_DENSITY:
      return false;
    case DENSITYGRIDFIELD_VELOCITIES:
      return false;
    case DENSITYGRIDFIELD_PRESSURE:
      return false;
    case DENSITYGRIDFIELD_MASS:
      return false;
    case DENSITYGRIDFIELD_TOTAL_ENERGY:
      return false;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return false;
    }
  }

  /**
   * @brief Is the given DensityGridField a heating term property?
   *
   * @param field_name DensityGridField.
   * @return True if the given field has a value for each heating term.
   */
  inline static bool is_heating_property(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELD_COORDINATES:
      return false;
    case DENSITYGRIDFIELD_NUMBER_DENSITY:
      return false;
    case DENSITYGRIDFIELD_TEMPERATURE:
      return false;
    case DENSITYGRIDFIELD_NEUTRAL_FRACTION:
      return false;
    case DENSITYGRIDFIELD_COOLING:
      return false;
    case DENSITYGRIDFIELD_PHOTOIONIZATION_RATE:
      return false;
    case DENSITYGRIDFIELD_HEATING_RATE:
      return true;
    case DENSITYGRIDFIELD_COSMIC_RAY_FACTOR:
      return false;
    case DENSITYGRIDFIELD_DENSITY:
      return false;
    case DENSITYGRIDFIELD_VELOCITIES:
      return false;
    case DENSITYGRIDFIELD_PRESSURE:
      return false;
    case DENSITYGRIDFIELD_MASS:
      return false;
    case DENSITYGRIDFIELD_TOTAL_ENERGY:
      return false;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return false;
    }
  }

  /**
   * @brief Get the value of the double scalar field with the given name.
   *
   * @param field_name DensityGridField.
   * @param it DensityGrid::iterator to a cell.
   * @return Value of the double scalar field.
   */
  inline static double
  get_scalar_double_value(const int_fast32_t field_name,
                          const DensityGrid::iterator &it) {
    switch (field_name) {
    case DENSITYGRIDFIELD_NUMBER_DENSITY:
      return it.get_ionization_variables().get_number_density();
    case DENSITYGRIDFIELD_TEMPERATURE:
      return it.get_ionization_variables().get_temperature();
    case DENSITYGRIDFIELD_COSMIC_RAY_FACTOR:
      return it.get_ionization_variables().get_cosmic_ray_factor();
    case DENSITYGRIDFIELD_DENSITY:
      return it.get_hydro_variables().get_primitives_density();
    case DENSITYGRIDFIELD_PRESSURE:
      return it.get_hydro_variables().get_primitives_pressure();
    case DENSITYGRIDFIELD_MASS:
      return it.get_hydro_variables().get_conserved_mass();
    case DENSITYGRIDFIELD_TOTAL_ENERGY:
      return it.get_hydro_variables().get_conserved_total_energy();
    default:
      cmac_error("Not a scalar DensityGridField: %" PRIiFAST32, field_name);
      return 0.;
    }
  }

  /**
   * @brief Get the value of the double vector field with the given name.
   *
   * @param field_name DensityGridField.
   * @param it DensityGrid::iterator to a cell.
   * @param box_anchor Anchor of the simulation box, used for vector corrections
   * (in m).
   * @return Value of the double vector field.
   */
  inline static CoordinateVector<>
  get_vector_double_value(const int_fast32_t field_name,
                          const DensityGrid::iterator &it,
                          const CoordinateVector<> box_anchor) {
    switch (field_name) {
    case DENSITYGRIDFIELD_COORDINATES:
      return it.get_cell_midpoint() - box_anchor;
    case DENSITYGRIDFIELD_VELOCITIES:
      return it.get_hydro_variables().get_primitives_velocity();
    default:
      cmac_error("Not a vector DensityGridField: %" PRIiFAST32, field_name);
      return CoordinateVector<>();
    }
  }

  /**
   * @brief Get the value of the double scalar ion field with the given name.
   *
   * @param field_name DensityGridField.
   * @param ion_name IonName.
   * @param it DensityGrid::iterator to a cell.
   * @return Value of the double scalar ion field.
   */
  inline static double
  get_scalar_double_ion_value(const int_fast32_t field_name,
                              const int_fast32_t ion_name,
                              const DensityGrid::iterator &it) {
    switch (field_name) {
    case DENSITYGRIDFIELD_NEUTRAL_FRACTION:
      return it.get_ionization_variables().get_ionic_fraction(ion_name);
    case DENSITYGRIDFIELD_COOLING:
#ifdef DO_OUTPUT_COOLING
      return it.get_ionization_variables().get_cooling(ion_name);
#else
      cmac_error("Code was not configured with cooling output support!");
      return 0.;
#endif
    case DENSITYGRIDFIELD_PHOTOIONIZATION_RATE:
      return it.get_ionization_variables().get_mean_intensity(ion_name);
    default:
      cmac_error("Not a scalar ion DensityGridField: %" PRIiFAST32, field_name);
      return 0.;
    }
  }

  /**
   * @brief Get the value of the double scalar heating property field with the
   * given name.
   *
   * @param field_name DensityGridField.
   * @param heating_property_name HeatingTermName.
   * @param it DensityGrid::iterator to a cell.
   * @return Value fo the double scalar heating property field.
   */
  inline static double
  get_scalar_double_heating_value(const int_fast32_t field_name,
                                  const int_fast32_t heating_property_name,
                                  const DensityGrid::iterator &it) {
    switch (field_name) {
    case DENSITYGRIDFIELD_HEATING_RATE:
      return it.get_ionization_variables().get_heating(heating_property_name);
    default:
      cmac_error("Not a scalar heating property DensityGridField: %" PRIiFAST32,
                 field_name);
      return false;
    }
  }

private:
  /*! @brief Field flag for each DensityGridField. For fields with multiple
   *  variables this flag encodes the variables that should be written to the
   *  output. */
  uint_least8_t _field_flag[DENSITYGRIDFIELD_NUMBER];

public:
  /**
   * @brief Default constructor.
   *
   * @param hydro Output hydro variables?
   */
  inline DensityGridWriterFields(const bool hydro) {

    _field_flag[DENSITYGRIDFIELD_COORDINATES] = 1;
    _field_flag[DENSITYGRIDFIELD_NUMBER_DENSITY] = 1;
    _field_flag[DENSITYGRIDFIELD_TEMPERATURE] = 0;
    _field_flag[DENSITYGRIDFIELD_NEUTRAL_FRACTION] = 1;
    _field_flag[DENSITYGRIDFIELD_COOLING] = 0;
    _field_flag[DENSITYGRIDFIELD_PHOTOIONIZATION_RATE] = 0;
    _field_flag[DENSITYGRIDFIELD_HEATING_RATE] = 0;
    _field_flag[DENSITYGRIDFIELD_COSMIC_RAY_FACTOR] = 0;
    _field_flag[DENSITYGRIDFIELD_DENSITY] = hydro;
    _field_flag[DENSITYGRIDFIELD_VELOCITIES] = hydro;
    _field_flag[DENSITYGRIDFIELD_PRESSURE] = hydro;
    _field_flag[DENSITYGRIDFIELD_MASS] = 0;
    _field_flag[DENSITYGRIDFIELD_TOTAL_ENERGY] = 0;
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  inline DensityGridWriterFields(ParameterFile &params) {

    for (uint_fast32_t property = 0; property < DENSITYGRIDFIELDTYPE_NUMBER;
         ++property) {
      const std::string prop_name = get_name(property);
      if (is_ion_property(property)) {
        _field_flag[property] = 0;
        for (uint_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
          const std::string name = prop_name + get_ion_name(ion);
          _field_flag[property] += params.get_value< uint_fast32_t >(
                                       "DensityGridWriterFields:" + name, false)
                                   << ion;
        }
      } else if (is_heating_property(property)) {
        _field_flag[property] = 0;
        for (uint_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
             ++heating) {
          const std::string name = prop_name + get_ion_name(heating);
          _field_flag[property] += params.get_value< uint_fast32_t >(
                                       "DensityGridWriterFields:" + name, false)
                                   << heating;
        }
      } else {
        // normal single variable property
        _field_flag[property] = params.get_value< bool >(
            "DensityGridWriterFields:" + prop_name, false);
      }
    }
  }
};

#endif // DENSITYGRIDWRITERFIELDS_HPP
