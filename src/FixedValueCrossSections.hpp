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
 * @file FixedValueCrossSections.hpp
 *
 * @brief CrossSections implementation that uses fixed values for all cross
 * sections.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FIXEDVALUECROSSSECTIONS_HPP
#define FIXEDVALUECROSSSECTIONS_HPP

#include "CrossSections.hpp"
#include "ParameterFile.hpp"

/**
 * @brief CrossSections implementation that uses fixed values for all cross
 * sections.
 */
class FixedValueCrossSections : public CrossSections {
private:
  /*! @brief Cross section per ion (in m^2). */
  const double _cross_sections[NUMBER_OF_IONNAMES];

public:
  /**
   * @brief Constructor.
   *
   * @param cross_section_H_n Photoionization cross section of neutral hydrogen
   * (in m^2).
   * @param cross_section_He_n Photoionization cross section of neutral helium
   * (in m^2).
   * @param cross_section_C_p1 Photoionization cross section of single ionized
   * carbon (in m^2).
   * @param cross_section_C_p2 Photoionization cross section of double ionized
   * carbon (in m^2).
   * @param cross_section_N_n Photoionization cross section of neutral nitrogen
   * (in m^2).
   * @param cross_section_N_p1 Photoionization cross section of single ionized
   * nitrogen (in m^2).
   * @param cross_section_N_p2 Photoionization cross section of double ionized
   * nitrogen (in m^2).
   * @param cross_section_O_n Photoionization cross section of neutral oxygen
   * (in m^2).
   * @param cross_section_O_p1 Photoionization cross section of single ionized
   * oxygen (in m^2).
   * @param cross_section_Ne_n Photoionization cross section of neutral neon
   * (in m^2).
   * @param cross_section_Ne_p1 Photoionization cross section of single ionized
   * neon (in m^2).
   * @param cross_section_S_p1 Photoionization cross section of single ionized
   * sulphur (in m^2).
   * @param cross_section_S_p2 Photoionization cross section of double ionized
   * sulphur (in m^2).
   * @param cross_section_S_p3 Photoionization cross section of triple ionized
   * sulphur (in m^2).
   */
  FixedValueCrossSections(double cross_section_H_n, double cross_section_He_n,
                          double cross_section_C_p1, double cross_section_C_p2,
                          double cross_section_N_n, double cross_section_N_p1,
                          double cross_section_N_p2, double cross_section_O_n,
                          double cross_section_O_p1, double cross_section_Ne_n,
                          double cross_section_Ne_p1, double cross_section_S_p1,
                          double cross_section_S_p2, double cross_section_S_p3)
      : _cross_sections{
            cross_section_H_n,  cross_section_He_n,  cross_section_C_p1,
            cross_section_C_p2, cross_section_N_n,   cross_section_N_p1,
            cross_section_N_p2, cross_section_O_n,   cross_section_O_p1,
            cross_section_Ne_n, cross_section_Ne_p1, cross_section_S_p1,
            cross_section_S_p2, cross_section_S_p3} {}

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  FixedValueCrossSections(ParameterFile &params)
      : FixedValueCrossSections(
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:hydrogen_0", "6.3e-18 cm^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:helium_0", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:carbon_1", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:carbon_2", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_0", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_1", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_2", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:oxygen_0", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:oxygen_1", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:neon_0", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:neon_1", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_1", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_2", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_3", "0. m^2")) {}

  /**
   * @brief Get the photoionization cross section for the given ion at the
   * given photon energy.
   *
   * @param ion IonName for a valid ion.
   * @param energy Photon frequency (in Hz).
   * @return Photoionization cross section (in m^2).
   */
  virtual double get_cross_section(IonName ion, double energy) const {
    return _cross_sections[ion];
  }
};

#endif // FIXEDVALUECROSSSECTIONS_HPP
