/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file BimodalCrossSections.hpp
 *
 * @brief CrossSections implementation that uses two sets of fixed values below
 * and above a given frequency limit.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BIMODALCROSSSECTIONS_HPP
#define BIMODALCROSSSECTIONS_HPP

#include "CrossSections.hpp"
#include "ParameterFile.hpp"

/**
 * @brief CrossSections implementation that uses fixed values for all cross
 * sections.
 */
class BiModalCrossSections : public CrossSections {
private:
  /*! @brief Frequency limit (in Hz). */
  const double _frequency_limit;

  /*! @brief Cross section per ion below the limit (in m^2). */
  const double _low_cross_sections[NUMBER_OF_IONNAMES];

  /*! @brief Cross section per ion above the limit (in m^2). */
  const double _high_cross_sections[NUMBER_OF_IONNAMES];

public:
  /**
   * @brief Constructor.
   *
   * @param frequency_limit Frequency limit (in Hz).
   * @param cross_section_H_n_low Photoionization cross section of neutral
   * hydrogen below the limit (in m^2).
   * @param cross_section_H_n_high Photoionization cross section of neutral
   * hydrogen above the limit (in m^2).
   * @param cross_section_He_n_low Photoionization cross section of neutral
   * helium below the limit (in m^2).
   * @param cross_section_He_n_high Photoionization cross section of neutral
   * helium above the limit (in m^2).
   * @param cross_section_C_p1_low Photoionization cross section of single
   * ionized carbon below the limit (in m^2).
   * @param cross_section_C_p1_high Photoionization cross section of single
   * ionized carbon above the limit (in m^2).
   * @param cross_section_C_p2_low Photoionization cross section of double
   * ionized carbon below the limit (in m^2).
   * @param cross_section_C_p2_high Photoionization cross section of double
   * ionized carbon above the limit (in m^2).
   * @param cross_section_N_n_low Photoionization cross section of neutral
   * nitrogen below the limit (in m^2).
   * @param cross_section_N_n_high Photoionization cross section of neutral
   * nitrogen above the limit (in m^2).
   * @param cross_section_N_p1_low Photoionization cross section of single
   * ionized nitrogen below the limit (in m^2).
   * @param cross_section_N_p1_high Photoionization cross section of single
   * ionized nitrogen above the limit (in m^2).
   * @param cross_section_N_p2_low Photoionization cross section of double
   * ionized nitrogen below the limit (in m^2).
   * @param cross_section_N_p2_high Photoionization cross section of double
   * ionized nitrogen above the limit (in m^2).
   * @param cross_section_O_n_low Photoionization cross section of neutral
   * oxygen below the limit (in m^2).
   * @param cross_section_O_n_high Photoionization cross section of neutral
   * oxygen above the limit (in m^2).
   * @param cross_section_O_p1_low Photoionization cross section of single
   * ionized oxygen below the limit (in m^2).
   * @param cross_section_O_p1_high Photoionization cross section of single
   * ionized oxygen above the limit (in m^2).
   * @param cross_section_Ne_n_low Photoionization cross section of neutral neon
   * below the limit (in m^2).
   * @param cross_section_Ne_n_high Photoionization cross section of neutral
   * neon above the limit (in m^2).
   * @param cross_section_Ne_p1_low Photoionization cross section of single
   * ionized neon below the limit (in m^2).
   * @param cross_section_Ne_p1_high Photoionization cross section of single
   * ionized neon above the limit (in m^2).
   * @param cross_section_S_p1_low Photoionization cross section of single
   * ionized sulphur below the limit (in m^2).
   * @param cross_section_S_p1_high Photoionization cross section of single
   * ionized sulphur above the limit (in m^2).
   * @param cross_section_S_p2_low Photoionization cross section of double
   * ionized sulphur below the limit (in m^2).
   * @param cross_section_S_p2_high Photoionization cross section of double
   * ionized sulphur above the limit (in m^2).
   * @param cross_section_S_p3_low Photoionization cross section of triple
   * ionized sulphur below the limit (in m^2).
   * @param cross_section_S_p3_high Photoionization cross section of triple
   * ionized sulphur above the limit (in m^2).
   */
  inline BiModalCrossSections(
      const double frequency_limit, const double cross_section_H_n_low,
      const double cross_section_H_n_high, const double cross_section_He_n_low,
      const double cross_section_He_n_high, const double cross_section_C_p1_low,
      const double cross_section_C_p1_high, const double cross_section_C_p2_low,
      const double cross_section_C_p2_high, const double cross_section_N_n_low,
      const double cross_section_N_n_high, const double cross_section_N_p1_low,
      const double cross_section_N_p1_high, const double cross_section_N_p2_low,
      const double cross_section_N_p2_high, const double cross_section_O_n_low,
      const double cross_section_O_n_high, const double cross_section_O_p1_low,
      const double cross_section_O_p1_high, const double cross_section_Ne_n_low,
      const double cross_section_Ne_n_high,
      const double cross_section_Ne_p1_low,
      const double cross_section_Ne_p1_high,
      const double cross_section_S_p1_low, const double cross_section_S_p1_high,
      const double cross_section_S_p2_low, const double cross_section_S_p2_high,
      const double cross_section_S_p3_low, const double cross_section_S_p3_high)
      : _frequency_limit(frequency_limit),
        _low_cross_sections{cross_section_H_n_low,
#ifdef HAS_HELIUM
                            cross_section_He_n_low,
#endif
#ifdef HAS_CARBON
                            cross_section_C_p1_low,  cross_section_C_p2_low,
#endif
#ifdef HAS_NITROGEN
                            cross_section_N_n_low,   cross_section_N_p1_low,
                            cross_section_N_p2_low,
#endif
#ifdef HAS_OXYGEN
                            cross_section_O_n_high,  cross_section_O_p1_low,
#endif
#ifdef HAS_NEON
                            cross_section_Ne_n_low,  cross_section_Ne_p1_low,
#endif
#ifdef HAS_SULPHUR
                            cross_section_S_p1_high, cross_section_S_p2_low,
                            cross_section_S_p3_low
#endif
        },
        _high_cross_sections{cross_section_H_n_high,
#ifdef HAS_HELIUM
                             cross_section_He_n_high,
#endif
#ifdef HAS_CARBON
                             cross_section_C_p1_high, cross_section_C_p2_high,
#endif
#ifdef HAS_NITROGEN
                             cross_section_N_n_high,  cross_section_N_p1_high,
                             cross_section_N_p2_high,
#endif
#ifdef HAS_OXYGEN
                             cross_section_O_n_low,   cross_section_O_p1_high,
#endif
#ifdef HAS_NEON
                             cross_section_Ne_n_high, cross_section_Ne_p1_high,
#endif
#ifdef HAS_SULPHUR
                             cross_section_S_p1_low,  cross_section_S_p2_high,
                             cross_section_S_p3_high
#endif
        } {
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  inline BiModalCrossSections(ParameterFile &params)
      : BiModalCrossSections(
            params.get_physical_value< QUANTITY_FREQUENCY >("frequency limit:",
                                                            "15. eV"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:hydrogen_0_low", "6.3e-18 cm^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:hydrogen_0_high", "6.3e-18 cm^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:helium_0_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:helium_0_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:carbon_1_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:carbon_1_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:carbon_2_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:carbon_2_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_0_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_0_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_1_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_1_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_2_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:nitrogen_2_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:oxygen_0_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:oxygen_0_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:oxygen_1_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:oxygen_1_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:neon_0_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:neon_0_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:neon_1_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:neon_1_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_1_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_1_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_2_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_2_high", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_3_low", "0. m^2"),
            params.get_physical_value< QUANTITY_SURFACE_AREA >(
                "CrossSections:sulphur_3_high", "0. m^2")) {}

  /**
   * @brief Get the photoionization cross section for the given ion at the
   * given photon energy.
   *
   * @param ion IonName for a valid ion.
   * @param energy Photon frequency (in Hz).
   * @return Photoionization cross section (in m^2).
   */
  virtual double get_cross_section(const int_fast32_t ion,
                                   const double energy) const {
    if (energy < _frequency_limit) {
      return _low_cross_sections[ion];
    } else {
      return _high_cross_sections[ion];
    }
  }
};

#endif // BIMODALCROSSSECTIONS_HPP
