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
 * @file SolarMetallicityAbundanceModel.hpp
 *
 * @brief AbundanceModel implementation that contains fixed values for all
 * abundances, based on a single metallicity parameter.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SOLARMETALLICITYABUNDANCEMODEL_HPP
#define SOLARMETALLICITYABUNDANCEMODEL_HPP

#include "AbundanceModel.hpp"
#include "ElementNames.hpp"
#include "ParameterFile.hpp"

/**
 * @brief AbundanceModel implementation that contains fixed values for all
 * abundances, based on a single metallicity parameter.
 *
 * The values are based on Asplund et al. (2009)
 * (https://ui.adsabs.harvard.edu/abs/2009ARA&A..47..481A/abstract). We
 * parametrise them using the oxygen abundance and scale all abundances (except
 * helium and nitrogren) up or down with the same factor.
 *
 * The helium abundance is kept fixed. The nitrogen abundance is varied
 * according to Stasinska et al. (2006)
 * (https://ui.adsabs.harvard.edu/abs/2006MNRAS.371..972S/abstract), eqs
 * (3)-(4):
 * @f[
 * \log_{10}(N/O) =
 * \begin{cases}
 *   -1.6 & \log_{10}(O/H) \leq{} -4, \\
 *   0.6 \left( \log_{10}(O/H) + 4 \right) -1.6 & \log_{10}(O/H) > -4.
 * \end{cases}
 * @f]
 */
class SolarMetallicityAbundanceModel : public AbundanceModel {
private:
  /*! @brief Per element abundance. */
  Abundances _abundances;

public:
  /**
   * @brief Constructor.
   *
   * @param metallicity Metallicity of the gas (given as the logarithm of the
   * oxygen abundance, @f$\log_{10}(O/H)@f$).
   */
  inline SolarMetallicityAbundanceModel(const double metallicity) {

    // based on a table provided by Natalia Vale Asari (private
    // communication)
    const double solar_He = -1.07;
    const double solar_C = -3.57;
    const double solar_N = -4.17;
    const double solar_O = -3.31;
    const double solar_Ne = -4.07;
    const double solar_S = -4.88;

    double actual_C, actual_N, actual_Ne, actual_S;
    if (metallicity != solar_O) {
      const double Odiff = metallicity - solar_O;
      actual_C = solar_C + Odiff;
      actual_Ne = solar_Ne + Odiff;
      actual_S = solar_S + Odiff;
      if (metallicity <= -4.) {
        actual_N = metallicity - 1.6;
      } else {
        actual_N = metallicity + 0.6 * (metallicity + 4.) - 1.6;
      }
    } else {
      actual_C = solar_C;
      actual_N = solar_N;
      actual_Ne = solar_Ne;
      actual_S = solar_S;
    }

#ifdef HAS_HELIUM
    _abundances.set_abundance(ELEMENT_He, std::pow(10., solar_He));
#endif
#ifdef HAS_CARBON
    _abundances.set_abundance(ELEMENT_C, std::pow(10., actual_C));
#endif
#ifdef HAS_NITROGEN
    _abundances.set_abundance(ELEMENT_N, std::pow(10., actual_N));
#endif
#ifdef HAS_OXYGEN
    _abundances.set_abundance(ELEMENT_O, std::pow(10., metallicity));
#endif
#ifdef HAS_NEON
    _abundances.set_abundance(ELEMENT_Ne, std::pow(10., actual_Ne));
#endif
#ifdef HAS_SULPHUR
    _abundances.set_abundance(ELEMENT_S, std::pow(10., actual_S));
#endif

    // make sure code compiles if no elements are active
    (void)solar_He;
    (void)actual_C;
    (void)actual_N;
    (void)actual_Ne;
    (void)actual_S;
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read from the file:
   *  - metallicity: Metallicity of the gas (default: -3.31).
   *
   * @param params ParameterFile to read from.
   */
  inline SolarMetallicityAbundanceModel(ParameterFile &params)
      : SolarMetallicityAbundanceModel(
            params.get_value< double >("AbundanceModel:metallicity", -3.31)) {}

  virtual ~SolarMetallicityAbundanceModel() {}

  /**
   * @brief Get the abundances for all cells in the simulation.
   *
   * @return Abundance values for all cells in the simulation.
   */
  virtual const Abundances get_abundances() const { return _abundances; }

  /**
   * @brief Get the abundance values for the given cell.
   *
   * @param cell Cell containing geometrical information about a cell.
   * @return Abundances for that cell.
   */
  virtual const Abundances get_abundances(const Cell &cell) const {
    return _abundances;
  }
};

#endif // SOLARMETALLICITYABUNDANCEMODEL_HPP
