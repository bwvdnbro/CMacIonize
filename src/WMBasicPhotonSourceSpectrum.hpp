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
 * @file WMBasicPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for the Pauldrach, Hoffmann &
 * Lennon (2001) stellar model spectra.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef WMBASICPHOTONSOURCESPECTRUM_HPP
#define WMBASICPHOTONSOURCESPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

#include <string>

class Log;
class ParameterFile;
class RandomGenerator;

/**
 * @brief Number of frequency bins used in the internal table.
 */
#define WMBASICPHOTONSOURCESPECTRUM_NUMFREQ 1000

/**
 * @brief PhotonSourceSpectrum implementation for the Pauldrach, Hoffmann &
 * Lennon (2001) (http://adsabs.harvard.edu/abs/2001A%26A...375..161P) stellar
 * model spectra.
 *
 * The data used comes from Sternberg, Hoffmann & Pauldrach (2003)
 * (http://adsabs.harvard.edu/abs/2003IAUS..212..206H), and was downloaded from
 * the link in that paper (ftp://wise3.tau.ac.il/pub/stars).
 *
 * The model spectra are based on detailed models of stellar atmospheres and
 * take into account key processes like radiation driven winds. They should be
 * representative for the atmospheres of young O stars.
 *
 * The spectra depend on two parameters: the effective temperature of the star,
 * and its surface gravity. We only have data tables for a limited number of
 * parameter values: temperatures in the range [25,000 K; 50,000 K] (every 1,000
 * K), and surface gravities in the range \f$\left[ \log_{10}\left(\frac{300}
 * {{\rm{} cm\,s^{-2}}}\right), \log_{10}\left(\frac{400}{{\rm{} cm\,s^{-2}}}
 * \right)\right]\f$ (every \f$\log_{10}\left(\frac{20}{{\rm{} cm\,s^{-2}}}
 * \right) \f$). Some combinations are missing. If tables are requested for a
 * non-existent combination of parameter values, the code will crash.
 *
 * We resample the spectra on a frequency grid of 1000 bins in the range
 * [13.6 eV, 54.4 eV]. This smooths out some of the features in the spectrum,
 * but works pretty well overall (as can be visually confirmed by looking at the
 * file wmbasic.txt produced by testWMBasicPhotonSourceSpectrum).
 */
class WMBasicPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  double _frequencies[WMBASICPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Cumulative distribution of the spectrum. */
  double _cumulative_distribution[WMBASICPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Total ionizing flux of the spectrum (in m^-2 s^-1). */
  double _total_flux;

public:
  WMBasicPhotonSourceSpectrum(double temperature, double surface_gravity,
                              Log *log = nullptr);

  WMBasicPhotonSourceSpectrum(std::string role, ParameterFile &params,
                              Log *log = nullptr);

  static std::string get_log_g_name(double surface_gravity);
  static std::string get_filename(double temperature, double surface_gravity);

  /**
   * @brief Virtual destructor.
   */
  virtual ~WMBasicPhotonSourceSpectrum() {}

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const;

  virtual double get_total_flux() const;
};

#endif // WMBASICPHOTONSOURCESPECTRUM_HPP
