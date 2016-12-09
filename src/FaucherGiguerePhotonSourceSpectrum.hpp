/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file FaucherGiguerePhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum implementation for the Faucher-Giguere UVB
 * spectrum.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FAUCHERGIGUEREPHOTONSOURCESPECTRUM_HPP
#define FAUCHERGIGUEREPHOTONSOURCESPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"

#include <string>

class Log;
class ParameterFile;
class RandomGenerator;

/*! @brief Number of frequency bins. */
#define FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ 100

/**
 * @brief PhotonSourceSpectrum implementation for the Faucher-Giguere UVB
 * spectrum.
 */
class FaucherGiguerePhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  double _frequencies[FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Cumulative distribution of the spectrum. */
  double _cumulative_distribution[FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ];

  /*! @brief Total ionizing flux of the spectrum (in m^-2 s^-1). */
  double _total_flux;

public:
  FaucherGiguerePhotonSourceSpectrum(double redshift, Log *log = nullptr);

  FaucherGiguerePhotonSourceSpectrum(ParameterFile &params, Log *log = nullptr);

  static std::string get_filename(double z);

  virtual double get_total_flux();

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.);
};

#endif // FAUCHERGIGUEREPHOTONSOURCESPECTRUM_HPP
