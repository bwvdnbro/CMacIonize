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
#include <vector>

class Log;
class ParameterFile;
class RandomGenerator;

/*! @brief Number of frequency bins. */
#define FAUCHERGIGUEREPHOTONSOURCESPECTRUM_NUMFREQ 100

/**
 * @brief PhotonSourceSpectrum implementation for the Faucher-Giguere UVB
 * spectrum.
 *
 * The spectrum data comes from Faucher-Giguère, C.-A., Lidz, A., Zaldarriaga,
 * M. & Hernquist, L., 2009, ApJ, 703, 1416
 * (http://adsabs.harvard.edu/abs/2009ApJ...703.1416F) and was downloaded from
 * http://galaxies.northwestern.edu/files/2013/09/fg_uvb_dec11.zip.
 *
 * The Faucher-Giguère et al. (2009) data contains non-normalized spectra for
 * all redshifts in the range [0. - 10.65] in increments of 0.05. We precompute
 * the spectrum for a specific redshift by linearly interpolating on these
 * spectra for the requested redshift and integrate it out to get the total
 * ionizing luminosity.
 */
class FaucherGiguerePhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  std::vector< double > _frequencies;

  /*! @brief Cumulative distribution of the spectrum. */
  std::vector< double > _cumulative_distribution;

  /*! @brief Total ionizing flux of the spectrum (in m^-2 s^-1). */
  double _total_flux;

public:
  FaucherGiguerePhotonSourceSpectrum(double redshift, Log *log = nullptr);

  FaucherGiguerePhotonSourceSpectrum(std::string role, ParameterFile &params,
                                     Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~FaucherGiguerePhotonSourceSpectrum() {}

  static std::string get_filename(double z);

  virtual double get_total_flux() const;

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const;
};

#endif // FAUCHERGIGUEREPHOTONSOURCESPECTRUM_HPP
