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
 * @file HeliumTwoPhotonContinuumSpectrum.hpp
 *
 * @brief Helium 2-photon continuum ionizing PhotonSourceSpectrum.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HELIUMTWOPHOTONCONTINUUMSPECTRUM_HPP
#define HELIUMTWOPHOTONCONTINUUMSPECTRUM_HPP

#include "PhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"
#include <vector>

/*! @brief Number of frequencies in the internal table. */
#define HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ 1000

/**
 * @brief Helium 2-photon continuum photoionization spectrum.
 */
class HeliumTwoPhotonContinuumSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins (in 13.6 eV). */
  double _frequency[HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ];

  /*! @brief Cumulative distribution function. */
  double _cumulative_distribution[HELIUMTWOPHOTONCONTINUUMSPECTRUM_NUMFREQ];

  /*! @brief RandomGenerator used to generate random numbers. */
  RandomGenerator &_random_generator;

public:
  HeliumTwoPhotonContinuumSpectrum(RandomGenerator &random_generator);

  void get_spectrum(std::vector< double > &yHe2q, std::vector< double > &AHe2q);
  double get_integral(std::vector< double > &yHe2q,
                      std::vector< double > &AHe2q);

  virtual double get_random_frequency(double temperature = 0.);

  virtual double get_total_flux();
};

#endif // HELIUMTWOPHOTONCONTINUUMSPECTRUM_HPP
