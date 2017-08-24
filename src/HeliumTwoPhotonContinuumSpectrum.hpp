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
 *
 * The data used comes from Drake, G. W. F., Victor, G. A. & Dalgarno, A. 1969,
 * PhRev, 180, 25 (http://adsabs.harvard.edu/abs/1969PhRv..180...25D), table II.
 * We took into account the symmetry in the variable \f$y\f$ to double the
 * number of rows in that table and stored the data values for helium in an
 * external data file.
 */
class HeliumTwoPhotonContinuumSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins (in 13.6 eV). */
  std::vector< double > _frequency;

  /*! @brief Cumulative distribution function. */
  std::vector< double > _cumulative_distribution;

public:
  HeliumTwoPhotonContinuumSpectrum();

  /**
   * @brief Virtual destructor.
   */
  virtual ~HeliumTwoPhotonContinuumSpectrum() {}

  void get_spectrum(std::vector< double > &yHe2q,
                    std::vector< double > &AHe2q) const;
  double get_integral(std::vector< double > &yHe2q,
                      std::vector< double > &AHe2q) const;

  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature = 0.) const;

  virtual double get_total_flux() const;
};

#endif // HELIUMTWOPHOTONCONTINUUMSPECTRUM_HPP
