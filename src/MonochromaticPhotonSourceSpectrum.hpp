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
 * @file MonochromaticPhotonSourceSpectrum.hpp
 *
 * @brief PhotonSourceSpectrum that returns photons of a single frequency.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MONOCHROMATICPHOTONSOURCESPECTRUM_HPP
#define MONOCHROMATICPHOTONSOURCESPECTRUM_HPP

#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceSpectrum.hpp"

/**
 * @brief PhotonSourceSpectrum that returns photons of a single frequency.
 */
class MonochromaticPhotonSourceSpectrum : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency of emitted photons (in Hz). */
  const double _frequency;

  /*! @brief Total flux of emitted photons (in m^-2 s^-1). */
  const double _total_flux;

public:
  /**
   * @brief Constructor
   *
   * @param frequency Frequency of emitted photons (in Hz).
   * @param total_flux Total flux of emitted photons (in m^-2 s^-1).
   * @param log Log to write logging info to.
   */
  MonochromaticPhotonSourceSpectrum(double frequency, double total_flux = -1,
                                    Log *log = nullptr)
      : _frequency(frequency), _total_flux(total_flux) {

    if (log) {
      log->write_status(
          "Created MonochromaticPhotonSourceSpectrum with frequency ",
          _frequency, " Hz and total flux ", _total_flux, " m^-2 s^-1.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - frequency: Constant monochromatic frequency of the spectrum (default:
   *    13.6 eV)
   *  - total flux: Total flux of the spectrum (default: none, and error when
   *    total flux is requested)
   *
   * @param role Role the spectrum will fulfil in the simulation. Parameters are
   * read from the corresponding block in the parameter file.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  MonochromaticPhotonSourceSpectrum(std::string role, ParameterFile &params,
                                    Log *log = nullptr)
      : MonochromaticPhotonSourceSpectrum(
            params.get_physical_value< QUANTITY_FREQUENCY >(role + ":frequency",
                                                            "13.6 eV"),
            params.get_physical_value< QUANTITY_FLUX >(role + ":total flux",
                                                       "-1. m^-2 s^-1"),
            log) {}

  /**
   * @brief Get a random frequency from this spectrum.
   *
   * Since we only have one frequency in our spectrum, this frequency will
   * always be the same "random" frequency.
   *
   * @param random_generator RandomGenerator that is not used.
   * @param temperature Temperature that is completely ignored.
   * @return Constant frequency of this spectrum.
   */
  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature) const {
    return _frequency;
  }

  /**
   * @brief Get the total ionizing flux emitted by this spectrum.
   *
   * @return Total ionizing flux (in m^-2 s^-1).
   */
  virtual double get_total_flux() const {
    if (_total_flux < 0.) {
      cmac_error("This function should not be used!");
    }
    return _total_flux;
  }
};

#endif // MONOCHROMATICPHOTONSOURCESPECTRUM_HPP
