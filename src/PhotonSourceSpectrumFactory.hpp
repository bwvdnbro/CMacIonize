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
 * @file PhotonSourceSpectrumFactory.hpp
 *
 * @brief Factory for PhotonSourceSpectrum instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCESPECTRUMFACTORY_HPP
#define PHOTONSOURCESPECTRUMFACTORY_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceSpectrum.hpp"

// implementations
#include "FaucherGiguerePhotonSourceSpectrum.hpp"
#include "MaskedPhotonSourceSpectrum.hpp"
#include "MonochromaticPhotonSourceSpectrum.hpp"
#include "Pegase3PhotonSourceSpectrum.hpp"
#include "PlanckPhotonSourceSpectrum.hpp"
#include "PopStarPhotonSourceSpectrum.hpp"
#include "UniformPhotonSourceSpectrum.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"

// HDF5 dependent implementations
#ifdef HAVE_HDF5
#include "CastelliKuruczPhotonSourceSpectrum.hpp"
#endif

/**
 * @brief Factory for PhotonSourceSpectrum instances.
 */
class PhotonSourceSpectrumFactory {
public:
  /**
   * @brief Method that checks if the requested PhotonSourceSpectrum
   * implementation requires HDF5.
   *
   * @param type Requested DensityFunction type.
   * @param log Log to write logging info to.
   */
  static void check_hdf5(std::string type, Log *log = nullptr) {
    if (type == "CastelliKurucz") {
      if (log) {
        log->write_error("Cannot create an instance of ", type,
                         "PhotonSourceSpectrum, since the code was "
                         "compiled without HDF5 support.");
      }
      cmac_error("A %sPhotonSourceSpectrum requires HDF5. However, the code "
                 "was compiled without HDF5 support!",
                 type.c_str());
    }
  }

  /**
   * @brief Generate a PhotonSourceSpectrum based on the given type.
   *
   * @param type Type of PhotonSourceSpectrum to generate.
   * @param role Role the PhotonSourceSpectrum will assume. Parameters will be
   * read from the corresponding parameter file block.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created PhotonSourceSpectrum instance. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  inline static PhotonSourceSpectrum *generate_from_type(std::string type,
                                                         std::string role,
                                                         ParameterFile &params,
                                                         Log *log = nullptr) {

#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif

    if (type == "FaucherGiguere") {
      return new FaucherGiguerePhotonSourceSpectrum(role, params, log);
    } else if (type == "Masked") {
      return new MaskedPhotonSourceSpectrum(role, params, log);
    } else if (type == "Monochromatic") {
      return new MonochromaticPhotonSourceSpectrum(role, params, log);
    } else if (type == "Pegase3") {
      return new Pegase3PhotonSourceSpectrum(role, params, log);
    } else if (type == "Planck") {
      return new PlanckPhotonSourceSpectrum(role, params, log);
    } else if (type == "PopStar") {
      return new PopStarPhotonSourceSpectrum(role, params, log);
    } else if (type == "Uniform") {
      return new UniformPhotonSourceSpectrum();
    } else if (type == "WMBasic") {
      return new WMBasicPhotonSourceSpectrum(role, params, log);
#ifdef HAVE_HDF5
    } else if (type == "CastelliKurucz") {
      return new CastelliKuruczPhotonSourceSpectrum(role, params, log);
#endif
    } else if (type == "None") {
      return nullptr;
    } else {
      cmac_error("Unknown PhotonSourceSpectrum type: \"%s\".", type.c_str());
      return nullptr;
    }
  }

  /**
   * @brief Generate a PhotonSourceSpectrum based on the type chosen in the
   * parameter file (corresponding to the given role).
   *
   * Supported types are (default: Planck):
   *  - CastelliKurucz: Castelli & Kurucz (2003) stellar atmosphere spectrum.
   *  - FaucherGiguere: Redshift dependent UVB spectrum of Faucher-Giguère et
   *    al. (2009)
   *  - Monochromatic: Monochromatic spectrum
   *  - Pegase3: Pegase 3 stellar spectra.
   *  - Planck: Black body spectrum
   *  - PopStar: PopStar stellar spectra.
   *  - WMBasic: Realistic stellar atmosphere spectrum of Sternberg, Hoffmann &
   *    Pauldrach (2003)
   *
   * @param role Role the PhotonSourceSpectrum will assume. Parameters will be
   * read from the corresponding parameter file block.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created PhotonSourceSpectrum instance. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  inline static PhotonSourceSpectrum *
  generate(std::string role, ParameterFile &params, Log *log = nullptr) {

    const std::string type =
        params.get_value< std::string >(role + ":type", "Monochromatic");
    if (log) {
      log->write_info("Requested PhotonSourceSpectrum for ", role, ": ", type);
    }
    return generate_from_type(type, role, params, log);
  }
};

#endif // PHOTONSOURCESPECTRUMFACTORY_HPP
