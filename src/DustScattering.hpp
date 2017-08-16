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
 * @file DustScattering.hpp
 *
 * @brief Dust scattering routines.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DUSTSCATTERING_HPP
#define DUSTSCATTERING_HPP

#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

#include <string>

class Photon;
class RandomGenerator;

/**
 * @brief Dust scattering routines.
 *
 * This file is based on:
 *  - White, R. L. 1979, ApJ, 229, 954
 *    (http://adsabs.harvard.edu/abs/1979ApJ...229..954W)
 *  - Code, Arthur D.; Whitney, Barbara A.  1995, ApJ, 441, 400
 *    (http://adsabs.harvard.edu/abs/1995ApJ...441..400C)
 *  - Witt, A. N. 1977, ApJS, 35, 1
 *    (http://adsabs.harvard.edu/abs/1977ApJS...35....1W)
 */
class DustScattering {
private:
  /*! @brief Henyey-Greenstein parameter @f$H_G@f$ for scattering. */
  const double _scattering_hgg;

  /*! @brief @f$H_G^2@f$. */
  const double _scattering_g2;

  /*! @brief @f$1-H_G^2@f$. */
  const double _scattering_omg2;

  /*! @brief @f$2H_G@f$. */
  const double _scattering_thgg;

  /*! @brief @f$1-H_G@f$. */
  const double _scattering_omhgg;

  /*! @brief @f$\frac{1}{2H_G}@f$. */
  const double _scattering_od2hgg;

  /*! @brief @f$1+H_G^2@f$. */
  const double _scattering_opg2;

  /*! @brief Peak linear polarization. */
  const double _scattering_pl;

  /*! @brief Asymmetry of the circular polarization. */
  const double _scattering_sc;

  /*! @brief Peak value of linear to circular polarization conversion. */
  const double _scattering_pc;

  /*! @brief Albedo of the dust. */
  const double _albedo;

  /*! @brief Dust attenuation coefficient (in m^2 kg^-1). */
  const double _kappa;

  /**
   * @brief Get the Henyey-Greenstein parameter for scattering corresponding to
   * the given band.
   *
   * @param band Band ("V" or "K").
   * @return Henyey-Greenstein parameter for scattering.
   */
  inline static double get_hgg_for_band(std::string band) {
    if (band == "V") {
      return 0.44;
    } else if (band == "K") {
      return 0.02;
    } else {
      cmac_error("Unknown band: %s!", band.c_str());
      return 0.;
    }
  }

  /**
   * @brief Get the peak linear polarization parameter corresponding to the
   * given band.
   *
   * @param band Band ("V" or "K").
   * @return Peak linear polarization parameter.
   */
  inline static double get_pl_for_band(std::string band) {
    if (band == "V") {
      return 0.43;
    } else if (band == "K") {
      return 0.93;
    } else {
      cmac_error("Unknown band: %s!", band.c_str());
      return 0.;
    }
  }

  /**
   * @brief Get the albedo corresponding to the given band.
   *
   * @param band Band ("V" or "K").
   * @return Albedo.
   */
  inline static double get_albedo_for_band(std::string band) {
    if (band == "V") {
      return 0.54;
    } else if (band == "K") {
      return 0.21;
    } else {
      cmac_error("Unknown band: %s!", band.c_str());
      return 0.;
    }
  }

  /**
   * @brief Get the dust attenuation coefficient corresponding to the given
   * band.
   *
   * @param band Band ("V" or "K").
   * @return Dust attenuation coefficient (in m^2 kg^-1).
   */
  inline static double get_kappa_for_band(std::string band) {
    if (band == "V") {
      return 21.9;
    } else if (band == "K") {
      return 2.;
    } else {
      cmac_error("Unknown band: %s!", band.c_str());
      return 0.;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param band Band in which photons are emitted ("V" or "K").
   * @param log Log to write logging info to.
   */
  inline DustScattering(std::string band, Log *log = nullptr)
      : _scattering_hgg(get_hgg_for_band(band)),
        _scattering_g2(_scattering_hgg * _scattering_hgg),
        _scattering_omg2(1. - _scattering_g2),
        _scattering_thgg(2. * _scattering_hgg),
        _scattering_omhgg(1. - _scattering_hgg),
        _scattering_od2hgg(0.5 / _scattering_hgg),
        _scattering_opg2(1. + _scattering_g2),
        _scattering_pl(get_pl_for_band(band)), _scattering_sc(1.),
        _scattering_pc(0.), _albedo(get_albedo_for_band(band)),
        _kappa(get_kappa_for_band(band)) {

    if (log) {
      log->write_status("Created DustScattering object for band ", band, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline DustScattering(ParameterFile &params, Log *log = nullptr)
      : DustScattering(params.get_value< std::string >("dust:band", "V"), log) {
  }

  /**
   * @brief Get the dust attenuation coefficient.
   *
   * @return Dust attenuation coefficient (in m^2 kg^-1).
   */
  inline double get_kappa() const { return _kappa; }

  /**
   * @brief Get the albedo.
   *
   * @return Albedo.
   */
  inline double get_albedo() const { return _albedo; }

  void scatter(Photon &photon, RandomGenerator &random_generator) const;
  double scatter_towards(Photon &photon, const CoordinateVector<> direction,
                         double sint, double cost, double phi, double sinp,
                         double cosp) const;
};

#endif // DUSTSCATTERING_HPP
