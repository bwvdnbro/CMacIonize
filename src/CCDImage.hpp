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
 * @file CCDImage.hpp
 *
 * @brief CCD detector image.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CCDIMAGE_HPP
#define CCDIMAGE_HPP

#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "Utilities.hpp"

#include <fstream>
#include <string>
#include <vector>

/*! @brief Maximum number of levels that can be distinguished in an image. */
#define CCDIMAGE_NUM_LEVELS 255

/**
 * @brief CCD detector image.
 */
class CCDImage {
private:
  /*! @brief Direction of the observer. */
  const CoordinateVector<> _direction;

  /*! @brief Combined flux image. */
  std::vector< double > _image_total;

  /*! @brief Polarized image: Q component. */
  std::vector< double > _image_Q;

  /*! @brief Polarized image: U component. */
  std::vector< double > _image_U;

  /*! @brief Resolution of the image. */
  const unsigned int _resolution[2];

  /*! @brief Lower left corner of the image box (in kpc). */
  const double _anchor[2];

  /*! @brief Side lengths of the image box (in kpc). */
  const double _sides[2];

public:
  /**
   * @brief Constructor.
   *
   * @param direction Direction of the observer.
   * @param resolution_x Number of pixels in the horizontal direction.
   * @param resolution_y Number of pixels in the vertical direction.
   * @param anchor_x Left corner of the image box (in kpc).
   * @param anchor_y Lower corner of the image box (in kpc).
   * @param sides_x Horizontal side length of the image box (in kpc).
   * @param sides_y Vertical side length of the image box (in kpc).
   */
  inline CCDImage(const CoordinateVector<> direction, unsigned int resolution_x,
                  unsigned int resolution_y, double anchor_x, double anchor_y,
                  double sides_x, double sides_y)
      : _direction(direction), _resolution{resolution_x, resolution_y},
        _anchor{anchor_x, anchor_y}, _sides{sides_x, sides_y} {

    cmac_assert(_direction.norm2() == 1.);

    const double npixel = _resolution[0] * _resolution[1];
    _image_total.resize(npixel, 0.);
    _image_Q.resize(npixel, 0.);
    _image_U.resize(npixel, 0.);
  }

  /**
   * @brief Get the direction of the observer.
   *
   * @return Direction of the observer.
   */
  inline const CoordinateVector<> get_direction() const { return _direction; }

  /**
   * @brief Add a photon emitted from the given position and with the given
   * weights.
   *
   * @param position Position where the photon is emitted.
   * @param weight_total Total weight of the photon.
   * @param weight_Q Weight of the Q polarized component of the photon.
   * @param weight_U Weight of the U polarized component of the photon.
   */
  inline void add_photon(const CoordinateVector<> position, double weight_total,
                         double weight_Q, double weight_U) {

    const double po = std::atan2(_direction.y(), _direction.x());
    const double cospo = std::cos(po);
    const double sinpo = std::sin(po);
    const double costo = _direction.z();
    double sinto = std::sqrt(1. - costo * costo);
    if (cospo * _direction.x() > 0) {
      if (sinto * _direction.x() < 0) {
        sinto = -sinto;
      }
    } else {
      if (sinto * _direction.x() > 0) {
        sinto = -sinto;
      }
    }

    double xphoton = position.y() * cospo - position.x() * sinpo;
    double yphoton = position.z() * sinto - position.y() * costo * sinpo -
                     position.x() * costo * cospo;

    if (xphoton >= _anchor[0] && yphoton >= _anchor[1]) {
      xphoton -= _anchor[0];
      yphoton -= _anchor[1];
      if (xphoton < _sides[0] && yphoton < _sides[1]) {
        const unsigned int ix = (_resolution[0] * xphoton / _sides[0]);
        const unsigned int iy = (_resolution[1] * yphoton / _sides[1]);

        _image_total[ix * _resolution[1] + iy] += weight_total;
        _image_Q[ix * _resolution[1] + iy] += weight_Q;
        _image_U[ix * _resolution[1] + iy] += weight_U;
      }
    }
  }

  /**
   * @brief Save the image as a PGM image with the given name.
   *
   * @param filename Name of the image file.
   */
  inline void save(std::string filename) const {

    if (!Utilities::string_ends_with(filename, ".pgm")) {
      filename += ".pgm";
    }

    double min_value = _image_total[0];
    double max_value = _image_total[0];
    for (unsigned int i = 1; i < _image_total.size(); ++i) {
      min_value = std::min(min_value, _image_total[i]);
      max_value = std::max(max_value, _image_total[i]);
    }
    max_value -= min_value;

    std::ofstream image_file(filename);
    image_file << "P2\n"
               << _resolution[0] << " " << _resolution[1] << "\n"
               << CCDIMAGE_NUM_LEVELS << "\n";
    for (unsigned int iy = 0; iy < _resolution[1]; ++iy) {
      unsigned int pixelvalue = 0;
      if (max_value > 0.) {
        pixelvalue = std::round(CCDIMAGE_NUM_LEVELS *
                                (_image_total[iy] - min_value) / max_value);
      }
      image_file << pixelvalue;
      for (unsigned int ix = 1; ix < _resolution[0]; ++ix) {
        pixelvalue = 0;
        if (max_value > 0.) {
          pixelvalue = std::round(
              CCDIMAGE_NUM_LEVELS *
              (_image_total[ix * _resolution[1] + iy] - min_value) / max_value);
        }
        image_file << " " << pixelvalue;
      }
      image_file << "\n";
    }
    image_file.close();
  }
};

#endif // CCDIMAGE_HPP
