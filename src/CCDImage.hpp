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
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Utilities.hpp"

#include <fstream>
#include <string>
#include <vector>

/*! @brief Maximum number of levels that can be distinguished in an image. */
#define CCDIMAGE_NUM_LEVELS 255

/**
 * @brief Types of CCD image output file.
 */
enum CCDImageType {
  /*! @brief Grayscale Netpbm format. */
  CCDIMAGETYPE_PGM,
  /*! @brief Binary floating point array. */
  CCDIMAGETYPE_BINARY_ARRAY
};

/**
 * @brief CCD detector image.
 */
class CCDImage {
private:
  /*! @brief Direction components: @f$\sin(\theta{})@f$, @f$\cos(\theta{})@f$,
   *  @f$\phi{}@f$ (in radians), @f$\sin(\phi{})@f$, and @f$\cos(\phi{})@f$. */
  const double _direction_parameters[5];

  /*! @brief Direction of the observer (@f$(\sin(\theta{})\cos(\phi{}),
   *  \sin(\theta{})\sin(\phi{}), \cos(\theta{}))@f$. */
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

  /*! @brief CCDImageType. */
  const CCDImageType _type;

  /*! @brief Name of the output file. */
  const std::string _filename;

  /**
   * @brief Get the CCDImageType corresponding to the given string.
   *
   * @param type std::string representation of a CCDImageType.
   * @return Corresponding CCDImageType.
   */
  inline static CCDImageType get_type(std::string type) {
    if (type == "PGM") {
      return CCDIMAGETYPE_PGM;
    } else if (type == "BinaryArray") {
      return CCDIMAGETYPE_BINARY_ARRAY;
    } else {
      cmac_error("Unknown image type: %s!", type.c_str());
      return CCDIMAGETYPE_BINARY_ARRAY;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param theta @f$\theta{}@f$ angle of the observer's direction (in radians).
   * @param phi @f$\phi{}@f$ angle of the observer's direction (in radians).
   * @param resolution_x Number of pixels in the horizontal direction.
   * @param resolution_y Number of pixels in the vertical direction.
   * @param anchor_x Left corner of the image box (in kpc).
   * @param anchor_y Lower corner of the image box (in kpc).
   * @param sides_x Horizontal side length of the image box (in kpc).
   * @param sides_y Vertical side length of the image box (in kpc).
   * @param type Type of image file.
   * @param filename Name of the output file.
   * @param output_folder Folder where the image is saved.
   * @param log Log to write logging info to.
   */
  inline CCDImage(double theta, double phi, unsigned int resolution_x,
                  unsigned int resolution_y, double anchor_x, double anchor_y,
                  double sides_x, double sides_y, std::string type,
                  std::string filename, std::string output_folder,
                  Log *log = nullptr)
      : _direction_parameters{std::sin(theta), std::cos(theta), phi,
                              std::sin(phi), std::cos(phi)},
        _direction(CoordinateVector<>(
            _direction_parameters[0] * _direction_parameters[4],
            _direction_parameters[0] * _direction_parameters[3],
            _direction_parameters[1])),
        _resolution{resolution_x, resolution_y}, _anchor{anchor_x, anchor_y},
        _sides{sides_x, sides_y}, _type(get_type(type)),
        _filename(Utilities::get_absolute_path(output_folder) + "/" +
                  filename) {

    const double npixel = _resolution[0] * _resolution[1];
    _image_total.resize(npixel, 0.);
    _image_Q.resize(npixel, 0.);
    _image_U.resize(npixel, 0.);

    if (log) {
      log->write_status(
          "Created emtpy CCDImage with a resolution of ", _resolution[0], "x",
          _resolution[1], ", with a viewing angle of (", theta, ", ", phi,
          "), and a viewport box with anchor [", _anchor[0], " m, ", _anchor[1],
          " m] and sides [", _sides[0], " m, ", _sides[1],
          " m]. The image will be saved under the name \"", _filename,
          "\", as an image of type ", type, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline CCDImage(ParameterFile &params, Log *log = nullptr)
      : CCDImage(
            params.get_physical_value< QUANTITY_ANGLE >("ccdimage:view_theta",
                                                        "89.7 degrees"),
            params.get_physical_value< QUANTITY_ANGLE >("ccdimage:view_phi",
                                                        "0. degrees"),
            params.get_value< unsigned int >("ccdimage:image_width", 200),
            params.get_value< unsigned int >("ccdimage:image_height", 200),
            params.get_physical_value< QUANTITY_LENGTH >("ccdimage:anchor_x",
                                                         "-12.1 kpc"),
            params.get_physical_value< QUANTITY_LENGTH >("ccdimage:anchor_y",
                                                         "-12.1 kpc"),
            params.get_physical_value< QUANTITY_LENGTH >("ccdimage:sides_x",
                                                         "24.2 kpc"),
            params.get_physical_value< QUANTITY_LENGTH >("ccdimage:sides_y",
                                                         "24.2 kpc"),
            params.get_value< std::string >("ccdimage:type", "BinaryArray"),
            params.get_value< std::string >("ccdimage:filename",
                                            "galaxy_image"),
            params.get_value< std::string >("output_folder", "."), log) {}

  /**
   * @brief Reset the image contents to zero.
   */
  inline void reset() {
    for (unsigned int i = 0; i < _image_total.size(); ++i) {
      _image_total[i] = 0.;
      _image_Q[i] = 0.;
      _image_U[i] = 0.;
    }
  }

  /**
   * @brief Get the direction of the observer (and the direction components).
   *
   * @param sin_theta Variable to store @f$\sin(\theta{})@f$ in.
   * @param cos_theta Variable to store @f$\cos(\theta{})@f$ in.
   * @param phi Variable to store @f$\phi{}@f$ in (in radians).
   * @param sin_phi Variable to store @f$\sin(\phi{})@f$ in.
   * @param cos_phi Variable to store @f$\cos(\phi{})@f$ in.
   * @return Direction of the observer.
   */
  inline const CoordinateVector<> get_direction(double &sin_theta,
                                                double &cos_theta, double &phi,
                                                double &sin_phi,
                                                double &cos_phi) const {
    sin_theta = _direction_parameters[0];
    cos_theta = _direction_parameters[1];
    phi = _direction_parameters[2];
    sin_phi = _direction_parameters[3];
    cos_phi = _direction_parameters[4];
    return _direction;
  }

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

    const double cospo = _direction_parameters[4];
    const double sinpo = _direction_parameters[3];
    const double costo = _direction_parameters[1];
    const double sinto = _direction_parameters[0];

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
   * @brief Add the given CCDImage to this image.
   *
   * @param image CCDImage to add.
   * @return Reference to the updated image.
   */
  inline CCDImage &operator+=(const CCDImage &image) {

    // make sure we are adding images of the same thing
    cmac_assert(_anchor[0] == image._anchor[0] &&
                _anchor[1] == image._anchor[1]);
    cmac_assert(_sides[0] == image._sides[0] && _sides[1] == image._sides[1]);
    cmac_assert(_resolution[0] == image._resolution[0] &&
                _resolution[1] == image._resolution[1]);
    cmac_assert(_direction == image._direction);
    cmac_assert(_image_total.size() == image._image_total.size());

    for (unsigned int i = 0; i < _image_total.size(); ++i) {
      _image_total[i] += image._image_total[i];
      _image_Q[i] += image._image_Q[i];
      _image_U[i] += image._image_U[i];
    }

    return *this;
  }

  /**
   * @brief Save the image as a file with the given name and type.
   *
   * @param normalization Normalization constant for the image.
   */
  inline void save(double normalization = 1.) const {

    std::string filename(_filename);

    if (_type == CCDIMAGETYPE_PGM) {

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
                (_image_total[ix * _resolution[1] + iy] - min_value) /
                max_value);
          }
          image_file << " " << pixelvalue;
        }
        image_file << "\n";
      }
      image_file.close();

    } else if (_type == CCDIMAGETYPE_BINARY_ARRAY) {

      if (!Utilities::string_ends_with(filename, ".dat")) {
        filename += ".dat";
      }

      std::vector< double > image_copy = _image_total;
      for (unsigned int i = 0; i < image_copy.size(); ++i) {
        image_copy[i] *= normalization;
      }

      std::ofstream array_file(filename, std::ios::binary);
      array_file.write(reinterpret_cast< const char * >(&image_copy[0]),
                       image_copy.size() * sizeof(double));
      array_file.close();

    } else {

      cmac_error("Unknown image type: %i!", _type);
    }
  }
};

#endif // CCDIMAGE_HPP
