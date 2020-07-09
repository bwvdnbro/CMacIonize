/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SpectrumTracker.hpp
 *
 * @brief Instrument that keeps track of the spectrum of photons that pass
 * through a cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPECTRUMTRACKER_HPP
#define SPECTRUMTRACKER_HPP

#include "Photon.hpp"
#include "PhotonPacket.hpp"
#include "Tracker.hpp"
#include "YAMLDictionary.hpp"

#include <cinttypes>
#include <fstream>
#include <string>
#include <vector>

/**
 * @brief Instrument that keeps track of the spectrum of photons that pass
 * through a cell.
 */
class SpectrumTracker : public Tracker {
private:
  /*! @brief Minimum frequency for the spectral bins (in Hz). */
  const double _minimum_frequency;

  /*! @brief Width of a single frequency bin (in Hz). */
  const double _frequency_width;

  /*! @brief Inverse width of a single frequency bin (in Hz^-1). */
  const double _inverse_frequency_width;

  /*! @brief Cosine of the opening angle needed for detection. */
  const double _cos_opening_angle;

  /*! @brief Reference direction for the detection opening angle. */
  const CoordinateVector<> _reference_direction;

  /*! @brief Number counts per bin for direct radiation from a source. */
  std::vector< uint_fast64_t > _number_counts_primary;

  /*! @brief Number counts per bin for diffuse hydrogen reemission. */
  std::vector< uint_fast64_t > _number_counts_diffuse_H;

  /*! @brief Number counts per bin for diffuse helium reemission. */
  std::vector< uint_fast64_t > _number_counts_diffuse_He;

  /**
   * @brief Normalize the given direction.
   *
   * If the given direction is [0., 0., 0.], we do nothing.
   *
   * @param direction Input direction.
   * @return Normalized output direction.
   */
  inline static CoordinateVector<>
  normalize_direction(const CoordinateVector<> direction) {
    const double nrm2 = direction.norm2();
    if (nrm2 > 0.) {
      return direction / std::sqrt(nrm2);
    } else {
      return direction;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_bins Number of bins to use.
   * @param opening_angle Opening angle for detection. Photon packets with
   * directions that are larger than this angle (w.r.t. the reference direction)
   * are not detected (in radians).
   * @param reference_direction Reference direction.
   */
  SpectrumTracker(
      const uint_fast32_t number_of_bins, const double opening_angle = 180.,
      const CoordinateVector<> reference_direction = CoordinateVector<>(0.))
      : _minimum_frequency(3.289e15),
        _frequency_width(3. * 3.289e15 / number_of_bins),
        _inverse_frequency_width(1. / _frequency_width),
        _cos_opening_angle(std::cos(opening_angle)),
        _reference_direction(normalize_direction(reference_direction)),
        _number_counts_primary(number_of_bins, 0),
        _number_counts_diffuse_H(number_of_bins, 0),
        _number_counts_diffuse_He(number_of_bins, 0) {}

  /**
   * @brief YAMLDictionary constructor.
   *
   * @param name Name of the block in the dictionary that contains additional
   * parameters for the spectrum tracker.
   * @param blocks YAMLDictionary that contains additional parameters.
   */
  SpectrumTracker(const std::string name, YAMLDictionary &blocks)
      : SpectrumTracker(
            blocks.get_value< uint_fast32_t >(name + "number of bins", 100),
            blocks.get_physical_value< QUANTITY_ANGLE >(name + "opening angle",
                                                        "180. degrees"),
            blocks.get_value< CoordinateVector<> >(name + "reference direction",
                                                   CoordinateVector<>(0.))) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~SpectrumTracker() {}

  /**
   * @brief Make a duplicate of the current tracker.
   *
   * @return Pointer to a new duplicate of the tracker.
   */
  virtual Tracker *duplicate() {
    return new SpectrumTracker(_number_counts_primary.size(),
                               std::acos(_cos_opening_angle),
                               _reference_direction);
  }

  /**
   * @brief Add the contribution from the given duplicate tracker to this
   * tracker.
   *
   * @param tracker Duplicate tracker (created using Tracker::duplicate()).
   */
  virtual void merge(Tracker *tracker) {
    const SpectrumTracker *other =
        reinterpret_cast< SpectrumTracker * >(tracker);
    for (uint_fast32_t i = 0; i < _number_counts_primary.size(); ++i) {
      _number_counts_primary[i] += other->_number_counts_primary[i];
      _number_counts_diffuse_H[i] += other->_number_counts_diffuse_H[i];
      _number_counts_diffuse_He[i] += other->_number_counts_diffuse_He[i];
    }
  }

  /**
   * @brief Add the contribution of the given photon packet to the bins.
   *
   * @param photon Photon to add.
   */
  virtual void count_photon(const Photon &photon) {

    if (_reference_direction.norm2() > 0.) {
      const double dot_product = CoordinateVector<>::dot_product(
          photon.get_direction(), _reference_direction);
      if (dot_product < _cos_opening_angle) {
        // photon is travelling in the wrong direction
        return;
      }
    }

    const double frequency = photon.get_energy();
    const uint_fast32_t index =
        (frequency - _minimum_frequency) * _inverse_frequency_width;
    if (index < _number_counts_primary.size()) {
      if (photon.get_type() == PHOTONTYPE_PRIMARY) {
        ++_number_counts_primary[index];
      }
      if (photon.get_type() == PHOTONTYPE_DIFFUSE_HI) {
        ++_number_counts_diffuse_H[index];
      }
      if (photon.get_type() == PHOTONTYPE_DIFFUSE_HeI) {
        ++_number_counts_diffuse_He[index];
      }
    }
  }

  /**
   * @brief Add the contribution of the given photon packet to the bins.
   *
   * @param photon Photon to add.
   * @param absorption Absorption counters within the cell for this photon
   * (in m^-1).
   */
  virtual void count_photon(const PhotonPacket &photon,
                            const double *absorption) {

    if (_reference_direction.norm2() > 0.) {
      const double dot_product = CoordinateVector<>::dot_product(
          photon.get_direction(), _reference_direction);
      if (dot_product < _cos_opening_angle) {
        // photon is travelling in the wrong direction
        return;
      }
    }

    const double frequency = photon.get_energy();
    const uint_fast32_t index =
        (frequency - _minimum_frequency) * _inverse_frequency_width;
    if (index < _number_counts_primary.size()) {
      if (photon.get_type() == PHOTONTYPE_PRIMARY) {
        ++_number_counts_primary[index];
      }
      if (photon.get_type() == PHOTONTYPE_DIFFUSE_HI) {
        ++_number_counts_diffuse_H[index];
      }
      if (photon.get_type() == PHOTONTYPE_DIFFUSE_HeI) {
        ++_number_counts_diffuse_He[index];
      }
    }
  }

  /**
   * @brief Output the spectrum to the file with the given name.
   *
   * @param filename Name of the output file.
   */
  virtual void output_tracker(const std::string filename) const {

    std::ofstream ofile(filename);
    ofile << "# frequency (Hz)\tprimary count\tdiffuse H count\tdiffuse He "
             "count\n";
    for (uint_fast32_t i = 0; i < _number_counts_primary.size(); ++i) {
      const double nu = _minimum_frequency + (i + 0.5) * _frequency_width;
      ofile << nu << "\t" << _number_counts_primary[i] << "\t"
            << _number_counts_diffuse_H[i] << "\t"
            << _number_counts_diffuse_He[i] << "\n";
    }
  }

  /**
   * @brief Describe the tracker in the given output stream, appending the given
   * prefix to each line of output.
   *
   * @param prefix Prefix to add to each output line.
   * @param stream std::ostream to write to.
   */
  virtual void describe(const std::string prefix, std::ostream &stream) const {

    stream << prefix << "type: Spectrum\n";
    stream << prefix << "number of bins: " << _number_counts_primary.size()
           << "\n";
    stream << prefix << "opening angle: " << std::acos(_cos_opening_angle)
           << " radians\n";
    stream << prefix << "reference direction: [" << _reference_direction.x()
           << ", " << _reference_direction.y() << ", "
           << _reference_direction.z() << "]\n";
  }
};

#endif // SPECTRUMTRACKER_HPP
