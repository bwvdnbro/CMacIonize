/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file WeightedSpectrumTracker.hpp
 *
 * @brief Instrument that keeps track of the spectrum of photons that pass
 * through a cell.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef WEIGHTEDSPECTRUMTRACKER_HPP
#define WEIGHTEDSPECTRUMTRACKER_HPP

#include "Cell.hpp"
#include "Photon.hpp"
#include "PhotonPacket.hpp"
#include "Tracker.hpp"
#include "YAMLDictionary.hpp"

#include <cinttypes>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Instrument that keeps track of the spectrum of photons that pass
 * through a cell.
 *
 * This version weighs all contributions with the projected surface area of the
 * incoming photon packets to provide a flux estimate.
 */
class WeightedSpectrumTracker : public Tracker {
private:
  /*! @brief Minimum frequency for the spectral bins (in Hz). */
  const double _minimum_frequency;

  /*! @brief Width of a single frequency bin (in Hz). */
  const double _frequency_width;

  /*! @brief Inverse width of a single frequency bin (in Hz^-1). */
  const double _inverse_frequency_width;

  /*! @brief Side length of a cell (in m). */
  double _side_length;

  /*! @brief Number counts per bin and per photon type. */
  std::vector< double > _number_counts[PHOTONTYPE_NUMBER];

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_bins Number of bins to use.
   * @param side_length Side length of a cell (in m).
   */
  WeightedSpectrumTracker(const uint_fast32_t number_of_bins,
                          const double side_length = 0.)
      : _minimum_frequency(3.289e15),
        _frequency_width(3. * 3.289e15 / number_of_bins),
        _inverse_frequency_width(1. / _frequency_width),
        _side_length(side_length) {

    for (int_fast32_t i = 0; i < PHOTONTYPE_NUMBER; ++i) {
      _number_counts[i].resize(number_of_bins, 0.);
    }
  }

  /**
   * @brief YAMLDictionary constructor.
   *
   * The following parameters are read:
   *  - number of bins: Number of bins in the spectral histrogram
   *    (default: 100).
   *
   * @param name Name of the block in the dictionary that contains additional
   * parameters for the spectrum tracker.
   * @param blocks YAMLDictionary that contains additional parameters.
   */
  WeightedSpectrumTracker(const std::string name, YAMLDictionary &blocks)
      : WeightedSpectrumTracker(
            blocks.get_value< uint_fast32_t >(name + "number of bins", 100)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~WeightedSpectrumTracker() {}

  /**
   * @brief Normalize the tracker for use in a cell with the given size.
   *
   * @param cell Cell the tracker is attached to.
   */
  virtual void normalize_for_cell(const Cell &cell) {
    _side_length = std::cbrt(cell.get_volume());
  }

  /**
   * @brief Normalize the tracker with the appropriate ionizing luminosity per
   * unit photon packet weight.
   *
   * @param luminosity_per_weight Ionizing luminosity per unit photon packet
   * weight (in s^-1).
   */
  virtual void normalize(const double luminosity_per_weight) {
    for (int_fast32_t i = 0; i < PHOTONTYPE_NUMBER; ++i) {
      for (uint_fast32_t j = 0; j < _number_counts[i].size(); ++j) {
        _number_counts[i][j] *= luminosity_per_weight;
      }
    }
  }

  /**
   * @brief Make a duplicate of the current tracker.
   *
   * @return Pointer to a new duplicate of the tracker.
   */
  virtual Tracker *duplicate() {
    return new WeightedSpectrumTracker(_number_counts[0].size(), _side_length);
  }

  /**
   * @brief Add the contribution from the given duplicate tracker to this
   * tracker.
   *
   * @param tracker Duplicate tracker (created using Tracker::duplicate()).
   */
  virtual void merge(const Tracker *tracker) {

    const WeightedSpectrumTracker *other =
        reinterpret_cast< const WeightedSpectrumTracker * >(tracker);
    for (int_fast32_t i = 0; i < PHOTONTYPE_NUMBER; ++i) {
      for (uint_fast32_t j = 0; j < _number_counts[i].size(); ++j) {
        _number_counts[i][j] += other->_number_counts[i][j];
      }
    }
  }

  /**
   * @brief Get the projected area of a cube with side length L, as seen from
   * the given direction.
   *
   * We compute the projection by making the assumption that from the given
   * direction, exactly three faces of the cube are visible. This is an upper
   * limit, as there are special direction from which only one or two faces are
   * visible. For these cases, one or two contributions in our calculation will
   * be zero.
   *
   * The faces we consider are oriented perpendicular to the @f$x@f$, @f$y@f$
   * and @f$z@f$ axis and correspond to the top face along that axis, i.e. the
   * face with the constant high value for the corresponding coordinate. If we
   * number the vertices according to @f$vXYZ@f$, where @f$X/Y/Z=0/1@f$
   * depending on wheter the corresponding coordinate is @f$0@f$ or @f$L@f$
   * (assuming a cube anchored at the origin), the the faces we project are
   * (we assume vertices that are ordered counterclockwise when viewing the
   * face from the top):
   * @f[
   *   f_x = (v100) - (v101) - (v111) - (v110),
   * @f]
   * @f[
   *   f_y = (v101) - (v011) - (v011) - (v111),
   * @f]
   * and
   * @f[
   *   f_z = (v110) - (v111) - (v011) - (v010).
   * @f]
   *
   * We first project these faces onto the plane through the midpoint @f$m=(L/2,
   * L/2, L/2)@f$ and perpendicular to the direction @f$n@f$ using
   * @f[
   *   pXYZ = vXYZ - \left[(vXYZ-m).n\right]n,
   * @f]
   * where all quantities involved are 3D vectors and @f$.@f$ is the dot
   * product of two vectors.
   *
   * Once we have the projections, the surface area of each projected surface is
   * computed by decomposing it into two triangles. For example:
   * @f[
   *   f'_x = (p100) - (p101) - (p111) - (p110) \rightarrow{}
   *     (p100) - (p101) - (p111) + (p100) - (p111) - (p110).
   * @f]
   * The surface area of a single triangle is computed using the standard
   * formulae, e.g.
   * @f[
   *   AB = p101 - p100,
   * @f]
   * @f[
   *   AC = p111 - p100,
   * @f]
   * @f[
   *   a_{x,1} = \frac{1}{2} \sqrt{AB^2 + AC^2 - \left(AB.AC\right)^2}.
   * @f]
   *
   * @param direction Viewing direction.
   * @return Projected area (in L^2).
   */
  inline static double get_projected_area(const CoordinateVector<> direction) {

    const CoordinateVector<> v100(1., 0., 0.);
    const CoordinateVector<> v010(0., 1., 0.);
    const CoordinateVector<> v001(0., 0., 1.);
    const CoordinateVector<> v110(1., 1., 0.);
    const CoordinateVector<> v101(1., 0., 1.);
    const CoordinateVector<> v011(0., 1., 1.);
    const CoordinateVector<> v111(1., 1., 1.);
    const CoordinateVector<> m(0.5, 0.5, 0.5);

    const CoordinateVector<> p100 =
        v100 - CoordinateVector<>::dot_product(v100 - m, direction) * direction;
    const CoordinateVector<> p010 =
        v010 - CoordinateVector<>::dot_product(v010 - m, direction) * direction;
    const CoordinateVector<> p001 =
        v001 - CoordinateVector<>::dot_product(v001 - m, direction) * direction;
    const CoordinateVector<> p110 =
        v110 - CoordinateVector<>::dot_product(v110 - m, direction) * direction;
    const CoordinateVector<> p101 =
        v101 - CoordinateVector<>::dot_product(v101 - m, direction) * direction;
    const CoordinateVector<> p011 =
        v011 - CoordinateVector<>::dot_product(v011 - m, direction) * direction;
    const CoordinateVector<> p111 =
        v111 - CoordinateVector<>::dot_product(v111 - m, direction) * direction;

    const CoordinateVector<> p100_101 = p101 - p100;
    const CoordinateVector<> p100_111 = p111 - p100;
    const CoordinateVector<> p100_110 = p110 - p100;
    const double p100_101_p100_111 =
        CoordinateVector<>::dot_product(p100_101, p100_111);
    const double p100_110_p100_111 =
        CoordinateVector<>::dot_product(p100_110, p100_111);
    const double ax1 = std::sqrt(p100_101.norm2() * p100_111.norm2() -
                                 p100_101_p100_111 * p100_101_p100_111);
    const double ax2 = std::sqrt(p100_110.norm2() * p100_111.norm2() -
                                 p100_110_p100_111 * p100_110_p100_111);

    const CoordinateVector<> p110_111 = p111 - p110;
    const CoordinateVector<> p110_011 = p011 - p110;
    const CoordinateVector<> p110_010 = p010 - p110;
    const double p110_111_p110_011 =
        CoordinateVector<>::dot_product(p110_111, p110_011);
    const double p110_011_p110_010 =
        CoordinateVector<>::dot_product(p110_011, p110_010);
    const double ay1 = std::sqrt(p110_111.norm2() * p110_011.norm2() -
                                 p110_111_p110_011 * p110_111_p110_011);
    const double ay2 = std::sqrt(p110_011.norm2() * p110_010.norm2() -
                                 p110_011_p110_010 * p110_011_p110_010);

    const CoordinateVector<> p101_001 = p001 - p101;
    const CoordinateVector<> p101_011 = p011 - p101;
    const CoordinateVector<> p101_111 = p111 - p101;
    const double p101_001_p101_011 =
        CoordinateVector<>::dot_product(p101_001, p101_011);
    const double p101_011_p101_111 =
        CoordinateVector<>::dot_product(p101_011, p101_111);
    const double az1 = std::sqrt(p101_001.norm2() * p101_011.norm2() -
                                 p101_001_p101_011 * p101_001_p101_011);
    const double az2 = std::sqrt(p101_011.norm2() * p101_111.norm2() -
                                 p101_011_p101_111 * p101_011_p101_111);

    return 0.5 * (ax1 + ax2 + ay1 + ay2 + az1 + az2);
  }

  /**
   * @brief Add the contribution of the given photon packet to the bins.
   *
   * @param photon Photon to add.
   */
  virtual void count_photon(const Photon &photon) {

    cmac_assert_message(_side_length > 0., "Tracker was not normalized!");

    const double frequency = photon.get_energy();
    const uint_fast32_t index =
        (frequency - _minimum_frequency) * _inverse_frequency_width;
    if (index < _number_counts[0].size()) {
      const double weight = get_projected_area(photon.get_direction()) *
                            _side_length * _side_length;
      _number_counts[photon.get_type()][index] += 1. / weight;
    }
  }

  /**
   * @brief Add the contribution of the given photon packet to the bins.
   *
   * @param photon Photon to add.
   */
  virtual void count_photon(const PhotonPacket &photon) {

    cmac_assert_message(_side_length > 0., "Tracker was not normalized!");

    const double frequency = photon.get_energy();
    const uint_fast32_t index =
        (frequency - _minimum_frequency) * _inverse_frequency_width;
    if (index < _number_counts[0].size()) {
      const double weight = get_projected_area(photon.get_direction()) *
                            _side_length * _side_length;
      _number_counts[photon.get_type()][index] += 1. / weight;
    }
  }

  /**
   * @brief Output the spectrum to the file with the given name.
   *
   * @param filename Name of the output file.
   */
  virtual void output_tracker(const std::string filename) const {

    std::ofstream ofile(filename);
    ofile << "# frequency (Hz)";
    for (int_fast32_t i = 0; i < PHOTONTYPE_NUMBER; ++i) {
      ofile << "\t" << get_photontype_name(i) << " flux (s^-1 m^-2)";
    }
    ofile << "\n";
    for (uint_fast32_t i = 0; i < _number_counts[0].size(); ++i) {
      const double nu = _minimum_frequency + (i + 0.5) * _frequency_width;
      ofile << nu;
      for (int_fast32_t j = 0; j < PHOTONTYPE_NUMBER; ++j) {
        ofile << "\t" << _number_counts[j][i];
      }
      ofile << "\n";
    }
  }

#ifdef HAVE_HDF5
  /**
   * @brief Output the tracker data to the given HDF5 group with the given name.
   *
   * @param group HDF5Group to write to.
   */
  virtual void output_tracker_to_hdf5(const HDF5Tools::HDF5Group group) {
    std::string type_string = "WeightedSpectrum";
    HDF5Tools::write_attribute< std::string >(group, "type", type_string);
    std::string unit_string = "s^-1";
    HDF5Tools::write_attribute< std::string >(group, "frequency unit",
                                              unit_string);
    unit_string = "m^-2 s^-1";
    HDF5Tools::write_attribute< std::string >(group, "flux unit", unit_string);
    std::vector< double > frequencies(_number_counts[0].size(), 0.);
    for (uint_fast32_t i = 0; i < frequencies.size(); ++i) {
      frequencies[i] = _minimum_frequency + (i + 0.5) * _frequency_width;
    }
    HDF5Tools::write_dataset(group, "frequencies", frequencies);
    for (int_fast32_t i = 0; i < PHOTONTYPE_NUMBER; ++i) {
      std::stringstream namestr;
      namestr << get_photontype_name(i) << " flux";
      HDF5Tools::write_dataset(group, namestr.str(), _number_counts[i]);
    }
  }
#endif

  /**
   * @brief Describe the tracker in the given output stream, appending the given
   * prefix to each line of output.
   *
   * @param prefix Prefix to add to each output line.
   * @param stream std::ostream to write to.
   */
  virtual void describe(const std::string prefix, std::ostream &stream) const {

    stream << prefix << "type: WeightedSpectrum\n";
    stream << prefix << "number of bins: " << _number_counts[0].size() << "\n";
  }
};

#endif // WEIGHTEDSPECTRUMTRACKER_HPP
