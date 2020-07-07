/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *                    Nina Sartorio (sartorio.nina@gmail.com)
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
 * @file AbsorptionTracker.hpp
 *
 * @brief Tracker that tracks the absorption within a cell.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Nina Sartorio (sartorio.nina@gmail.com)
 */
#ifndef ABSORPTIONTRACKER_HPP
#define ABSORPTIONTRACKER_HPP

#include "Configuration.hpp"
#include "ElementNames.hpp"
#include "Error.hpp"
#include "PhotonPacket.hpp"
#include "PhotonType.hpp"
#include "Tracker.hpp"
#include "YAMLDictionary.hpp"

#ifdef HAVE_HDF5
#include "HDF5Tools.hpp"
#endif

#include <fstream>
#include <typeinfo>

/**
 * @brief General interface for trackers that record photon properties.
 */
class AbsorptionTracker : public Tracker {
private:
  /*! @brief the bins for absorption. One for each photon type (source, H
   *  reemission, He reemission) and each ion. Units of m^3. */
  double _absorption_bins[PHOTONTYPE_NUMBER][NUMBER_OF_IONNAMES];

public:
  /**
   * @brief Constructor.
   */
  inline AbsorptionTracker() {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; i++) {
      for (int_fast32_t j = 0; j < PHOTONTYPE_NUMBER; j++) {
        _absorption_bins[j][i] = 0.;
      }
    }
  }

  /**
   * @brief YAMLDictionary constructor.
   *
   * @param name Name of the block in the dictionary that contains additional
   * parameters for the spectrum tracker.
   * @param blocks YAMLDictionary that contains additional parameters.
   */
  AbsorptionTracker(const std::string name, YAMLDictionary &blocks)
      : AbsorptionTracker() {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~AbsorptionTracker() {}

  /**
   * @brief Normalize the tracker with the appropriate ionizing luminosity per
   * unit photon packet weight.
   *
   * @param luminosity_per_weight Ionizing luminosity per unit photon packet
   * weight (in s^-1).
   */
  virtual void normalize(const double luminosity_per_weight) {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; i++) {
      for (int_fast32_t j = 0; j < PHOTONTYPE_NUMBER; j++) {
        _absorption_bins[j][i] *= luminosity_per_weight;
      }
    }
  }

  /**
   * @brief Make a duplicate of the current tracker.
   *
   * @return Pointer to a new duplicate of the tracker.
   */
  virtual Tracker *duplicate() { return new AbsorptionTracker(); }

  /**
   * @brief Add the contribution from the given duplicate tracker to this
   * tracker.
   *
   * @param tracker Duplicate tracker (created using Tracker::duplicate()).
   */
  virtual void merge(const Tracker *tracker) {
    const AbsorptionTracker *other_tracker =
        static_cast< const AbsorptionTracker * >(tracker);
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; i++) {
      for (int_fast32_t j = 0; j < PHOTONTYPE_NUMBER; j++) {
        _absorption_bins[j][i] += other_tracker->_absorption_bins[j][i];
      }
    }
  }

  /**
   * @brief Add the contribution of the given photon packet to the tracker.
   *
   * @param photon Photon to add.
   */
  virtual void count_photon(const Photon &photon) {
    cmac_error("Function not in use in current version of code")
  }

  /**
   * @brief Add the contribution of the given photon packet to the tracker.
   *
   * Note this is the absorption volume. In order to get the number of photons
   * absorbed one has to multiply by the number density of hydrogen (as cross
   * sections are abundance weighted).
   *
   * @param photon Photon to add.
   * @param absorption Absorption counters within the cell for this photon
   * (in m^3).
   */
  virtual void count_photon(const PhotonPacket &photon,
                            const double *absorption) {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; i++) {
      _absorption_bins[i][photon.get_type()] += absorption[i];
    }
  }

  /**
   * @brief Output the tracker data to the file with the given name.
   *
   * @param filename Name of the output file.
   */
  virtual void output_tracker(const std::string filename) const {
    std::ofstream ofile(filename);
    ofile << "# Ion ";
    for (int_fast32_t j = 0; j < PHOTONTYPE_NUMBER; j++) {
      ofile << "\t" << get_photontype_name(j);
    }
    ofile << "\n";
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; i++) {
      ofile << get_ion_name(i);
      for (int_fast32_t j = 0; j < PHOTONTYPE_NUMBER; j++) {
        ofile << "\t" << _absorption_bins[j][i];
      }
      ofile << "\n";
    }
  }

#ifdef HAVE_HDF5
  /**
   * @brief Does the given tracker belong to the same group as this tracker?
   *
   * @param tracker Other tracker.
   * @return True if both trackers belong to the same group.
   */
  virtual bool same_group(const Tracker *tracker) const {
    return typeid(*tracker).hash_code() == typeid(*this).hash_code();
  }

  /**
   * @brief Create the header and shared datasets for an HDF5 group containing
   * one or multiple trackers of this type.
   *
   * @param group HDF5Group to write to.
   * @param group_size Number of trackers in the group.
   */
  virtual void create_group(const HDF5Tools::HDF5Group group,
                            const uint_fast32_t group_size) {
    std::vector< std::string > ion_names(NUMBER_OF_IONNAMES);
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; i++) {
      ion_names[i] = get_ion_name(i);
    }
    HDF5Tools::write_dataset(group, "ion name", ion_names);

    for (int_fast32_t i = 0; i < PHOTONTYPE_NUMBER; i++) {
      std::stringstream mystring;
      mystring << get_photontype_name(i) << " absorption";
      HDF5Tools::create_datatable< double >(group, mystring.str(), group_size,
                                            NUMBER_OF_IONNAMES);
    }
  }

  /**
   * @brief Append the tracker to the given group.
   *
   * We assume the group was already properly initialised using create_group().
   *
   * @param group HDF5Group to write to.
   * @param group_index Index of this particular tracker within the group.
   */
  virtual void append_to_group(const HDF5Tools::HDF5Group group,
                               const uint_fast32_t group_index) {
    for (int_fast32_t i = 0; i < PHOTONTYPE_NUMBER; i++) {
      std::stringstream mystring;
      mystring << get_photontype_name(i) << " absorption";
      HDF5Tools::fill_row(group, mystring.str(), group_index,
                          _absorption_bins[i]);
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
    stream << prefix << "type: AbsorptionTracker\n";
  }
};

#endif // ABSORPTIONTRACKER_HPP
