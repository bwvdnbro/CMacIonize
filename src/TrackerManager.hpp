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
 * @file TrackerManager.hpp
 *
 * @brief Class that manages Trackers in the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TRACKERMANAGER_HPP
#define TRACKERMANAGER_HPP

#include "Configuration.hpp"
#include "DensityGrid.hpp"
#include "DensitySubGridCreator.hpp"
#include "ParameterFile.hpp"
#include "TrackerFactory.hpp"
#include "YAMLDictionary.hpp"

#ifdef HAVE_HDF5
#include "HDF5Tools.hpp"
#endif

#include <fstream>
#include <vector>

/**
 * @brief Class that manages Trackers in the grid.
 */
class TrackerManager {
private:
  /*! @brief List of tracker positions. */
  std::vector< CoordinateVector<> > _tracker_positions;

  /*! @brief List of trackers. */
  std::vector< Tracker * > _trackers;

  /*! @brief List of multi-trackers (for cells with more than one tracker). */
  std::vector< MultiTracker * > _multi_trackers;

  /*! @brief Indices of the original tracker corresponding to each copy. */
  std::vector< size_t > _originals;

  /*! @brief Indices of the first copy of each tracker. */
  std::vector< size_t > _copies;

  /*! @brief List of output file names. */
  std::vector< std::string > _output_names;

  /*! @brief Number of photon packets to use during the tracking step. */
  const uint_fast64_t _number_of_photons;

  /*! @brief Output all trackers to a single HDF5 file? */
  const bool _hdf5_output;

  /*! @brief Name of the HDF5 output file. */
  const std::string _hdf5_name;

#ifdef HAVE_HDF5
  /*! @brief Index of the first tracker in each tracker group. */
  std::vector< size_t > _tracker_groups;

  /*! @brief Group size for each tracker group. */
  std::vector< size_t > _group_size;

  /*! @brief Tracker group a each tracker belongs to. */
  std::vector< size_t > _group_index;
#endif

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the file that contains the positions of the
   * trackers.
   * @param number_of_photons Number of photon packets to use during the
   * tracking step.
   * @param hdf5_output Output all trackers to a single HDF5 file?
   * @param hdf5_name Name of the HDF5 output file.
   */
  TrackerManager(const std::string filename,
                 const uint_fast64_t number_of_photons,
                 const bool hdf5_output = false,
                 const std::string hdf5_name = "")
      : _number_of_photons(number_of_photons), _hdf5_output(hdf5_output),
        _hdf5_name(hdf5_name) {

#ifndef HAVE_HDF5
    if (hdf5_output) {
      cmac_error("HDF5 output requested, but HDF5 not enabled!");
    }
#endif

    std::ifstream file(filename);

    if (!file) {
      cmac_error("Error while opening file \"%s\"!", filename.c_str());
    }

    YAMLDictionary blocks(file);
    const uint_fast32_t number_of_trackers =
        blocks.get_value< uint_fast32_t >("number of trackers");
    _tracker_positions.resize(number_of_trackers);
    _trackers.resize(number_of_trackers, nullptr);
    _copies.resize(number_of_trackers, 0xffffffff);
    _output_names.resize(number_of_trackers);
    for (uint_fast32_t i = 0; i < number_of_trackers; ++i) {
      std::stringstream blockname;
      blockname << "tracker[" << i << "]:";

      _tracker_positions[i] = blocks.get_physical_vector< QUANTITY_LENGTH >(
          blockname.str() + "position");

      _trackers[i] = TrackerFactory::generate(blockname.str(), blocks);

      std::stringstream default_name;
      default_name << "Tracker" << i;
      if (!_hdf5_output) {
        default_name << ".txt";
      }
      _output_names[i] = blocks.get_value< std::string >(
          blockname.str() + "output name", default_name.str());
    }

#ifdef HAVE_HDF5
    if (_hdf5_output) {
      _group_index.resize(number_of_trackers, 0);
      // we skip the first tracker, it trivially belongs to group 0
      _tracker_groups.push_back(0);
      _group_size.push_back(1);
      for (uint_fast32_t i = 1; i < number_of_trackers; ++i) {
        uint_fast32_t group_id = 0;
        while (
            group_id < _tracker_groups.size() &&
            !_trackers[_tracker_groups[group_id]]->same_group(_trackers[i])) {
          ++group_id;
        }
        if (group_id == _tracker_groups.size()) {
          _tracker_groups.push_back(i);
          _group_size.push_back(0);
        }
        _group_index[i] = group_id;
        ++_group_size[group_id];
      }
    }
#endif

    std::ofstream ofile(filename + ".used-values");
    blocks.print_contents(ofile, true);
    ofile.close();
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read from the parameter file:
   *  - filename: Name of the file that contains the tracker positions
   *    (required)
   *  - minimum number of photon packets: Minimum number of photon packets to
   *    use during the spectrum tracking step (default: 0, meaning we do not
   *    use a different number for the spectrum tracking step)
   *  - HDF5 output: Output all trackers to a single HDF5 file? (default: false)
   *  - HDF5 output name: Name of the HDF5 output file (default: trackers.hdf5).
   *
   * @param params ParameterFile to read from.
   */
  TrackerManager(ParameterFile &params)
      : TrackerManager(
            params.get_filename("TrackerManager:filename"),
            params.get_value< uint_fast64_t >(
                "TrackerManager:minimum number of photon packets", 0),
            params.get_value< bool >("TrackerManager:HDF5 output", false),
            params.get_value< std::string >("TrackerManager:HDF5 output name",
                                            "trackers.hdf5")) {}

  /**
   * @brief Destructor.
   */
  ~TrackerManager() {
    for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
      delete _trackers[i];
    }
    for (uint_fast32_t i = 0; i < _multi_trackers.size(); ++i) {
      _multi_trackers[i]->erase_trackers();
      delete _multi_trackers[i];
    }
  }

  /**
   * @brief Add the trackers to the given DensityGrid.
   *
   * @param grid Grid to add trackers to.
   */
  inline void add_trackers(DensityGrid &grid) {
    for (uint_fast32_t i = 0; i < _tracker_positions.size(); ++i) {
      if (!grid.get_box().inside(_tracker_positions[i])) {
        cmac_error("Tracker is not inside grid!");
      }
      auto it = grid.get_cell(_tracker_positions[i]);
      _trackers[i]->normalize_for_cell(it);
      IonizationVariables &ionization_variables = it.get_ionization_variables();
      if (ionization_variables.get_tracker() != nullptr) {
        MultiTracker *multi_tracker;
        if (ionization_variables.get_tracker()->is_multi_tracker()) {
          multi_tracker =
              static_cast< MultiTracker * >(ionization_variables.get_tracker());
        } else {
          multi_tracker = new MultiTracker();
          multi_tracker->add_tracker(ionization_variables.get_tracker());
          ionization_variables.add_tracker(multi_tracker);
          _multi_trackers.push_back(multi_tracker);
        }
        multi_tracker->add_tracker(_trackers[i]);
      } else {
        ionization_variables.add_tracker(_trackers[i]);
      }
    }
  }

  /**
   * @brief Add the trackers to the given distributed grid.
   *
   * @param grid Grid to add trackers to.
   */
  template < class _subgrid_type_ >
  inline void add_trackers(DensitySubGridCreator< _subgrid_type_ > &grid) {
    for (uint_fast32_t i = 0; i < _tracker_positions.size(); ++i) {
      if (!grid.get_box().inside(_tracker_positions[i])) {
        cmac_error("Tracker is not inside grid!");
      }
      auto gridit = grid.get_subgrid(_tracker_positions[i]);
      {
        auto cellit = (*gridit).get_cell(_tracker_positions[i]);
        _trackers[i]->normalize_for_cell(cellit);
        IonizationVariables &ionization_variables =
            cellit.get_ionization_variables();
        if (ionization_variables.get_tracker() != nullptr) {
          MultiTracker *multi_tracker;
          if (ionization_variables.get_tracker()->is_multi_tracker()) {
            multi_tracker = static_cast< MultiTracker * >(
                ionization_variables.get_tracker());
          } else {
            multi_tracker = new MultiTracker();
            multi_tracker->add_tracker(ionization_variables.get_tracker());
            ionization_variables.add_tracker(multi_tracker);
            _multi_trackers.push_back(multi_tracker);
          }
          multi_tracker->add_tracker(_trackers[i]);
        } else {
          ionization_variables.add_tracker(_trackers[i]);
        }
      }
      auto copies = gridit.get_copies();
      bool first = true;
      for (auto copyit = copies.first; copyit != copies.second; ++copyit) {
        Tracker *new_tracker = _trackers[i]->duplicate();
        _trackers.push_back(new_tracker);
        _originals.push_back(i);
        if (first) {
          _copies[i] = _trackers.size() - 1;
          first = false;
        }
        IonizationVariables &ionization_variables =
            (*copyit)
                .get_cell(_tracker_positions[i])
                .get_ionization_variables();
        if (ionization_variables.get_tracker() != nullptr) {
          MultiTracker *multi_tracker;
          if (ionization_variables.get_tracker()->is_multi_tracker()) {
            multi_tracker = static_cast< MultiTracker * >(
                ionization_variables.get_tracker());
          } else {
            multi_tracker = new MultiTracker();
            multi_tracker->add_tracker(ionization_variables.get_tracker());
            ionization_variables.add_tracker(multi_tracker);
            _multi_trackers.push_back(multi_tracker);
          }
          multi_tracker->add_tracker(_trackers[i]);
        } else {
          ionization_variables.add_tracker(_trackers[i]);
        }
      }
    }
  }

  /**
   * @brief Normalize the trackers based on the actual physical weight of each
   * photon packet.
   *
   * This function also merges copies made by add_trackers().
   *
   * @param luminosity_per_weight Ionizing luminosity per unit photon packet
   * weight used to normalize photon packet contributions (in s^-1).
   */
  inline void normalize(const double luminosity_per_weight) {
    for (uint_fast32_t i = 0; i < _tracker_positions.size(); ++i) {
      if (_copies[i] != 0xffffffff) {
        size_t copy = _copies[i] - _tracker_positions.size();
        while (copy < _originals.size() && _originals[copy] == i) {
          _trackers[i]->merge(_trackers[copy + _tracker_positions.size()]);
          ++copy;
        }
      }
      _trackers[i]->normalize(luminosity_per_weight);
    }
  }

  /**
   * @brief Output the tracker information.
   */
  inline void output_trackers() const {
    if (_hdf5_output) {
#ifdef HAVE_HDF5
      HDF5Tools::HDF5File file =
          HDF5Tools::open_file(_hdf5_name, HDF5Tools::HDF5FILEMODE_WRITE);
      for (uint_fast32_t igroup = 0; igroup < _tracker_groups.size();
           ++igroup) {
        std::stringstream groupname;
        groupname << "Group" << igroup;
        HDF5Tools::HDF5Group group =
            HDF5Tools::create_group(file, groupname.str());
        std::string unit_string = "m";
        const uint_fast32_t group_size = _group_size[igroup];
        _trackers[_tracker_groups[igroup]]->create_group(group, group_size);
        std::vector< CoordinateVector<> > positions(group_size);
        std::vector< std::string > labels(group_size);
        uint_fast32_t group_id = 0;
        for (uint_fast32_t i = 0; i < _trackers.size(); ++i) {
          if (_group_index[i] == igroup) {
            positions[group_id] = _tracker_positions[i];
            labels[group_id] = _output_names[i];
            _trackers[i]->append_to_group(group, group_id);
            ++group_id;
            if (group_id == group_size) {
              break;
            }
          }
        }
        HDF5Tools::write_dataset(group, "positions", positions);
        HDF5Tools::write_dataset(group, "tracker labels", labels);
        HDF5Tools::write_attribute< std::string >(group, "position unit",
                                                  unit_string);
        HDF5Tools::close_group(group);
      }
      HDF5Tools::close_file(file);
#endif
    } else {
      for (uint_fast32_t i = 0; i < _tracker_positions.size(); ++i) {
        _trackers[i]->output_tracker(_output_names[i]);
      }
    }
  }

  /**
   * @brief Get the number of photon packets to use during the tracking step.
   *
   * @return Number of photons to use during the tracking step.
   */
  inline uint_fast64_t get_number_of_photons() const {
    return _number_of_photons;
  }
};

#endif // TRACKERMANAGER_HPP
