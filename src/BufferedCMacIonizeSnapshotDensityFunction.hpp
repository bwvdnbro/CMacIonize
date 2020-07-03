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
 * @file BufferedCMacIonizeSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads a density grid from a task-based CMacIonize
 * snapshot in a buffered fashion.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef BUFFEREDCMACIONIZESNAPSHOTDENSITYFUNCTION_HPP
#define BUFFEREDCMACIONIZESNAPSHOTDENSITYFUNCTION_HPP

#include "Box.hpp"
#include "CPUCycle.hpp"
#include "DensityFunction.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "ThreadLock.hpp"

/**
 * @brief DensityFunction that reads a density grid from a task-based CMacIonize
 * snapshot in a buffered fashion.
 *
 * The snapshot is read per block, with blocks being kept in internal buffers
 * as long as there is space. If buffer space runs out, the oldest block is
 * discarded. This should allow reading snapshots that are much bigger than the
 * memory of the machine running the code.
 *
 * On top of the improved memory usage, this variant of
 * CMacIonizeSnapshotDensityFunction also ensures mass conservation when
 * degrading the grid resolution.
 */
class BufferedCMacIonizeSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief Number of subgrids that can be buffered. */
  const uint_fast32_t _buffer_size;

  /*! @brief Snapshot file. */
  HDF5Tools::HDF5File _file;

  /*! @brief Particle group within the snapshot file. */
  HDF5Tools::HDF5Group _particle_group;

  /*! @brief Number of old cells per new cell (in 1 coordinate dimension). If
   *  larger than one, a mass-conserving mapping is used during read. */
  uint_fast32_t _number_of_old_cells_per_new_cell_1D;

  /*! @brief Anchor of the old simulation box (in m). */
  CoordinateVector<> _old_anchor;

  /*! @brief Width of a subgrid in the old simulation (in m). */
  CoordinateVector<> _subgrid_width;

  /*! @brief Number of cells in 1D in an original subgrid. */
  CoordinateVector< uint_fast32_t > _original_subgrid_ncell;

  /*! @brief Number of cells in an original subgrid. */
  uint_fast32_t _original_subgrid_size;

  /*! @brief Number of cells in 1D in a remapped subgrid. */
  CoordinateVector< uint_fast32_t > _mapped_subgrid_ncell;

  /*! @brief Number of cells in a remapped subgrid. */
  uint_fast32_t _mapped_subgrid_size;

  /*! @brief Number of subgrids in the original snapshot in each dimension. */
  CoordinateVector< uint_fast32_t > _number_of_subgrids;

  /*! @brief Subgrid buffer. */
  std::vector< DensityValues > _buffer;

  /*! @brief Last usage timestamp for each buffer. */
  std::vector< uint_fast64_t > _buffer_timestamps;

  /*! @brief Index of each subgrid in the buffer (if buffered). */
  std::vector< uint_fast32_t > _buffer_indices;

  /*! @brief Index of the subgrid that has been buffered (to make sure another
   *  thread did not swap it out). */
  std::vector< uint_fast32_t > _buffer_subgrid_indices;

  /*! @brief Lock protecting the buffer and HDF5 file. */
  ThreadLock _buffer_lock;

  /*! @brief Locks per buffer element. */
  std::vector< ThreadLock > _buffer_element_locks;

  /*! @brief Read the number density (true) or density (false)? */
  bool _read_number_density;

  /*! @brief Read the temperature (true) or pressure (false)? */
  bool _read_temperature;

  /*! @brief Flags indicating which ionic fractions to read. */
  bool _read_ionic_fraction[NUMBER_OF_IONNAMES];

  /*! @brief Read the velocity? */
  bool _read_velocity;

  /*! @brief Log to write logging info to. */
  Log *_log;

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the snapshot file.
   * @param buffer_size Number of subgrids that can be buffered.
   * @param new_box Simulation box of the current simulation (in m).
   * @param new_ncell Number of cells in the current simulation.
   * @param log Log to write logging info to.
   */
  BufferedCMacIonizeSnapshotDensityFunction(
      const std::string filename, const uint_fast32_t buffer_size,
      const Box<> new_box, const CoordinateVector< uint_fast32_t > new_ncell,
      Log *log = nullptr)
      : _buffer_size(buffer_size), _buffer_timestamps(buffer_size, 0),
        _buffer_subgrid_indices(buffer_size),
        _buffer_element_locks(buffer_size), _log(log) {

    // check that the file can be opened
    std::ifstream file(filename);
    if (!file.is_open()) {
      cmac_error("Could not open file \"%s\"!", filename.c_str());
    }

    // open the file using HDF5
    _file = HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

    // parse the parameter block
    HDF5Tools::HDF5Group group = HDF5Tools::open_group(_file, "/Parameters");
    std::vector< std::string > parameternames =
        HDF5Tools::get_attribute_names(group);
    ParameterFile parameters;
    for (auto it = parameternames.begin(); it != parameternames.end(); ++it) {
      std::string attname = *it;
      std::string attvalue =
          HDF5Tools::read_attribute< std::string >(group, attname);
      parameters.add_value(attname, attvalue);
    }
    HDF5Tools::close_group(group);

    // check that the file is a task-based snapshot
    if (!parameters.has_value("DensitySubGridCreator:number of subgrids")) {
      cmac_error("A BufferedCMacIonizeSnapshotDensityFunction can only be used "
                 "to read task-based CMacIonize snapshots!");
    }

    // now get the parameters for the original task-based grid
    _number_of_subgrids =
        parameters.get_value< CoordinateVector< uint_fast32_t > >(
            "DensitySubGridCreator:number of subgrids");
    Box<> old_box(parameters.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:anchor"),
                  parameters.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:sides"));
    CoordinateVector< uint_fast32_t > old_ncell =
        parameters.get_value< CoordinateVector< uint_fast32_t > >(
            "DensityGrid:number of cells");

    // check that the new simulation box is contained within the old one
    const CoordinateVector<> anchor_in_old_box =
        new_box.get_anchor() - old_box.get_anchor();
    const CoordinateVector<> available_sides =
        new_box.get_top_anchor() - old_box.get_anchor();
    if (log) {
      log->write_info("Anchor in old box: [", anchor_in_old_box.x(), " m, ",
                      anchor_in_old_box.y(), " m, ", anchor_in_old_box.z(),
                      " m].");
      log->write_info("Available sides: [", available_sides.x(), " m, ",
                      available_sides.y(), " m, ", available_sides.z(), " m].");
    }
    if (anchor_in_old_box.x() < 0. || anchor_in_old_box.y() < 0. ||
        anchor_in_old_box.z() < 0. ||
        available_sides.x() < new_box.get_sides().x() ||
        available_sides.y() < new_box.get_sides().y() ||
        available_sides.z() < new_box.get_sides().z()) {
      cmac_error("New simulation box is not inside old simulation box!");
    }

    // check that the new box is compatible with the old resolution
    const CoordinateVector<> old_cell_size(
        old_box.get_sides().x() / old_ncell.x(),
        old_box.get_sides().y() / old_ncell.y(),
        old_box.get_sides().z() / old_ncell.z());
    const CoordinateVector<> cell_offset_float(
        anchor_in_old_box.x() / old_cell_size.x(),
        anchor_in_old_box.y() / old_cell_size.y(),
        anchor_in_old_box.z() / old_cell_size.z());
    const CoordinateVector< uint_fast32_t > cell_offset(
        static_cast< uint_fast32_t >(std::round(cell_offset_float.x())),
        static_cast< uint_fast32_t >(std::round(cell_offset_float.y())),
        static_cast< uint_fast32_t >(std::round(cell_offset_float.z())));
    const CoordinateVector<> cell_sides_float = {
        available_sides.x() / old_cell_size.x(),
        available_sides.y() / old_cell_size.y(),
        available_sides.z() / old_cell_size.z()};
    const CoordinateVector< uint_fast32_t > cell_sides = {
        static_cast< uint_fast32_t >(std::round(cell_sides_float.x())),
        static_cast< uint_fast32_t >(std::round(cell_sides_float.y())),
        static_cast< uint_fast32_t >(std::round(cell_sides_float.z()))};
    if (log) {
      log->write_info("Old cell size: [", old_cell_size.x(), " m, ",
                      old_cell_size.y(), " m, ", old_cell_size.z(), " m].");
      log->write_info("Cell offset float: [", cell_offset_float.x(), ", ",
                      cell_offset_float.y(), ", ", cell_offset_float.z(), "].");
      log->write_info("Cell offset: [", cell_offset.x(), ", ", cell_offset.y(),
                      ", ", cell_offset.z(), "].");
      log->write_info("Cell sides float: [", cell_sides_float.x(), ", ",
                      cell_sides_float.y(), ", ", cell_sides_float.z(), "].");
      log->write_info("Cell sides: [", cell_sides.x(), ", ", cell_sides.y(),
                      ", ", cell_sides.z(), "].");
    }
    if ((std::abs(cell_offset_float.x()) - cell_offset.x()) > 1.e-10 ||
        (std::abs(cell_offset_float.y()) - cell_offset.y()) > 1.e-10 ||
        (std::abs(cell_offset_float.z()) - cell_offset.z()) > 1.e-10 ||
        (std::abs(cell_sides_float.x()) - cell_sides.x()) > 1.e-10 ||
        (std::abs(cell_sides_float.y()) - cell_sides.y()) > 1.e-10 ||
        (std::abs(cell_sides_float.z()) - cell_sides.z()) > 1.e-10) {
      cmac_error("New box not compatible with old resolution!");
    }

    // check that the old and new grid resolution are compatible
    // we currently limit ourselves to square cells and boxes
    const CoordinateVector<> new_cell_size(
        new_box.get_sides().x() / new_ncell.x(),
        new_box.get_sides().y() / new_ncell.y(),
        new_box.get_sides().z() / new_ncell.z());
    if (log) {
      log->write_info("New cell size: [", new_cell_size.x(), " m, ",
                      new_cell_size.y(), " m, ", new_cell_size.z(), " m].");
    }
    if (std::abs(old_box.get_sides().x() - old_box.get_sides().y()) > 1.e-10 ||
        std::abs(old_box.get_sides().x() - old_box.get_sides().z()) > 1.e-10 ||
        std::abs(new_box.get_sides().x() - new_box.get_sides().y()) > 1.e-10 ||
        std::abs(new_box.get_sides().x() - new_box.get_sides().z()) > 1.e-10 ||
        std::abs(old_cell_size.x() - old_cell_size.y()) > 1.e-10 ||
        std::abs(old_cell_size.x() - old_cell_size.z()) > 1.e-10 ||
        std::abs(new_cell_size.x() - new_cell_size.y()) > 1.e-10 ||
        std::abs(new_cell_size.x() - new_cell_size.z()) > 1.e-10) {
      cmac_error("Buffered snapshot reading currently only works for square "
                 "boxes and cells!");
    }
    if (cell_sides.x() <= new_ncell.x()) {
      // there is a many/one to one mapping between new and old cells, and we
      // don't need to do anything special when reading cells
      _number_of_old_cells_per_new_cell_1D = 1;
    } else {
      if (cell_sides.x() % new_ncell.x() != 0) {
        cmac_error("New resolution not compatible with old resolution!");
      }
      _number_of_old_cells_per_new_cell_1D = cell_sides.x() / new_ncell.x();
    }

    _original_subgrid_ncell = CoordinateVector< uint_fast32_t >(
        old_ncell.x() / _number_of_subgrids.x(),
        old_ncell.y() / _number_of_subgrids.y(),
        old_ncell.z() / _number_of_subgrids.z());
    _subgrid_width =
        CoordinateVector<>(old_box.get_sides().x() / _number_of_subgrids.x(),
                           old_box.get_sides().y() / _number_of_subgrids.y(),
                           old_box.get_sides().z() / _number_of_subgrids.z());

    if (log) {
      log->write_info("Number of old cells per new cell 1D: ",
                      _number_of_old_cells_per_new_cell_1D);
      log->write_info("Original subgrid ncell: [", _original_subgrid_ncell.x(),
                      ", ", _original_subgrid_ncell.y(), ", ",
                      _original_subgrid_ncell.z(), "].");
      log->write_info("Subgrid width: [", _subgrid_width.x(), " m, ",
                      _subgrid_width.y(), " m, ", _subgrid_width.z(), " m].");
    }

    // we also make sure that we don't need to degrade across subgrid boundaries
    if (_number_of_old_cells_per_new_cell_1D > 1) {
      const CoordinateVector<> subgrid_offset_float(
          anchor_in_old_box.x() / _subgrid_width.x(),
          anchor_in_old_box.y() / _subgrid_width.y(),
          anchor_in_old_box.z() / _subgrid_width.z());
      const CoordinateVector< uint_fast32_t > subgrid_offset(
          static_cast< uint_fast32_t >(std::round(subgrid_offset_float.x())),
          static_cast< uint_fast32_t >(std::round(subgrid_offset_float.y())),
          static_cast< uint_fast32_t >(std::round(subgrid_offset_float.z())));
      if (log) {
        log->write_info("Subgrid offset float: [", subgrid_offset_float.x(),
                        ", ", subgrid_offset_float.y(), ", ",
                        subgrid_offset_float.z(), "].");
        log->write_info("Subgrid offset: [", subgrid_offset.x(), ", ",
                        subgrid_offset.y(), ", ", subgrid_offset.z(), "].");
      }
      if (std::abs(subgrid_offset_float.x() - subgrid_offset.x()) > 1.e-10 ||
          std::abs(subgrid_offset_float.y() - subgrid_offset.y()) > 1.e-10 ||
          std::abs(subgrid_offset_float.z() - subgrid_offset.z()) > 1.e-10 ||
          _original_subgrid_ncell.x() % _number_of_old_cells_per_new_cell_1D !=
              0) {
        cmac_error("Degrading resolution across subgrid boundaries not yet "
                   "supported!");
      }
    }

    // now compute the relevant quantities for locating blocks
    _old_anchor = old_box.get_anchor();

    _original_subgrid_size = _original_subgrid_ncell.x() *
                             _original_subgrid_ncell.y() *
                             _original_subgrid_ncell.z();
    _mapped_subgrid_ncell =
        _original_subgrid_ncell / _number_of_old_cells_per_new_cell_1D;
    _mapped_subgrid_size = _mapped_subgrid_ncell.x() *
                           _mapped_subgrid_ncell.y() *
                           _mapped_subgrid_ncell.z();

    // open the particle group
    _particle_group = HDF5Tools::open_group(_file, "PartType0");

    // figure out which values to read
    if (HDF5Tools::group_exists(_particle_group, "NumberDensity")) {
      _read_number_density = true;
    } else {
      if (!HDF5Tools::group_exists(_particle_group, "Density")) {
        cmac_error("No density variable present in snapshot file!");
      }
      _read_number_density = false;
    }
    if (HDF5Tools::group_exists(_particle_group, "Temperature")) {
      _read_temperature = true;
    } else {
      if (!HDF5Tools::group_exists(_particle_group, "Pressure")) {
        cmac_error("No temperature variable present in snapshot file!");
      }
      _read_temperature = false;
    }
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _read_ionic_fraction[i] = HDF5Tools::group_exists(
          _particle_group, "NeutralFraction" + get_ion_name(i));
    }
    _read_velocity = HDF5Tools::group_exists(_particle_group, "Velocities");

    if (log) {
      log->write_info("Old anchor: [", _old_anchor.x(), " m, ", _old_anchor.y(),
                      " m, ", _old_anchor.z(), " m].");
      log->write_info("Original subgrid size: ", _original_subgrid_size);
      log->write_info("Mapped subgrid ncell: [", _mapped_subgrid_ncell.x(),
                      ", ", _mapped_subgrid_ncell.y(), ", ",
                      _mapped_subgrid_ncell.z(), "].");
      log->write_info("Mapped subgrid size: ", _mapped_subgrid_size);
      log->write_info("Read number density: ", _read_number_density);
      log->write_info("Read temperature: ", _read_temperature);
      for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        log->write_info("Read ionic fraction ", get_ion_name(i), ": ",
                        _read_ionic_fraction[i]);
      }
      log->write_info("Read velocity: ", _read_velocity);
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read:
   *  - filename: Name of the snapshot file containing the initial condition
   *    (required).
   *  - buffer size: Number of subgrids that can be stored in the internal
   *    buffer (default: 100).
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  BufferedCMacIonizeSnapshotDensityFunction(ParameterFile &params,
                                            Log *log = nullptr)
      : BufferedCMacIonizeSnapshotDensityFunction(
            params.get_filename("DensityFunction:filename"),
            params.get_value< uint_fast32_t >("DensityFunction:buffer size",
                                              100),
            Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:anchor"),
                  params.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:sides")),
            params.get_value< CoordinateVector< uint_fast32_t > >(
                "DensityGrid:number of cells"),
            log) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~BufferedCMacIonizeSnapshotDensityFunction() {}

  /**
   * @brief Initialize the internal buffer.
   */
  virtual void initialize() {
    _buffer.resize(_buffer_size * _mapped_subgrid_size);
    _buffer_indices.resize(_number_of_subgrids.x() * _number_of_subgrids.y() *
                               _number_of_subgrids.z(),
                           0xffffffff);
  }

  /**
   * @brief Close the HDF5 file and free the buffer.
   */
  virtual void free() {
    HDF5Tools::close_group(_particle_group);
    HDF5Tools::close_file(_file);
    _buffer.clear();
    _buffer_timestamps.clear();
    _buffer_element_locks.clear();
    _buffer_subgrid_indices.clear();
    _buffer_indices.clear();
  }

  /**
   * @brief Buffer the subgrid with the given index.
   *
   * This function uses its own lock to ensure thread safe access to the HDF5
   * file.
   *
   * @param subgrid_index Index of the subgrid to buffer.
   * @return Index within the buffer of the buffered subgrid. The corresponding
   * element is locked and cannot be altered by any other thread until
   * unlock_buffer_element() is called.
   */
  inline uint_fast32_t buffer_subgrid(const uint_fast32_t subgrid_index) {

    // sort the buffers according to their last access time
    const std::vector< uint_fast32_t > buffer_timestamps_copy(
        _buffer_timestamps);
    const std::vector< uint_fast32_t > timesort =
        Utilities::argsort(buffer_timestamps_copy);
    // try to lock an old buffer
    uint_fast32_t ibuffer = 0;
    while (ibuffer < timesort.size() &&
           !_buffer_element_locks[timesort[ibuffer]].try_lock()) {
      ++ibuffer;
    }
    if (ibuffer == timesort.size()) {
      cmac_error("Unable to obtain a free subgrid buffer!");
    }

    // buffer_index is now locked and can be overwritten
    const uint_fast32_t buffer_index = timesort[ibuffer];

    // we are going to read the HDF5 file, so from this point we need to be
    // thread-safe
    _buffer_lock.lock();

    if (_log) {
      _log->write_info("Reading subgrid ", subgrid_index,
                       " into buffer element ", buffer_index);
    }

    const uint_fast32_t subgrid_offset = subgrid_index * _original_subgrid_size;

    std::vector< double > number_density;
    if (_read_number_density) {
      number_density = HDF5Tools::read_dataset_part< double >(
          _particle_group, "NumberDensity", subgrid_offset,
          _original_subgrid_size);
    } else {
      number_density = HDF5Tools::read_dataset_part< double >(
          _particle_group, "Density", subgrid_offset, _original_subgrid_size);
    }
    std::vector< double > temperature;
    if (_read_temperature) {
      temperature = HDF5Tools::read_dataset_part< double >(
          _particle_group, "Temperature", subgrid_offset,
          _original_subgrid_size);
    } else {
      temperature = HDF5Tools::read_dataset_part< double >(
          _particle_group, "Pressure", subgrid_offset, _original_subgrid_size);
    }
    std::vector< std::vector< double > > neutral_fractions(
        NUMBER_OF_IONNAMES,
        std::vector< double >(_original_subgrid_size, 1.e-6));
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      // skip ionic fractions that do not exist
      if (_read_ionic_fraction[i]) {
        neutral_fractions[i] = HDF5Tools::read_dataset_part< double >(
            _particle_group, "NeutralFraction" + get_ion_name(i),
            subgrid_offset, _original_subgrid_size);
      }
    }
    std::vector< CoordinateVector<> > velocities(_original_subgrid_size);
    if (_read_velocity) {
      cmac_warning("Not reading velocities for now!");
    }

    // we are done reading the file, unlock the file so that another thread
    // can access it
    _buffer_lock.unlock();

    const double norm = 1. / (_number_of_old_cells_per_new_cell_1D *
                              _number_of_old_cells_per_new_cell_1D *
                              _number_of_old_cells_per_new_cell_1D);
    for (uint_fast32_t ix = 0; ix < _mapped_subgrid_ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _mapped_subgrid_ncell.y(); ++iy) {
        for (uint_fast32_t iz = 0; iz < _mapped_subgrid_ncell.z(); ++iz) {
          const uint_fast32_t mapped_subgrid_index =
              ix * _mapped_subgrid_ncell.y() * _mapped_subgrid_ncell.z() +
              iy * _mapped_subgrid_ncell.z() + iz;
          double cell_number_density = 0.;
          double cell_temperature = 0.;
          double cell_ionic_fraction[NUMBER_OF_IONNAMES] = {0.};
          CoordinateVector<> cell_velocity;
          for (uint_fast32_t oix = 0;
               oix < _number_of_old_cells_per_new_cell_1D; ++oix) {
            for (uint_fast32_t oiy = 0;
                 oiy < _number_of_old_cells_per_new_cell_1D; ++oiy) {
              for (uint_fast32_t oiz = 0;
                   oiz < _number_of_old_cells_per_new_cell_1D; ++oiz) {
                const uint_fast32_t original_subgrid_index =
                    (_number_of_old_cells_per_new_cell_1D * ix + oix) *
                        _original_subgrid_ncell.y() *
                        _original_subgrid_ncell.z() +
                    (_number_of_old_cells_per_new_cell_1D * iy + oiy) *
                        _original_subgrid_ncell.z() +
                    (_number_of_old_cells_per_new_cell_1D * iz + oiz);
                if (_read_number_density) {
                  cell_number_density += number_density[original_subgrid_index];
                } else {
                  cell_number_density +=
                      number_density[original_subgrid_index] /
                      PhysicalConstants::get_physical_constant(
                          PHYSICALCONSTANT_PROTON_MASS);
                }
                if (_read_temperature) {
                  cell_temperature += temperature[original_subgrid_index];
                } else {
                  const double kB = PhysicalConstants::get_physical_constant(
                      PHYSICALCONSTANT_BOLTZMANN);
                  const double mu =
                      0.5 *
                      (1. + neutral_fractions[ION_H_n][original_subgrid_index]);
                  cell_temperature +=
                      mu / (number_density[original_subgrid_index] * kB);
                }
                for (int_fast32_t j = 0; j < NUMBER_OF_IONNAMES; ++j) {
                  cell_ionic_fraction[j] +=
                      neutral_fractions[j][original_subgrid_index];
                }
                cell_velocity += velocities[original_subgrid_index];
              }
            }
          }
          DensityValues &cell = _buffer[buffer_index * _mapped_subgrid_size +
                                        mapped_subgrid_index];
          cell.set_number_density(cell_number_density * norm);
          cell.set_temperature(cell_temperature * norm);
          for (int_fast32_t j = 0; j < NUMBER_OF_IONNAMES; ++j) {
            cell.set_ionic_fraction(j, cell_ionic_fraction[j] * norm);
          }
          cell.set_velocity(cell_velocity * norm);
        }
      }
    }
    // invalidate the old subgrid pointer, this buffer will no longer contain
    // that subgrid
    _buffer_indices[_buffer_subgrid_indices[buffer_index]] = 0xffffffff;
    // point the new subgrid to this buffer
    _buffer_indices[subgrid_index] = buffer_index;
    // make sure the cross check works
    _buffer_subgrid_indices[buffer_index] = subgrid_index;

    return buffer_index;
  }

  /**
   * @brief Obtain the index within the buffer of the subgrid with the given
   * index.
   *
   * @param subgrid_index Subgrid index.
   * @return Index within the buffer. The corresponding element is locked and
   * cannot be altered by any other thread until unlock_buffer_element() is
   * called.
   */
  inline uint_fast32_t get_buffer_index(const uint_fast32_t subgrid_index) {

    // first retrieve the index of the subgrid in the buffer
    uint_fast32_t buffer_index = _buffer_indices[subgrid_index];
    // check if the subgrid was buffered
    if (buffer_index == 0xffffffff) {
      // subgrid was not buffered, buffer it
      buffer_index = buffer_subgrid(subgrid_index);
    } else {
      // subgrid might be buffered
      // we need to lock it and check that is wasn't swapped out for another
      // subgrid before we obtained the lock
      _buffer_element_locks[buffer_index].lock();
      if (_buffer_subgrid_indices[buffer_index] != subgrid_index) {
        // too bad, another thread swapped it out! We need to release the lock
        // and buffer it again
        _buffer_element_locks[buffer_index].unlock();
        buffer_index = buffer_subgrid(subgrid_index);
      }
    }
    // update the access time for the buffer element
    cpucycle_tick(_buffer_timestamps[buffer_index]);
    return buffer_index;
  }

  /**
   * @brief Unlock the given buffer element.
   *
   * @param buffer_index Index of a buffer element.
   */
  inline void unlock_buffer_element(const uint_fast32_t buffer_index) {
    // update the access time for the buffer element
    cpucycle_tick(_buffer_timestamps[buffer_index]);
    _buffer_element_locks[buffer_index].unlock();
  }

  /**
   * @brief Get the DensityValues for the given cell.
   *
   * This function will first determine in which block the given cell resides.
   * It will then try to obtain that block and use it to satisfy the query. If
   * the requested block is not available, it will be buffered.
   * This function is thread safe. If the block is not available, it might
   * however prevent other blocks from continuing by locking the buffer.
   *
   * @param cell Cell for which the density needs to be computed.
   * @return Initial values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {

    // find the (original) subgrid that contains the cell midpoint position
    const CoordinateVector<> p = cell.get_cell_midpoint();
    const uint_fast32_t six = static_cast< uint_fast32_t >(
        (p.x() - _old_anchor.x()) / _subgrid_width.x());
    const uint_fast32_t siy = static_cast< uint_fast32_t >(
        (p.y() - _old_anchor.y()) / _subgrid_width.y());
    const uint_fast32_t siz = static_cast< uint_fast32_t >(
        (p.z() - _old_anchor.z()) / _subgrid_width.z());
    const uint_fast32_t subgrid_index =
        six * _number_of_subgrids.y() * _number_of_subgrids.z() +
        siy * _number_of_subgrids.z() + siz;

    //    if(_log){
    //      _log->write_info("p: ", p.x(), " ", p.y(), " ", p.z());
    //      _log->write_info("si: ", six, " ", siy, " ", siz);
    //      _log->write_info("subgrid_index: ", subgrid_index);
    //    }

    // obtain the buffer that contains the subgrid (and lock it)
    const uint_fast32_t buffer_index = get_buffer_index(subgrid_index);

    //    if(_log){
    //      _log->write_info("buffer index: ", buffer_index);
    //    }

    // find the index of the cell within the (mapped) buffer
    const CoordinateVector<> subgrid_anchor =
        _old_anchor + CoordinateVector<>(six * _subgrid_width.x(),
                                         siy * _subgrid_width.y(),
                                         siz * _subgrid_width.z());
    const uint_fast32_t cix = static_cast< uint_fast32_t >(
        (p.x() - subgrid_anchor.x()) / _subgrid_width.x() *
        _mapped_subgrid_ncell.x());
    const uint_fast32_t ciy = static_cast< uint_fast32_t >(
        (p.y() - subgrid_anchor.y()) / _subgrid_width.y() *
        _mapped_subgrid_ncell.y());
    const uint_fast32_t ciz = static_cast< uint_fast32_t >(
        (p.z() - subgrid_anchor.z()) / _subgrid_width.z() *
        _mapped_subgrid_ncell.z());

    //    if(_log){
    //      _log->write_info("subgrid anchor: ", subgrid_anchor.x(), " ",
    //                       subgrid_anchor.y(), " ",
    //                       subgrid_anchor.z());
    //      _log->write_info("ci: ", cix, " ", ciy, " ", ciz);
    //    }

    // obtain the cell values
    const DensityValues values =
        _buffer[buffer_index * _mapped_subgrid_size +
                cix * _mapped_subgrid_ncell.y() * _mapped_subgrid_ncell.z() +
                ciy * _mapped_subgrid_ncell.z() + ciz];

    // unlock the buffer
    unlock_buffer_element(buffer_index);

    return values;
  }
};

#endif // BUFFEREDCMACIONIZESNAPSHOTDENSITYFUNCTION_HPP
