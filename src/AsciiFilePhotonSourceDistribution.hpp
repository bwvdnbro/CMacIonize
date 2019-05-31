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
 * @file AsciiFilePhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution read from an ASCII text file in YAML format.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ASCIIFILEPHOTONSOURCEDISTRIBUTION_HPP
#define ASCIIFILEPHOTONSOURCEDISTRIBUTION_HPP

#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "YAMLDictionary.hpp"

class Log;

/**
 * @brief PhotonSourceDistribution read from an ASCII text file in YAML format.
 */
class AsciiFilePhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Positions of the sources (in m). */
  std::vector< CoordinateVector<> > _source_positions;

  /*! @brief Ionising luminosities of the sources (in s^-1). */
  std::vector< double > _source_luminosities;

  /*! @brief Total ionising luminosity of all sources (in s^-1). */
  double _total_luminosity;

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the file that contains the source information.
   */
  inline AsciiFilePhotonSourceDistribution(const std::string filename) {

    std::ifstream file(filename);

    if (!file) {
      cmac_error("Error while opening file \"%s\"!", filename.c_str());
    }

    YAMLDictionary blocks(file);
    const uint_fast32_t number_of_sources =
        blocks.get_value< uint_fast32_t >("number of sources");
    _source_positions.resize(number_of_sources);
    _source_luminosities.resize(number_of_sources);
    _total_luminosity = 0.;
    for (uint_fast32_t i = 0; i < number_of_sources; ++i) {
      std::stringstream blockname;
      blockname << "source[" << i << "]:";

      _source_positions[i] = blocks.get_physical_vector< QUANTITY_LENGTH >(
          blockname.str() + "position");
      _source_luminosities[i] = blocks.get_physical_value< QUANTITY_FREQUENCY >(
          blockname.str() + "luminosity");
      _total_luminosity += _source_luminosities[i];
    }

    std::ofstream ofile(filename + ".used-values");
    blocks.print_contents(ofile, true);
    ofile.close();
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline AsciiFilePhotonSourceDistribution(ParameterFile &params, Log *log)
      : AsciiFilePhotonSourceDistribution(params.get_filename(
            "PhotonSourceDistribution:filename", "sources.yml")) {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_sources() const {
    return _source_positions.size();
  }

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {
    return _source_positions[index];
  }

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the photon source, used to determine how many photons are
   * emitted from this particular source.
   */
  virtual double get_weight(photonsourcenumber_t index) const {
    return _source_luminosities[index] / _total_luminosity;
  }

  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */
  virtual double get_total_luminosity() const { return _total_luminosity; }

  /**
   * @brief Write the distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    const size_t number_of_sources = _source_positions.size();
    restart_writer.write(number_of_sources);
    for (size_t i = 0; i < number_of_sources; ++i) {
      _source_positions[i].write_restart_file(restart_writer);
      restart_writer.write(_source_luminosities[i]);
    }
    restart_writer.write(_total_luminosity);
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline AsciiFilePhotonSourceDistribution(RestartReader &restart_reader) {

    const size_t number_of_sources = restart_reader.read< size_t >();
    _source_positions.resize(number_of_sources);
    _source_luminosities.resize(number_of_sources);
    for (size_t i = 0; i < number_of_sources; ++i) {
      _source_positions[i] = CoordinateVector<>(restart_reader);
      _source_luminosities[i] = restart_reader.read< double >();
    }
    _total_luminosity = restart_reader.read< double >();
  }
};

#endif // ASCIIFILEPHOTONSOURCEDISTRIBUTION_HPP
