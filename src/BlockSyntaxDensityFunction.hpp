/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file BlockSyntaxDensityFunction.hpp
 *
 * @brief DensityFunction implementation that constructs a density field based
 * on geometrical building blocks specified in a YAML file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BLOCKSYNTAXDENSITYFUNCTION_HPP
#define BLOCKSYNTAXDENSITYFUNCTION_HPP

#include "BlockSyntaxBlock.hpp"
#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "YAMLDictionary.hpp"

#include <fstream>
#include <sstream>

/**
 * @brief DensityFunction implementation that constructs a density field based
 * on geometrical building blocks specified in a YAML file.
 */
class BlockSyntaxDensityFunction : public DensityFunction {
private:
  /*! @brief Geometrical building blocks. */
  std::vector< BlockSyntaxBlock > _blocks;

  /**
   * @brief Get the exponent corresponding to a given block type.
   *
   * @param type Type of block.
   * @return Exponent corresponding to that type.
   */
  inline static double get_exponent(std::string type) {
    if (type == "rhombus") {
      return 1.;
    } else if (type == "sphere") {
      return 2.;
    } else if (type == "cube") {
      return 10.;
    } else {
      cmac_error("Unknown block type: \"%s\"!", type.c_str());
      return 0.;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the YAML file containing the block information.
   * @param log Log to write logging info to.
   */
  BlockSyntaxDensityFunction(std::string filename, Log *log = nullptr) {
    std::ifstream file(filename);

    if (!file) {
      cmac_error("Error while opening file \"%s\"!", filename.c_str());
    }

    YAMLDictionary blockfile(file);

    const int numblock = blockfile.get_value< int >("number of blocks");
    for (int i = 0; i < numblock; ++i) {
      std::stringstream blockname;
      blockname << "block[" << i << "]:";
      CoordinateVector<> origin =
          blockfile.get_physical_vector< QUANTITY_LENGTH >(blockname.str() +
                                                           "origin");
      CoordinateVector<> sides =
          blockfile.get_physical_vector< QUANTITY_LENGTH >(blockname.str() +
                                                           "sides");
      std::string type =
          blockfile.get_value< std::string >(blockname.str() + "type");
      double exponent = get_exponent(type);
      double density = blockfile.get_physical_value< QUANTITY_NUMBER_DENSITY >(
          blockname.str() + "number density");
      double temperature = blockfile.get_physical_value< QUANTITY_TEMPERATURE >(
          blockname.str() + "initial temperature");
      CoordinateVector<> velocity =
          blockfile.get_physical_vector< QUANTITY_VELOCITY >(
              blockname.str() + "initial velocity",
              "[0. m s^-1, 0. m s^-1, 0. m s^-1]");
      if (density < 0.) {
        cmac_error("Negative density (%g) given for block %i!", density, i);
      }
      if (temperature < 0.) {
        cmac_error("Negative temperature (%g) given for block %i!", temperature,
                   i);
      }
      _blocks.push_back(BlockSyntaxBlock(origin, sides, exponent, density,
                                         temperature, velocity));
    }

    if (log) {
      log->write_status("Created BlockSyntaxDensityFunction with ", numblock,
                        " blocks.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read.
   * @param log Log to write logging info to.
   */
  BlockSyntaxDensityFunction(ParameterFile &params, Log *log = nullptr)
      : BlockSyntaxDensityFunction(
            params.get_value< std::string >("densityfunction:filename"), log) {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * Due to the way this function is written, the values for the last block
   * containing the given position are used. This means the order in which
   * nested blocks are given is important!
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) const {
    DensityValues values;

    const CoordinateVector<> position = cell.get_cell_midpoint();

    double density = -1.;
    double temperature = -1.;
    CoordinateVector<> velocity;
    for (unsigned int i = 0; i < _blocks.size(); ++i) {
      if (_blocks[i].is_inside(position)) {
        density = _blocks[i].get_number_density();
        temperature = _blocks[i].get_temperature();
        velocity = _blocks[i].get_velocity();
      }
    }
    if (density < 0.) {
      cmac_error("No block found containing position [%g m, %g m, %g m]!",
                 position.x(), position.y(), position.z());
    }
    if (temperature < 0.) {
      cmac_error("No block found containing position [%g m, %g m, %g m]!",
                 position.x(), position.y(), position.z());
    }

    values.set_number_density(density);
    values.set_temperature(temperature);
    values.set_ionic_fraction(ION_H_n, 1.e-6);
    values.set_ionic_fraction(ION_He_n, 1.e-6);
    values.set_velocity(velocity);
    return values;
  }
};

#endif // BLOCKSYNTAXDENSITYFUNCTION_HPP
