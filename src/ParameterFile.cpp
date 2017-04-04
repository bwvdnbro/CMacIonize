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
 * @file ParameterFile.cpp
 *
 * @brief Parameter file: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "ParameterFile.hpp"
#include "Error.hpp"
#include <fstream>

/**
 * @brief Constructor.
 *
 * Reads in the contents of the file and stores them in the internal dictionary.
 *
 * @param filename Name of the parameter file.
 */
ParameterFile::ParameterFile(std::string filename) {
  std::ifstream file(filename);

  if (!file) {
    cmac_error("Failed to open parameter file \"%s\"", filename.c_str());
  }

  _yaml_dictionary = YAMLDictionary(file);
}

/**
 * @brief Print the contents of the internal dictionary to the given stream.
 *
 * This routine is meant to reproduce the parameter file that was actually used,
 * containing all parameters, and not only the ones that were present in the
 * original parameter file.
 *
 * This function calls the print_contents() routine of the underlying YAML
 * dictionary, but adds a time stamp to the first line of the file.
 *
 * @param stream std::ostream to write to.
 */
void ParameterFile::print_contents(std::ostream &stream) const {

  stream << "# file written on " << Utilities::get_timestamp() << ".\n";

  _yaml_dictionary.print_contents(stream, true);
}
