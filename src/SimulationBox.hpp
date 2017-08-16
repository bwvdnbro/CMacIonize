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
 * @file SimulationBox.hpp
 *
 * @brief Simulation box.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SIMULATIONBOX_HPP
#define SIMULATIONBOX_HPP

#include "Box.hpp"
#include "ParameterFile.hpp"

/**
 * @brief Simulation box.
 */
class SimulationBox {
private:
  /*! @brief Simulation box (in m). */
  const Box<> _box;

  /*! @brief Periodicity flags. */
  const CoordinateVector< bool > _periodicity;

public:
  /**
   * @brief Constructor.
   *
   * @param box Simulation box (in m).
   * @param periodicity Periodicity flags.
   */
  SimulationBox(const Box<> box, const CoordinateVector< bool > periodicity)
      : _box(box), _periodicity(periodicity) {}

  /**
   * @brief ParameterFile constructor.
   *
   * @param parameter_file ParameterFile to read from.
   */
  SimulationBox(ParameterFile &parameter_file)
      : SimulationBox(
            Box<>(parameter_file.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:anchor", "[-5. pc, -5. pc, -5. pc]"),
                  parameter_file.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:sides", "[10. pc, 10. pc, 10. pc]")),
            parameter_file.get_value< CoordinateVector< bool > >(
                "SimulationBox:periodicity", CoordinateVector< bool >(false))) {
  }

  /**
   * @brief Get the simulation box.
   *
   * @return Simulation box (in m).
   */
  const Box<> &get_box() const { return _box; }

  /**
   * @brief Get the periodicity flags.
   *
   * @return Periodicity flags.
   */
  const CoordinateVector< bool > &get_periodicity() const {
    return _periodicity;
  }
};

#endif // SIMULATIONBOX_HPP
