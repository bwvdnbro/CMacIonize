/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Signals.hpp
 *
 * @brief Operating system independent signal handlers.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SIGNALS_HPP
#define SIGNALS_HPP

#include "Error.hpp"

/**
 * @brief Operating system independent signal handlers.
 */
namespace Signals {

/**
 * @brief Signal handler for SIGINT (interrupt signal).
 */
inline void signal_interrupt_handler() {
  cmac_error("\nCTRL+C interrupt detected!");
}
} // namespace Signals

#endif // SIGNALS_HPP
