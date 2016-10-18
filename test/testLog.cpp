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
 * @file testLog.cpp
 *
 * @brief Unit test for the Log class and its implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "FileLog.hpp"
#include "TerminalLog.hpp"

/**
 * @brief Unit test for the Log class and its implementations.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  TerminalLog tlog(LOGLEVEL_INFO);

  tlog.info("Info message.");
  tlog.status("Status message.");
  tlog.warning("Warning message.");
  tlog.error("Error message.");

  FileLog flog("test.log", LOGLEVEL_INFO);

  flog.info("Info message.");
  flog.status("Status message.");
  flog.warning("Warning message.");
  flog.error("Error message.");

  return 0;
}
