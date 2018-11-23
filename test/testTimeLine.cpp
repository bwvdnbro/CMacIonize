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
 * @file testTimeLine.cpp
 *
 * @brief Unit test for the TimeLine class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "TimeLine.hpp"

#include <cinttypes>

/**
 * @brief Unit test for the TimeLine class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// basic test
  {
    cmac_status("Basic time line test...");
    TimeLine timeline(0., 1., 0.01, 0.1);
    double actual_timestep, current_time;
    uint_fast32_t numstep = 1;
    while (timeline.advance(0.02, actual_timestep, current_time)) {
      cmac_status("Step: %g %g.", actual_timestep, current_time);
      assert_condition(actual_timestep == 0.015625);
      assert_condition(current_time == numstep * 0.015625);
      ++numstep;
    }
    cmac_status("Step: %g %g.", actual_timestep, current_time);
    assert_condition(actual_timestep == 0.015625);
    assert_condition(current_time == 1.);
    cmac_status("Took %" PRIuFAST32 " steps.", numstep);
    assert_condition(numstep == 64);
    cmac_status("Done.");
  }

  /// more advanced test:
  /// same as basic test, but with offset
  {
    cmac_status("Advanced time line test with offset...");
    TimeLine timeline(0.77, 1.77, 0.01, 0.1);
    double actual_timestep, current_time;
    uint_fast32_t numstep = 1;
    while (timeline.advance(0.02, actual_timestep, current_time)) {
      cmac_status("Step: %g %g.", actual_timestep, current_time);
      assert_condition(actual_timestep == 0.015625);
      assert_condition(current_time == 0.77 + numstep * 0.015625);
      ++numstep;
    }
    cmac_status("Step: %g %g.", actual_timestep, current_time);
    assert_condition(actual_timestep == 0.015625);
    assert_condition(current_time == 1.77);
    cmac_status("Took %" PRIuFAST32 " steps.", numstep);
    assert_condition(numstep == 64);
    cmac_status("Done.");
  }

  /// more advanced test:
  /// same as basic test, but with variable requested time step
  {
    cmac_status("Advanced time line test with variable time step...");
    TimeLine timeline(0., 1., 0.005, 0.1);

    /// restart test
    {
      RestartWriter restart_writer("timeline.dump");
      timeline.write_restart_file(restart_writer);
    }

    double actual_timestep;
    double current_time = 0.;
    uint_fast32_t numstep = 1;
    while (timeline.advance(0.02 + 0.015 * std::sin(2. * M_PI * current_time),
                            actual_timestep, current_time)) {
      cmac_status("Step: %g %g.", actual_timestep, current_time);
      ++numstep;
    }
    cmac_status("Step: %g %g.", actual_timestep, current_time);
    assert_condition(current_time == 1.);
    cmac_status("Took %" PRIuFAST32 " steps.", numstep);
    assert_condition(numstep == 107);
    cmac_status("Done.");

    /// restart test
    {
      RestartReader restart_reader("timeline.dump");
      TimeLine restart_timeline(restart_reader);

      current_time = 0.;
      uint_fast32_t numstep = 1;
      while (restart_timeline.advance(
          0.02 + 0.015 * std::sin(2. * M_PI * current_time), actual_timestep,
          current_time)) {
        ++numstep;
      }
      assert_condition(numstep == 107);
    }
  }

  return 0;
}
