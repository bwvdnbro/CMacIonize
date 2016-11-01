################################################################################
# This file is part of CMacIonize
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CMacIonize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMacIonize is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
################################################################################

execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
                OUTPUT_VARIABLE GIT_BUILD_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# We want to get the entire date in one command, to have a single time stamp.
set(DATE_STRING
    "+Day: %-d, Month: %-m, Year: %Y, Hour: %-H, Minutes: %-M, Seconds: %-S")
execute_process(COMMAND date ${DATE_STRING}
                OUTPUT_VARIABLE FULL_DATE
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# We extract the various components and store them in separate variables
string(REGEX REPLACE ".*Day: ([0-9]*).*" "\\1" COMPILATION_TIME_DAY
       ${FULL_DATE})
string(REGEX REPLACE ".*Month: ([0-9]*).*" "\\1" COMPILATION_TIME_MONTH
       ${FULL_DATE})
string(REGEX REPLACE ".*Year: ([0-9]*).*" "\\1" COMPILATION_TIME_YEAR
       ${FULL_DATE})
string(REGEX REPLACE ".*Hour: ([0-9]*).*" "\\1" COMPILATION_TIME_HOUR
       ${FULL_DATE})
string(REGEX REPLACE ".*Minutes: ([0-9]*).*" "\\1" COMPILATION_TIME_MINUTES
       ${FULL_DATE})
string(REGEX REPLACE ".*Seconds: ([0-9]*).*" "\\1" COMPILATION_TIME_SECONDS
       ${FULL_DATE})

# Gather OS information
execute_process(COMMAND uname --operating-system
                OUTPUT_VARIABLE OS_NAME
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND uname --kernel-name
                OUTPUT_VARIABLE OS_KERNEL_NAME
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND uname --kernel-release
                OUTPUT_VARIABLE OS_KERNEL_RELEASE
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND uname --kernel-version
                OUTPUT_VARIABLE OS_KERNEL_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND uname --machine
                OUTPUT_VARIABLE OS_HARDWARE_NAME
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND uname --nodename
                OUTPUT_VARIABLE OS_HOST_NAME
                OUTPUT_STRIP_TRAILING_WHITESPACE)

configure_file(${PROJECT_SOURCE_DIR}/src/CompilerInfo.cpp.in
               ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp @ONLY)
