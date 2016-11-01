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

configure_file(${PROJECT_SOURCE_DIR}/src/CompilerInfo.cpp.in
               ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp @ONLY)
