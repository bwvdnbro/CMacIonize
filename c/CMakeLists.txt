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

# Copy the library header file
configure_file(${PROJECT_SOURCE_DIR}/c/cmi_c_library.h
               ${PROJECT_BINARY_DIR}/include/cmi_c_library.h COPYONLY)

# Create files that contain the necessary includes and linker dependencies to
# compile a program that uses the CMI C library
# An example compilation could look like this:
#  gcc -o test `cat ${PROJECT_BINARY_DIR}/compilation/cmi_c_libs.txt` test.c \
#    `cat ${PROJECT_BINARY_DIR}/compilation/cmi_c_includes.txt`
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(CLIBS hdf5 mpi-cxx)
  string(REPLACE ";" " " CLIBS_LDFLAGS_STRING "${CLIBS_STATIC_LDFLAGS}")
  configure_file(${PROJECT_SOURCE_DIR}/c/cmi_c_libs.txt
                 ${PROJECT_BINARY_DIR}/compilation/cmi_c_libs.txt @ONLY)
  configure_file(${PROJECT_SOURCE_DIR}/c/cmi_c_includes.txt
                 ${PROJECT_BINARY_DIR}/compilation/cmi_c_includes.txt @ONLY)
endif(PKG_CONFIG_FOUND)
