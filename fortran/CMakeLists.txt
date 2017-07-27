################################################################################
# This file is part of CMacIonize
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

set(LIBCMIFORTRANLIBRARY_SOURCES
    cmi_fortran_library.f90
)

add_library(CMIFortranLibrary ${LIBCMIFORTRANLIBRARY_SOURCES})
target_link_libraries(CMIFortranLibrary CMILibrary)

# Create files that contain the necessary includes and linker dependencies to
# compile a program that uses the CMI Fortran library
# An example compilation could look like this:
# gcc -o test `cat ${PROJECT_BINARY_DIR}/compilation/cmi_fortran_includes.txt` \
#   test.c `cat ${PROJECT_BINARY_DIR}/compilation/cmi_fortran_libs.txt`
string(REPLACE ";" " " HDF5LIBS_STRING "${HDF5_LIBRARIES}")
string(REPLACE ";" " " MPILIBS_STRING  "${MPI_CXX_LIBRARIES}")
set(FortranLIBS_LDFLAGS_STRING "${HDF5LIBS_STRING} ${MPILIBS_STRING}")
configure_file(${PROJECT_SOURCE_DIR}/fortran/cmi_fortran_libs.txt
               ${PROJECT_BINARY_DIR}/compilation/cmi_fortran_libs.txt @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/fortran/cmi_fortran_includes.txt
               ${PROJECT_BINARY_DIR}/compilation/cmi_fortran_includes.txt
               @ONLY)
