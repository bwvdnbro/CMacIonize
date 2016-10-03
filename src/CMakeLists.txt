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

# Rules for files that need to be configured.
message(STATUS "Configuring external data files...")
# LineCoolingDataLocation.hpp
configure_file(${PROJECT_SOURCE_DIR}/data/atom4.dat
               ${PROJECT_BINARY_DIR}/data/atom4.dat COPYONLY)
set(LINECOOLINGDATALOCATION ${PROJECT_BINARY_DIR}/data/atom4.dat)
message(STATUS "LineCoolingDataLocation.hpp: ${LINECOOLINGDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/LineCoolingDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/LineCoolingDataLocation.hpp @ONLY)
# CrossSectionsDataLocation.hpp
configure_file(${PROJECT_SOURCE_DIR}/data/verner_data.txt
               ${PROJECT_BINARY_DIR}/data/verner_data.txt COPYONLY)
set(CROSSSECTIONSDATALOCATION ${PROJECT_BINARY_DIR}/data/verner_data.txt)
message(STATUS "CrossSectionsDataLocation.hpp: ${CROSSSECTIONSDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/CrossSectionsDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/CrossSectionsDataLocation.hpp @ONLY)

set(CMACIONIZE_SOURCES
    CMacIonize.cpp
    LineCoolingData.cpp

    LineCoolingData.hpp
    LineCoolingDataLocation.hpp
)

add_executable(CMacIonize ${CMACIONIZE_SOURCES})
