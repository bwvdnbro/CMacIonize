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
 * @file MPIUtilities.hpp
 *
 * @brief General MPI utility methods.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MPIUTILITIES_HPP
#define MPIUTILITIES_HPP

#include "Configuration.hpp"

#ifdef HAVE_MPI
#include <mpi.h>

/**
 * @brief General MPI utility methods.
 */
namespace MPIUtilities {
/**
 * @brief Template function that returns the MPI_Datatype corresponding to the
 * given template data type.
 *
 * This function needs to be specialized for every data type.
 *
 * @return MPI_Datatype for the given template data type.
 */
template < typename _datatype_ > static inline MPI_Datatype get_datatype();

/**
 * @brief Template function that returns the MPI_Datatype corresponding to the
 * given template data type.
 *
 * Specialization for a double precision floating point value.
 *
 * @return MPI_DOUBLE.
 */
template <> inline MPI_Datatype get_datatype< double >() { return MPI_DOUBLE; }

/**
 * @brief Template function that returns the MPI_Datatype corresponding to the
 * given template data type.
 *
 * Specialization for an unsigned integer value.
 *
 * @return MPI_UNSIGNED.
 */
template <> inline MPI_Datatype get_datatype< unsigned int >() {
  return MPI_UNSIGNED;
}

/**
 * @brief Template function that returns the MPI_Datatype corresponding to the
 * given template data type.
 *
 * Specialization for an unsigned long integer value.
 *
 * @return MPI_UNSIGNED_LONG.
 */
template <> inline MPI_Datatype get_datatype< unsigned long >() {
  return MPI_UNSIGNED_LONG;
}

/**
 * @brief Template function that returns the MPI_Datatype corresponding to the
 * given template data type.
 *
 * Specialization for a signed integer value.
 *
 * @return MPI_INT.
 */
template <> inline MPI_Datatype get_datatype< int >() { return MPI_INT; }
}

#endif // HAVE_MPI

#endif // MPIUTILITIES_HPP
