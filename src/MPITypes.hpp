/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file MPITypes.hpp
 *
 * @brief Macros to deal with C++ integer types.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MPITYPES_HPP
#define MPITYPES_HPP

#include "Configuration.hpp"

#ifdef HAVE_MPI

#include <cinttypes>
#include <mpi.h>

// char types

/*! @brief Macro wrapper providing the right MPI_Datatype for an int_fast8_t. */
#define MPI_INT_FAST8_T (MPITypes::get_datatype< int_fast8_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for an
 *  int_least8_t. */
#define MPI_INT_LEAST8_T (MPITypes::get_datatype< int_least8_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for a uint_fast8_t. */
#define MPI_UINT_FAST8_T (MPITypes::get_datatype< uint_fast8_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for a
 *  uint_least8_t. */
#define MPI_UINT_LEAST8_T (MPITypes::get_datatype< uint_least8_t >())

// integer types

/*! @brief Macro wrapper providing the right MPI_Datatype for an
 *  int_fast32_t. */
#define MPI_INT_FAST32_T (MPITypes::get_datatype< int_fast32_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for an
 *  int_least32_t. */
#define MPI_INT_LEAST32_T (MPITypes::get_datatype< int_least32_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for a
 *  uint_fast32_t. */
#define MPI_UINT_FAST32_T (MPITypes::get_datatype< uint_fast32_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for a
 *  uint_least32_t. */
#define MPI_UINT_LEAST32_T (MPITypes::get_datatype< uint_least32_t >())

// long integer types

/*! @brief Macro wrapper providing the right MPI_Datatype for an
 *  int_fast64_t. */
#define MPI_INT_FAST64_T (MPITypes::get_datatype< int_fast64_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for an
 *  int_least64_t. */
#define MPI_INT_LEAST64_T (MPITypes::get_datatype< int_least64_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for a
 *  uint_fast64_t. */
#define MPI_UINT_FAST64_T (MPITypes::get_datatype< uint_fast64_t >())

/*! @brief Macro wrapper providing the right MPI_Datatype for a
 *  uint_least64_t. */
#define MPI_UINT_LEAST64_T (MPITypes::get_datatype< uint_least64_t >())

/**
 * @brief Auxiliary functions to obtain the right MPI_Datatype for special
 * C++ integer types.
 */
namespace MPITypes {

/**
 * @brief Auxiliary template function prototype.
 *
 * @return The MPI_Datatype for the given template _integer_type_.
 */
template < typename _integer_type_ > inline MPI_Datatype get_datatype();

/**
 * @brief MPITypes::get_datatype() specialization for an unsigned char.
 *
 * @return MPI_UNSIGNED_CHAR.
 */
template <> inline MPI_Datatype get_datatype< unsigned char >() {
  return MPI_UNSIGNED_CHAR;
}

/**
 * @brief MPITypes::get_datatype() specialization for a signed char.
 *
 * @return MPI_SIGNED_CHAR.
 */
template <> inline MPI_Datatype get_datatype< char >() {
  return MPI_SIGNED_CHAR;
}

/**
 * @brief MPITypes::get_datatype() specialization for a short.
 *
 * @return MPI_SHORT.
 */
template <> inline MPI_Datatype get_datatype< short >() { return MPI_SHORT; }

/**
 * @brief MPITypes::get_datatype() specialization for an unsigned short.
 *
 * @return MPI_UNSIGNED_SHORT.
 */
template <> inline MPI_Datatype get_datatype< unsigned short >() {
  return MPI_UNSIGNED_SHORT;
}

/**
 * @brief MPITypes::get_datatype() specialization for an int.
 *
 * @return MPI_INT.
 */
template <> inline MPI_Datatype get_datatype< int >() { return MPI_INT; }

/**
 * @brief MPITypes::get_datatype() specialization for an unsigned int.
 *
 * @return MPI_UNSIGNED.
 */
template <> inline MPI_Datatype get_datatype< unsigned int >() {
  return MPI_UNSIGNED;
}

/**
 * @brief MPITypes::get_datatype() specialization for a long.
 *
 * @return MPI_LONG.
 */
template <> inline MPI_Datatype get_datatype< long >() { return MPI_LONG; }

/**
 * @brief MPITypes::get_datatype() specialization for an unsigned long.
 *
 * @return MPI_UNSIGNED_LONG.
 */
template <> inline MPI_Datatype get_datatype< unsigned long >() {
  return MPI_UNSIGNED_LONG;
}

/**
 * @brief MPITypes::get_datatype() specialization for a long long.
 *
 * @return MPI_LONG_LONG.
 */
template <> inline MPI_Datatype get_datatype< long long >() {
  return MPI_LONG_LONG;
}

/**
 * @brief MPITypes::get_datatype() specialization for an unsigned long long.
 *
 * @return MPI_UNSIGNED_LONG_LONG.
 */
template <> inline MPI_Datatype get_datatype< unsigned long long >() {
  return MPI_UNSIGNED_LONG_LONG;
}

} // namespace MPITypes

#endif // HAVE_MPI

#endif // MPITYPES_HPP
