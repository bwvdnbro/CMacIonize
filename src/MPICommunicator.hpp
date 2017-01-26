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
 * @file MPICommunicator.hpp
 *
 * @brief C++ wrapper around basic MPI functions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MPICOMMUNICATOR_HPP
#define MPICOMMUNICATOR_HPP

#include "Error.hpp"

#include <mpi.h>
#include <vector>

#if defined(__clang__)
// header containing functions to temporarily disable the address sanitizer
#include <sanitizer/lsan_interface.h>
#endif

/*! @brief Temporarily disable instrumentation for memory allocations done
 *  between the call to this macro, and the call to NO_LEAK_CHECK_END. */
#if defined(__clang__)
#define NO_LEAK_CHECK_BEGIN __lsan_disable();
#else
#define NO_LEAK_CHECK_BEGIN
#endif

/*! @brief Re-enable instrumentation after it was disabled by
 *  NO_LEAK_CHECK_BEGIN. */
#if defined(__clang__)
#define NO_LEAK_CHECK_END __lsan_enable();
#else
#define NO_LEAK_CHECK_END
#endif

/*! @brief Default size of the reduction buffer, if no size is provided. */
#define MPICOMMUNICATOR_DEFAULT_BUFFERSIZE 1000000

/**
 * @brief Aliases for MPI_Op types.
 */
enum MPIOperatorType {
  /*! @brief Take the sum of a variable across all processes. */
  MPI_SUM_OF_ALL_PROCESSES = 0
};

/**
 * @brief C++ wrapper around basic MPI functions.
 */
class MPICommunicator {
private:
  /*! @brief Rank of the local MPI process. */
  int _rank;

  /*! @brief Total number of MPI processes. */
  int _size;

public:
  /**
   * @brief Constructor.
   *
   * Calls MPI_Init(), sets up custom error handling, and initializes rank and
   * size variables.
   *
   * @param argc Number of command line arguments passed on to the main program.
   * @param argv Command line arguments passed on to the main program.
   */
  inline MPICommunicator(int &argc, char **argv) {
    // MPI_Init is known to cause memory leak detections, so we disable the
    // address sanitizer for all allocations made by it
    NO_LEAK_CHECK_BEGIN
    int status = MPI_Init(&argc, &argv);
    NO_LEAK_CHECK_END
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to initialize MPI!");
    }

    // make sure errors are handled by us, not by the MPI library
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    // get the size and rank
    status = MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to obtain rank of MPI process!");
    }
    status = MPI_Comm_size(MPI_COMM_WORLD, &_size);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to obtain number of MPI processes!");
    }
  }

  /**
   * @brief Destructor.
   *
   * Calls MPI_Finalize().
   */
  inline ~MPICommunicator() {
    int status = MPI_Finalize();
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to clean up MPI!");
    }
  }

  /**
   * @brief Get the rank of the local MPI process.
   *
   * @return Rank of the local MPI process.
   */
  inline int get_rank() const { return _rank; }

  /**
   * @brief Get the total number of MPI processes.
   *
   * @return Total number of MPI processes.
   */
  inline int get_size() const { return _size; }

  /**
   * @brief Distribute the given number across all processes, so that the sum of
   * the returned values across all processes is the given number.
   *
   * @param number Number to distribute.
   * @return Part of the number that is assigned to this process.
   */
  inline unsigned int distribute(const unsigned int number) const {
    unsigned int quotient = number / _size;
    int remainder = number % _size;
    // all processes with a rank smaller than the remainder get one element
    // extra (since _rank < remainder evaluates to either 0 or 1)
    return quotient + (_rank < remainder);
  }

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
   * @brief Function that returns the MPI_Op corresponding to the given
   * MPIOperatorType.
   *
   * @param type MPIOperatorType.
   * @return MPI_Op corresponding to the given MPIOperatorType.
   */
  inline static MPI_Op get_operator(MPIOperatorType type) {
    switch (type) {
    case MPI_SUM_OF_ALL_PROCESSES:
      return MPI_SUM;
    default:
      cmac_error("Unknown MPIOperatorType: %i!", type);
      return 0;
    }
  }

  /**
   * @brief Reduce the elements of the given std::vector of objects, using the
   * given getter member function to obtain an object data member to reduce, and
   * the given setter member function to set the result of the reduction.
   *
   * @param v std::vector containing objects of a given template class type.
   * @param getter Member function of the given template class type that returns
   * a value of the given template data type that will be reduced across all
   * processes.
   * @param setter Member function of the given template class type that takes
   * a value of the given template data type and stores it in the class type
   * object after the reduction.
   */
  template < MPIOperatorType _operatortype_, typename _datatype_,
             typename _classtype_ >
  void reduce(std::vector< _classtype_ > &v,
              _datatype_ (_classtype_::*getter)() const,
              void (_classtype_::*setter)(_datatype_)) const {
    // we only communicate if there are multiple processes
    if (_size > 1) {
      std::vector< _datatype_ > sendbuffer(v.size());
      for (unsigned int i = 0; i < v.size(); ++i) {
        sendbuffer[i] = (v[i].*(getter))();
      }
      MPI_Datatype dtype = get_datatype< _datatype_ >();
      MPI_Op otype = get_operator(_operatortype_);
      int status =
          MPI_Allreduce(MPI_IN_PLACE, &sendbuffer[0], sendbuffer.size(), dtype,
                        otype, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
      for (unsigned int i = 0; i < v.size(); ++i) {
        (v[i].*(setter))(sendbuffer[i]);
      }
    }
  }

  /**
   * @brief Reduce the elements pointed to by the given begin and end iterator,
   * using the given getter member function to obtain an object data member to
   * reduce, and the given setter member function to set the result of the
   * reduction.
   *
   * This routine also accepts extra arguments that will be passed on to the
   * getter and setter. The setter function call should accept these arguments
   * after the value of the reduced variable, like this:
   * \code{.cpp}
   *   setter(reduced_value, args)
   * \endcode
   *
   * @param begin Iterator to the first element that should be reduced.
   * @param end Iterator to the first element that should not be reduced, or the
   * end of the list.
   * @param getter Member function of the given template class type that returns
   * a value of the given template data type that will be reduced across all
   * processes.
   * @param setter Member function of the given template class type that takes
   * a value of the given template data type and stores it in the class type
   * object after the reduction.
   * @param size Number of elements to reduce in a single MPI communication.
   * This value sets the memory size of the buffer that is used internally.
   * @param args Extra arguments passed on to the getter and setter.
   */
  template < MPIOperatorType _operatortype_, typename _datatype_,
             typename _classtype_, typename _iteratortype_,
             typename... _arguments_ >
  void reduce(_iteratortype_ begin, _iteratortype_ end,
              _datatype_ (_classtype_::*getter)(_arguments_...) const,
              void (_classtype_::*setter)(_datatype_, _arguments_...),
              unsigned int size, _arguments_... args) const {
    if (_size > 1) {
      if (size == 0) {
        size = MPICOMMUNICATOR_DEFAULT_BUFFERSIZE;
      }
      std::vector< _datatype_ > sendbuffer(size);
      MPI_Datatype dtype = get_datatype< _datatype_ >();
      MPI_Op otype = get_operator(_operatortype_);
      // we loop over all elements of the iterator
      _iteratortype_ it = begin;
      while (it != end) {
        // depending on the given buffer size, the reduction might be split in a
        // number of blocks
        _iteratortype_ blockit = it;
        unsigned int i = 0;
        // this loop ends if there are no more elements to reduce, or if the
        // communication buffer size is reached
        while (blockit != end && i < size) {
          sendbuffer[i] = ((*blockit).*(getter))(args...);
          ++i;
          ++blockit;
        }
        // we reduce i elements, since that is the number of elements we added
        // to the send buffer above
        // by providing MPI_IN_PLACE as sendbuffer argument, we tell MPI to the
        // an in place reduction, reading and storing from the recvbuffer (which
        // is our sendbuffer).
        int status = MPI_Allreduce(MPI_IN_PLACE, &sendbuffer[0], i, dtype,
                                   otype, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS) {
          cmac_error("Error in MPI_Allreduce!");
        }
        // we reset the counters to the same values used in the first loop
        i = 0;
        blockit = it;
        // the condition here is exactly the same as the first loop, so the same
        // elements will be traversed in the same order
        while (blockit != end && i < size) {
          ((*blockit).*(setter))(sendbuffer[i], args...);
          ++i;
          ++blockit;
        }
        // make sure the next block starts where the current one ended
        it = blockit;
      }
    }
  }

  /**
   * @brief Reduce the given variable across all processes.
   *
   * @param value Variable to reduce.
   * @return Reduced variable.
   */
  template < MPIOperatorType _operatortype_, typename _datatype_ >
  _datatype_ reduce(_datatype_ value) const {
    if (_size > 1) {
      MPI_Datatype dtype = get_datatype< _datatype_ >();
      MPI_Op otype = get_operator(_operatortype_);
      int status =
          MPI_Allreduce(MPI_IN_PLACE, &value, 1, dtype, otype, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
    }
    return value;
  }

  /**
   * @brief Reduce the given array of the given template size across all
   * processes.
   *
   * The given array is replaced with the resulting reduced array (the reduction
   * is actually done in place).
   *
   * @param array Array to reduce.
   */
  template < MPIOperatorType _operatortype_, unsigned int _size_,
             typename _datatype_ >
  void reduce(_datatype_ *array) const {
    if (_size > 1) {
      MPI_Datatype dtype = get_datatype< _datatype_ >();
      MPI_Op otype = get_operator(_operatortype_);
      int status = MPI_Allreduce(MPI_IN_PLACE, array, _size_, dtype, otype,
                                 MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
    }
  }

  /**
   * @brief Reduce the given std::vector across all processes.
   *
   * @param vector std::vector to reduce.
   */
  template < MPIOperatorType _operatortype_, typename _datatype_ >
  void reduce(std::vector< _datatype_ > &vector) {
    if (_size > 1) {
      MPI_Datatype dtype = get_datatype< _datatype_ >();
      MPI_Op otype = get_operator(_operatortype_);
      int status = MPI_Allreduce(MPI_IN_PLACE, vector.data(), vector.size(),
                                 dtype, otype, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
    }
  }
};

/**
 * @brief Template function that returns the MPI_Datatype corresponding to the
 * given template data type.
 *
 * Specialization for a double precision floating point value.
 *
 * @return MPI_DOUBLE.
 */
template <> inline MPI_Datatype MPICommunicator::get_datatype< double >() {
  return MPI_DOUBLE;
}

/**
 * @brief Template function that returns the MPI_Datatype corresponding to the
 * given template data type.
 *
 * Specialization for an unsigned integer value.
 *
 * @return MPI_UNSIGNED.
 */
template <>
inline MPI_Datatype MPICommunicator::get_datatype< unsigned int >() {
  return MPI_UNSIGNED;
}

#endif // MPICOMMUNICATOR_HPP
