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
 * This is the only file that is aware of the existence of MPI.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MPICOMMUNICATOR_HPP
#define MPICOMMUNICATOR_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "MPIMessage.hpp"
#include "MPIMessageBox.hpp"
#include "MPIUtilities.hpp"

#include <vector>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_POSIX
#include <fstream>
#include <sstream>
#include <unistd.h>
#endif

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
  int_fast32_t _rank;

  /*! @brief Total number of MPI processes. */
  int_fast32_t _size;

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

#ifdef HAVE_MPI
    // MPI_Init is known to cause memory leak detections, so we disable the
    // address sanitizer for all allocations made by it
    NO_LEAK_CHECK_BEGIN
    int_fast32_t status = MPI_Init(&argc, &argv);
    NO_LEAK_CHECK_END
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to initialize MPI!");
    }

    // make sure errors are handled by us, not by the MPI library
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    // get the size and rank (these are integers)
    int rank, size;
    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to obtain rank of MPI process!");
    }
    status = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to obtain number of MPI processes!");
    }
    _rank = rank;
    _size = size;
#else

    // check that we do not accidentally run the program in parallel
    // we first try an OpenMP specific check
    char *env = getenv("OMPI_COMM_WORLD_SIZE");
    if (env != nullptr) {
      const int_fast32_t size = atoi(env);
      if (size > 1) {
        cmac_error("Program is being run in parallel using MPI, but MPI was "
                   "disabled at compile time! Please compile again with MPI "
                   "support!");
      }
    } else {
// the OMPI environment variable was not found
// this means that either we are not running in MPI mode (which is fine)
// OR we are not using OpenMPI, or an old version that does not set the
// environment variable
// we try to see if the program has a parent process that might indicate
// we are running in MPI mode
// the method below will only work for POSIX systems
#ifdef HAVE_POSIX
      const pid_t ppid = getppid();
      std::stringstream procfilename;
      procfilename << "/proc/" << ppid << "/cmdline";
      std::ifstream procfile(procfilename.str());
      std::string procname;
      procfile >> procname;
      if (procname.find("mpirun") != std::string::npos ||
          procname.find("mpiexec") != std::string::npos) {
        cmac_error("Program is being run in parallel using MPI, but MPI was "
                   "disabled at compile time! Please compile again with MPI "
                   "support!");
      }
#else
      cmac_warning("Unable to check if non MPI build is being run in parallel. "
                   "Are you sure you want to run a serial version of the "
                   "code?");
#endif
    }
    // no MPI: we assume a single process with rank 0
    _rank = 0;
    _size = 1;
#endif
  }

  /**
   * @brief Destructor.
   *
   * Calls MPI_Finalize().
   */
  inline ~MPICommunicator() {
#ifdef HAVE_MPI
    const int_fast32_t status = MPI_Finalize();
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to clean up MPI!");
    }
#endif
  }

  /**
   * @brief Get the rank of the local MPI process.
   *
   * @return Rank of the local MPI process.
   */
  inline int_fast32_t get_rank() const { return _rank; }

  /**
   * @brief Get the total number of MPI processes.
   *
   * @return Total number of MPI processes.
   */
  inline int_fast32_t get_size() const { return _size; }

  /**
   * @brief Distribute the given number across all processes, so that the sum of
   * the returned values across all processes is the given number.
   *
   * @param number Number to distribute.
   * @return Part of the number that is assigned to this process.
   */
  inline uint_fast32_t distribute(uint_fast32_t number) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      const uint_fast32_t quotient = number / _size;
      const int_fast32_t remainder = number % _size;
      // all processes with a rank smaller than the remainder get one element
      // extra (since _rank < remainder evaluates to either 0 or 1)
      return quotient + (_rank < remainder);
    } else {
      return number;
    }
#else
    // why bother doing a calculation if you already know the result?
    return number;
#endif
  }

  /**
   * @brief Distribute the continuous block of indices with given begin and end
   * index across a given number of processes, and get the part for the process
   * with the given rank.
   *
   * @param rank Rank for which we want the part.
   * @param size Total number of processes.
   * @param begin Start of the block.
   * @param end End of the block.
   * @return std::pair containing the begin and end index for the part of the
   * block that is assigned to process rank.
   */
  static inline std::pair< size_t, size_t > distribute_block(int_fast32_t rank,
                                                             int_fast32_t size,
                                                             size_t begin,
                                                             size_t end) {

    const size_t block_size = end - begin;
    size_t block_begin;
    size_t block_end;
    const size_t quotient = block_size / size;
    const int_fast32_t remainder = block_size % size;
    // here is the logic: the start of the local block should be the sum of
    // all blocks on processes with ranks lower than this block
    // these have size quotient + (_rank < remainder), which means they have
    // total size _rank*quotient + the number of blocks with 1 element more
    // the end of the block is the beginning of the next block, which means we
    // have the same logic for _rank+1
    block_begin = rank * quotient + std::min(rank, remainder);
    block_end = (rank + 1) * quotient + std::min(rank + 1, remainder);
    return std::make_pair(block_begin, block_end);
  }

  /**
   * @brief Distribute the continuous block of indices with given begin and end
   * index across all processes.
   *
   * @param begin First index of the block.
   * @param end First index not part of the block.
   * @return std::pair containing the begin and end index for this MPI process.
   */
  inline std::pair< size_t, size_t > distribute_block(size_t begin,
                                                      size_t end) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      return distribute_block(_rank, _size, begin, end);
    } else {
      return std::make_pair(begin, end);
    }
#else
    return std::make_pair(begin, end);
#endif
  }

#ifdef HAVE_MPI
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
#endif

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

#ifdef HAVE_MPI
    // we only communicate if there are multiple processes
    if (_size > 1) {
      std::vector< _datatype_ > sendbuffer(v.size());
      for (size_t i = 0; i < v.size(); ++i) {
        sendbuffer[i] = (v[i].*(getter))();
      }
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();
      const MPI_Op otype = get_operator(_operatortype_);
      const int_fast32_t status =
          MPI_Allreduce(MPI_IN_PLACE, &sendbuffer[0], sendbuffer.size(), dtype,
                        otype, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
      for (size_t i = 0; i < v.size(); ++i) {
        (v[i].*(setter))(sendbuffer[i]);
      }
    }
#endif
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
              size_t size, _arguments_... args) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      if (size == 0) {
        size = MPICOMMUNICATOR_DEFAULT_BUFFERSIZE;
      }
      std::vector< _datatype_ > sendbuffer(size);
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();
      const MPI_Op otype = get_operator(_operatortype_);
      // we loop over all elements of the iterator
      _iteratortype_ it = begin;
      while (it != end) {
        // depending on the given buffer size, the reduction might be split in a
        // number of blocks
        _iteratortype_ blockit = it;
        size_t i = 0;
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
        const int_fast32_t status = MPI_Allreduce(
            MPI_IN_PLACE, &sendbuffer[0], i, dtype, otype, MPI_COMM_WORLD);
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
#endif
  }

  /**
   * @brief Reduce the given variable across all processes.
   *
   * @param value Variable to reduce.
   * @return Reduced variable.
   */
  template < MPIOperatorType _operatortype_, typename _datatype_ >
  void reduce(_datatype_ &value) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();
      const MPI_Op otype = get_operator(_operatortype_);
      const int_fast32_t status =
          MPI_Allreduce(MPI_IN_PLACE, &value, 1, dtype, otype, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
    }
#endif
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

#ifdef HAVE_MPI
    if (_size > 1) {
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();
      const MPI_Op otype = get_operator(_operatortype_);
      const int_fast32_t status = MPI_Allreduce(MPI_IN_PLACE, array, _size_,
                                                dtype, otype, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
    }
#endif
  }

  /**
   * @brief Reduce the given std::vector across all processes.
   *
   * @param vector std::vector to reduce.
   */
  template < MPIOperatorType _operatortype_, typename _datatype_ >
  void reduce(std::vector< _datatype_ > &vector) {

#ifdef HAVE_MPI
    if (_size > 1) {
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();
      const MPI_Op otype = get_operator(_operatortype_);
      const int_fast32_t status =
          MPI_Allreduce(MPI_IN_PLACE, vector.data(), vector.size(), dtype,
                        otype, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        cmac_error("Error in MPI_Allreduce!");
      }
    }
#endif
  }

  /**
   * @brief Reduce the elements pointed to by the given begin and end iterator,
   * using the given template property accessor to get and set the relevant
   * values for each element.
   *
   * @param begin Iterator to the first element that should be reduced.
   * @param end Iterator to the first element that should not be reduced, or the
   * end of the list.
   * @param size Number of elements to reduce in a single MPI communication.
   * This value sets the memory size of the buffer that is used internally.
   */
  template < MPIOperatorType _operatortype_, typename _datatype_,
             typename _PropertyAccessor_, typename _iteratortype_ >
  void reduce(_iteratortype_ begin, _iteratortype_ end, size_t size) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      if (size == 0) {
        size = MPICOMMUNICATOR_DEFAULT_BUFFERSIZE;
      }
      std::vector< _datatype_ > sendbuffer(size);
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();
      const MPI_Op otype = get_operator(_operatortype_);
      // we loop over all elements of the iterator
      _iteratortype_ it = begin;
      while (it != end) {
        // depending on the given buffer size, the reduction might be split in a
        // number of blocks
        _iteratortype_ blockit = it;
        size_t i = 0;
        // this loop ends if there are no more elements to reduce, or if the
        // communication buffer size is reached
        while (blockit != end && i < size) {
          sendbuffer[i] = _PropertyAccessor_::get_value(blockit);
          ++i;
          ++blockit;
        }
        // we reduce i elements, since that is the number of elements we added
        // to the send buffer above
        // by providing MPI_IN_PLACE as sendbuffer argument, we tell MPI to the
        // an in place reduction, reading and storing from the recvbuffer (which
        // is our sendbuffer).
        const int_fast32_t status = MPI_Allreduce(
            MPI_IN_PLACE, &sendbuffer[0], i, dtype, otype, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS) {
          cmac_error("Error in MPI_Allreduce!");
        }
        // we reset the counters to the same values used in the first loop
        i = 0;
        blockit = it;
        // the condition here is exactly the same as the first loop, so the same
        // elements will be traversed in the same order
        while (blockit != end && i < size) {
          _PropertyAccessor_::set_value(blockit, sendbuffer[i]);
          ++i;
          ++blockit;
        }
        // make sure the next block starts where the current one ended
        it = blockit;
      }
    }
#endif
  }

  /**
   * @brief Ensure the given std::vector is up to date on all processes,
   * assuming that MPI process i holds the block returned by
   * distribute_block(i, 0, vector.size()).
   *
   * @param vector std::vector to gather.
   */
  template < typename _datatype_ >
  void gather(std::vector< _datatype_ > &vector) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();
      const std::pair< size_t, size_t > local_block =
          distribute_block(0, vector.size());
      // do a complicated communication ring:
      // we do a loop with _size steps; each process sends to process
      // _rank+step, and receives from process _rank-step
      // since we do not require _size to be even, there is no way of knowing
      // up front in which order to send and receive, which means we have to do
      // at least a non-blocking send
      // we cannot do a sendrecv, since the number of elements sent and received
      // can differ by one
      // for the same reason, we cannot use allgather, which essentially does
      // what we do here, but for equal sizes on each process
      for (int_fast32_t step = 0; step < _size; ++step) {
        // the %_size is necessary to wrap our ranks in the range [0, _size[
        // the +_size in recv_block is not really necessary, but it makes sure
        // the rank is always positive (and would be necessary if _rank was an
        // unsigned integer)
        const int_fast32_t sendrank = (_rank + step) % _size;
        const int_fast32_t recvrank = (_rank + _size - step) % _size;
        const std::pair< size_t, size_t > recv_block =
            distribute_block(recvrank, _size, 0, vector.size());
        MPI_Request request;
        int_fast32_t status = MPI_Isend(
            &vector[local_block.first], local_block.second - local_block.first,
            dtype, sendrank, 0, MPI_COMM_WORLD, &request);
        if (status != MPI_SUCCESS) {
          cmac_error("Failed to issue a non-blocking send!");
        }
        MPI_Status recvstatus;
        status = MPI_Recv(&vector[recv_block.first],
                          recv_block.second - recv_block.first, dtype, recvrank,
                          0, MPI_COMM_WORLD, &recvstatus);
        if (status != MPI_SUCCESS) {
          cmac_error("Failed to receive vector block!");
        }
        status = MPI_Wait(&request, &recvstatus);
        if (status != MPI_SUCCESS) {
          cmac_error("Failed to send non-blocking vector block!");
        }
      }
    }
#endif
  }

  /**
   * @brief Gather the elements pointed to by the given begin and end iterator,
   * using the given getter member function to obtain an object data member to
   * gather, and the given setter member function to set the result of the
   * gathering.
   *
   * This routine also accepts extra arguments that will be passed on to the
   * getter and setter. The setter function call should accept these arguments
   * after the value of the reduced variable, like this:
   * \code{.cpp}
   *   setter(reduced_value, args)
   * \endcode
   *
   * This routine assumes that the process with the lowest rank holds the lowest
   * iterator values, and so on. This will be the case if the iterator values
   * where distributed using distribute_block.
   *
   * @param global_begin Iterator to the first element that should be gathered.
   * @param global_end Iterator to the first element that should not be
   * gathered, or the end of the list.
   * @param local_begin Iterator to the first local element that should be send.
   * @param local_end Iterator to the first local element that should not be
   * send, or the end of the list.
   * @param getter Member function of the given template class type that returns
   * a value of the given template data type that will be gathered across all
   * processes.
   * @param setter Member function of the given template class type that takes
   * a value of the given template data type and stores it in the class type
   * object after the gathering.
   * @param size Number of elements to gather in a single MPI communication.
   * This value sets the memory size of the buffer that is used internally.
   * @param args Extra arguments passed on to the getter and setter.
   */
  template < typename _datatype_, typename _classtype_, typename _iteratortype_,
             typename... _arguments_ >
  void gather(_iteratortype_ global_begin, _iteratortype_ global_end,
              _iteratortype_ local_begin, _iteratortype_ local_end,
              _datatype_ (_classtype_::*getter)(_arguments_...) const,
              void (_classtype_::*setter)(_datatype_, _arguments_...),
              size_t size, _arguments_... args) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      if (size == 0) {
        size = MPICOMMUNICATOR_DEFAULT_BUFFERSIZE;
      }
      std::vector< _datatype_ > sendbuffer(size);
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();

      _iteratortype_ it = global_begin;
      // we loop over all processes; if the active processor rank irank matches
      // the local rank, we send. If it does not match, we receive.
      for (int_fast32_t irank = 0; irank < _size; ++irank) {

        const bool is_local = (irank == _rank);

        // note that we use a single loop for both the local and the other
        // processes
        // this way, we can use MPI_Bcast instead of an ordinary send/recv,
        // which is more efficient for large processor numbers
        // we cannot use fast integer types, as we need to send 'done' and 'i'
        int done = false;
        if (is_local) {
          // at this point, we should have reached the local part of the
          // iterator, check this to make sure
          cmac_assert(it == local_begin);

          done = (it == local_end);
        }

        while (!done) {
          unsigned int i = 0;
          if (is_local) {
            // first fill the buffer until there are no more local elements, or
            // until the buffer is full
            while (it != local_end && i < size) {
              sendbuffer[i] = ((*it).*(getter))(args...);
              ++i;
              ++it;
            }
            done = (it == local_end);
          }

          // first communicate the number of elemens that we send
          MPI_Bcast(&i, 1, MPI_UNSIGNED, irank, MPI_COMM_WORLD);
          // now communicate the elements themselves
          MPI_Bcast(sendbuffer.data(), i, dtype, irank, MPI_COMM_WORLD);
          // finally communicate the done flag
          MPI_Bcast(&done, 1, MPI_INT, irank, MPI_COMM_WORLD);

          if (!is_local) {
            // now fill the non local buffers with the result of the
            // communication
            for (uint_fast32_t j = 0; j < i; ++j) {
              ((*it).*(setter))(sendbuffer[j], args...);
              ++it;
            }
          }
        }
      }

      // check that we did all elements
      cmac_assert(it == global_end);
    }
#endif
  }

  /**
   * @brief Gather the elements pointed to by the given begin and end iterator,
   * using the given PropertyAccessor to access iterator values.
   *
   * This routine assumes that the process with the lowest rank holds the lowest
   * iterator values, and so on. This will be the case if the iterator values
   * where distributed using distribute_block.
   *
   * @param global_begin Iterator to the first element that should be gathered.
   * @param global_end Iterator to the first element that should not be
   * gathered, or the end of the list.
   * @param local_begin Iterator to the first local element that should be send.
   * @param local_end Iterator to the first local element that should not be
   * send, or the end of the list.
   * @param size Number of elements to gather in a single MPI communication.
   * This value sets the memory size of the buffer that is used internally.
   */
  template < typename _datatype_, typename _PropertyAccessor_,
             typename _iteratortype_ >
  void gather(_iteratortype_ global_begin, _iteratortype_ global_end,
              _iteratortype_ local_begin, _iteratortype_ local_end,
              size_t size) const {

#ifdef HAVE_MPI
    if (_size > 1) {
      if (size == 0) {
        size = MPICOMMUNICATOR_DEFAULT_BUFFERSIZE;
      }
      std::vector< _datatype_ > sendbuffer(size);
      const MPI_Datatype dtype = MPIUtilities::get_datatype< _datatype_ >();

      _iteratortype_ it = global_begin;
      // we loop over all processes; if the active processor rank irank matches
      // the local rank, we send. If it does not match, we receive.
      for (int_fast32_t irank = 0; irank < _size; ++irank) {

        const bool is_local = (irank == _rank);

        // note that we use a single loop for both the local and the other
        // processes
        // this way, we can use MPI_Bcast instead of an ordinary send/recv,
        // which is more efficient for large processor numbers
        int done = false;
        if (is_local) {
          // at this point, we should have reached the local part of the
          // iterator, check this to make sure
          cmac_assert(it == local_begin);

          done = (it == local_end);
        }

        while (!done) {
          unsigned int i = 0;
          if (is_local) {
            // first fill the buffer until there are no more local elements, or
            // until the buffer is full
            while (it != local_end && i < size) {
              sendbuffer[i] = _PropertyAccessor_::get_value(it);
              ++i;
              ++it;
            }
            done = (it == local_end);
          }

          // first communicate the number of elemens that we send
          MPI_Bcast(&i, 1, MPI_UNSIGNED, irank, MPI_COMM_WORLD);
          // now communicate the elements themselves
          MPI_Bcast(sendbuffer.data(), i, dtype, irank, MPI_COMM_WORLD);
          // finally communicate the done flag
          MPI_Bcast(&done, 1, MPI_INT, irank, MPI_COMM_WORLD);

          if (!is_local) {
            // now fill the non local buffers with the result of the
            // communication
            for (uint_fast32_t j = 0; j < i; ++j) {
              _PropertyAccessor_::set_value(it, sendbuffer[j]);
              ++it;
            }
          }
        }
      }

      // check that we did all elements
      cmac_assert(it == global_end);
    }
#endif
  }

  /**
   * @brief Send the given message to the given process.
   *
   * @param destination MPI process that should receive the message.
   * @param message Pointer to the message that should be sent. Memory
   * management for this pointer will be taken over by the MPIMessageBox.
   * @param message_box MessageBox to use.
   */
  inline void send_message(int_fast32_t destination, MPIMessage *message,
                           MPIMessageBox &message_box) const {

#ifdef HAVE_MPI
    MPITagSizeMessage *announcement = message_box.get_announcement(message);
    int_fast32_t status =
        MPI_Isend(announcement->get_array_handle(), 3, MPI_INT, destination, 0,
                  MPI_COMM_WORLD, announcement->get_request_handle());
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to send MPITagSizeMessage!");
    }
    status =
        MPI_Isend(message->get_buffer_handle(), message->get_buffer_size(),
                  message->get_datatype(), destination, announcement->get_tag(),
                  MPI_COMM_WORLD, message->get_request_handle());
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to send MPIMessage!");
    }
#endif
  }

  /**
   * @brief Send a semi ready MPITagSizeMessage to the given process.
   *
   * @param destination MPI process that should receive the announcement.
   * @param message_box MessageBox to use.
   */
  inline void send_semi_ready(int_fast32_t destination,
                              MPIMessageBox &message_box) const {

#ifdef HAVE_MPI
    MPITagSizeMessage *announcement = message_box.get_announcement_semi_ready();
    const int_fast32_t status =
        MPI_Isend(announcement->get_array_handle(), 3, MPI_INT, destination, 0,
                  MPI_COMM_WORLD, announcement->get_request_handle());
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to send semi ready MPITagSizeMessage!");
    }
#endif
  }

  /**
   * @brief Send a ready MPITagSizeMessage to the given process.
   *
   * @param destination MPI process that should receive the announcement.
   * @param message_box MessageBox to use.
   */
  inline void send_ready(int_fast32_t destination,
                         MPIMessageBox &message_box) const {

#ifdef HAVE_MPI
    MPITagSizeMessage *announcement = message_box.get_announcement_ready();
    const int_fast32_t status =
        MPI_Isend(announcement->get_array_handle(), 3, MPI_INT, destination, 0,
                  MPI_COMM_WORLD, announcement->get_request_handle());
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to send ready MPITagSizeMessage!");
    }
#endif
  }

  /**
   * @brief Check if there is a pending message that can be received.
   *
   * @param message_box MPIMessageBox that generates a draft message for a given
   * message tag.
   */
  inline void check_for_message(MPIMessageBox &message_box) const {

#ifdef HAVE_MPI
    MPI_Status probestatus;
    int flag;
    int_fast32_t status =
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &probestatus);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to probe for incoming message!");
    }
    if (flag > 0) {
      const int_fast32_t source = probestatus.MPI_SOURCE;
      MPITagSizeMessage announcement;
      status = MPI_Recv(announcement.get_array_handle(), 3, MPI_INT, source, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (status != MPI_SUCCESS) {
        cmac_error("Failed to receive incoming MPITagSizeMessage!");
      }
      MPIMessage *message = message_box.get_message(source, announcement);
      if (message != nullptr) {
        status =
            MPI_Irecv(message->get_buffer_handle(), message->get_buffer_size(),
                      message->get_datatype(), source, announcement.get_tag(),
                      MPI_COMM_WORLD, message->get_request_handle());
        if (status != MPI_SUCCESS) {
          cmac_error("Failed to set up non-blocking receive for MPIMessage!");
        }
      }
    }
#endif
  }
};

#endif // MPICOMMUNICATOR_HPP
