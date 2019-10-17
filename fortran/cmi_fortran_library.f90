!-------------------------------------------------------------------------------
! This file is part of CMacIonize
! Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
!
! CMacIonize is free software: you can redistribute it and/or modify
! it under the terms of the GNU Affero General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CMacIonize is distributed in the hope that it will be useful,
! but WITOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Affero General Public License for more details.
!
! You should have received a copy of the GNU Affero General Public License
! along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------

!-
!> @file cmi_fortran_library.f90
!>
!> @brief Fortran interface to the CMI library.
!>
!> @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
!-

!-
!> @brief Fortran exposure of the CMI library.
!-
module cmi_fortran_library

  use iso_c_binding
  implicit none

  !-
  !> @brief Subroutines declared in the CMI library that need to be exposed to
  !> Fortran.
  !-
  interface c_subroutines

    !-
    !> @brief Fortran interface for CMILibrary::cmi_init().
    !>
    !> @param parameter_file Name of the parameter file to use.
    !> @param num_thread Number of shared memory parallel threads to use.
    !> @param unit_length_in_SI Internal length unit (in m).
    !> @param unit_mass_in_SI Internal mass unit (in kg).
    !> @param mapping_type Type of density mapping to use.
    !> @param talk Print output to the terminal window?
    !-
    subroutine cmi_init_c(parameter_file, num_thread, unit_length_in_SI, &
                          unit_mass_in_SI, mapping_type, talk) &
      bind(C, name = "cmi_init")

      use iso_c_binding
      implicit none

      character (kind = c_char), intent(in) :: parameter_file
      integer (kind = c_int), intent(in), value :: num_thread
      real (kind = c_double), intent(in), value :: unit_length_in_SI
      real (kind = c_double), intent(in), value :: unit_mass_in_SI
      character (kind = c_char), intent(in) :: mapping_type
      integer (kind = c_int), intent(in), value :: talk

    end subroutine cmi_init_c

    !-
    !> @brief Fortran interface for CMILibrary::cmi_init_periodic_dp().
    !>
    !> Double precision version.
    !>
    !> @param parameter_file Name of the parameter file to use.
    !> @param num_thread Number of shared memory parallel threads to use.
    !> @param unit_length_in_SI Internal length unit (in m).
    !> @param unit_mass_in_SI Internal mass unit (in kg).
    !> @param box_anchor Coordinates of the left front bottom corner of the
    !> simulation box (in internal length units).
    !> @param box_sides Side lengths of the simulation box (in internal length
    !> units).
    !> @param mapping_type Type of density mapping to use.
    !> @param talk Print output to the terminal window?
    !-
    subroutine cmi_init_periodic_dp_c(parameter_file, num_thread, &
                                      unit_length_in_SI, unit_mass_in_SI, &
                                      box_anchor, box_sides, mapping_type, &
                                      talk) &
      bind(C, name = "cmi_init_periodic_dp")

      use iso_c_binding
      implicit none

      character (kind = c_char), intent(in) :: parameter_file
      integer (kind = c_int), intent(in), value :: num_thread
      real (kind = c_double), intent(in), value :: unit_length_in_SI
      real (kind = c_double), intent(in), value :: unit_mass_in_SI
      real (kind = c_double), intent(in) :: box_anchor(3)
      real (kind = c_double), intent(in) :: box_sides(3)
      character (kind = c_char), intent(in) :: mapping_type
      integer (kind = c_int), intent(in), value :: talk

    end subroutine cmi_init_periodic_dp_c

    !-
    !> @brief Fortran interface for CMILibrary::cmi_init_periodic_dp().
    !>
    !> Single precision version.
    !>
    !> @param parameter_file Name of the parameter file to use.
    !> @param num_thread Number of shared memory parallel threads to use.
    !> @param unit_length_in_SI Internal length unit (in m).
    !> @param unit_mass_in_SI Internal mass unit (in kg).
    !> @param box_anchor Coordinates of the left front bottom corner of the
    !> simulation box (in internal length units).
    !> @param box_sides Side lengths of the simulation box (in internal length
    !> units).
    !> @param mapping_type Type of density mapping to use.
    !> @param talk Print output to the terminal window?
    !-
    subroutine cmi_init_periodic_sp_c(parameter_file, num_thread, &
                                      unit_length_in_SI, unit_mass_in_SI, &
                                      box_anchor, box_sides, mapping_type, &
                                      talk) &
      bind(C, name = "cmi_init_periodic_sp")

      use iso_c_binding
      implicit none

      character (kind = c_char), intent(in) :: parameter_file
      integer (kind = c_int), intent(in), value :: num_thread
      real (kind = c_double), intent(in), value :: unit_length_in_SI
      real (kind = c_double), intent(in), value :: unit_mass_in_SI
      real (kind = c_float), intent(in) :: box_anchor(3)
      real (kind = c_float), intent(in) :: box_sides(3)
      character (kind = c_char), intent(in) :: mapping_type
      integer (kind = c_int), intent(in), value :: talk

    end subroutine cmi_init_periodic_sp_c

    !-
    !> @brief Fortran interface for CMILibrary::cmi_destroy().
    !-
    subroutine cmi_destroy() bind(C, name = "cmi_destroy")
    end subroutine cmi_destroy

    !-
    !> @brief Fortran interface for
    !> CMILibrary::cmi_compute_neutral_fraction_dp().
    !>
    !> @param x X coordinates (in internal length units).
    !> @param y Y coordinates (in internal length units).
    !> @param z Z coordinates (in internal length units).
    !> @param h Smoothing lengths (in internal length units).
    !> @param m Masses (in internal mass units).
    !> @param nH Neutral fraction array to compute.
    !> @param N Size of all arrays.
    !-
    subroutine cmi_compute_neutral_fraction_dp(x, y, z, h, m, nH, N) &
      bind(C, name = "cmi_compute_neutral_fraction_dp")

      use iso_c_binding
      implicit none

      real (kind = c_double), intent(in) :: x(N)
      real (kind = c_double), intent(in) :: y(N)
      real (kind = c_double), intent(in) :: z(N)
      real (kind = c_double), intent(in) :: h(N)
      real (kind = c_double), intent(in) :: m(N)
      real (kind = c_double), intent(inout) :: nH(N)
      integer (kind = c_size_t), intent(in), value :: N

    end subroutine cmi_compute_neutral_fraction_dp

    !-
    !> @brief Fortran interface for
    !> CMILibrary::cmi_compute_neutral_fraction_mp().
    !>
    !> @param x X coordinates (in internal length units).
    !> @param y Y coordinates (in internal length units).
    !> @param z Z coordinates (in internal length units).
    !> @param h Smoothing lengths (in internal length units).
    !> @param m Masses (in internal mass units).
    !> @param nH Neutral fraction array to compute.
    !> @param N Size of all arrays.
    !-
    subroutine cmi_compute_neutral_fraction_mp(x, y, z, h, m, nH, N) &
      bind(C, name = "cmi_compute_neutral_fraction_mp")

      use iso_c_binding
      implicit none

      real (kind = c_double), intent(in) :: x(N)
      real (kind = c_double), intent(in) :: y(N)
      real (kind = c_double), intent(in) :: z(N)
      real (kind = c_float), intent(in) :: h(N)
      real (kind = c_float), intent(in) :: m(N)
      real (kind = c_float), intent(inout) :: nH(N)
      integer (kind = c_size_t), intent(in), value :: N

    end subroutine cmi_compute_neutral_fraction_mp

    !-
    !> @brief Fortran interface for
    !> CMILibrary::cmi_compute_neutral_fraction_sp().
    !>
    !> @param x X coordinates (in internal length units).
    !> @param y Y coordinates (in internal length units).
    !> @param z Z coordinates (in internal length units).
    !> @param h Smoothing lengths (in internal length units).
    !> @param m Masses (in internal mass units).
    !> @param nH Neutral fraction array to compute.
    !> @param N Size of all arrays.
    !-
    subroutine cmi_compute_neutral_fraction_sp(x, y, z, h, m, nH, N) &
      bind(C, name = "cmi_compute_neutral_fraction_sp")

      use iso_c_binding
      implicit none

      real (kind = c_float), intent(in) :: x(N)
      real (kind = c_float), intent(in) :: y(N)
      real (kind = c_float), intent(in) :: z(N)
      real (kind = c_float), intent(in) :: h(N)
      real (kind = c_float), intent(in) :: m(N)
      real (kind = c_float), intent(inout) :: nH(N)
      integer (kind = c_size_t), intent(in), value :: N

    end subroutine cmi_compute_neutral_fraction_sp

  end interface c_subroutines

  contains

    !-
    !> @brief Wrapper for cmi_init_c().
    !>
    !> We need to convert the Fortran string into a C string before we can pass
    !> it on to the C function.
    !>
    !> @param parameter_file Name of the parameter file to use.
    !> @param num_thread Number of shared memory parallel threads to use.
    !> @param unit_length_in_SI Internal length unit (in m).
    !> @param unit_mass_in_SI Internal mass unit (in kg).
    !> @param mapping_type Type of density mapping to use.
    !> @param talk Print output to the terminal window?
    !-
    subroutine cmi_init(parameter_file, num_thread, unit_length_in_SI, &
                        unit_mass_in_SI, mapping_type, talk)

      use iso_c_binding
      implicit none

      character (len = *), intent(in) :: parameter_file
      integer, intent(in) :: num_thread
      real*8, intent(in) :: unit_length_in_SI
      real*8, intent(in) :: unit_mass_in_SI
      character (len = *), intent(in) :: mapping_type
      integer, intent(in) :: talk

      call cmi_init_c(trim(parameter_file)//C_NULL_CHAR, num_thread, &
                      unit_length_in_SI, unit_mass_in_SI, &
                      trim(mapping_type)//C_NULL_CHAR, talk)

    end subroutine cmi_init

    !-
    !> @brief Wrapper for cmi_init_periodic_dp_c().
    !>
    !> We need to convert the Fortran string into a C string before we can pass
    !> it on to the C function.
    !>
    !> @param parameter_file Name of the parameter file to use.
    !> @param num_thread Number of shared memory parallel threads to use.
    !> @param unit_length_in_SI Internal length unit (in m).
    !> @param unit_mass_in_SI Internal mass unit (in kg).
    !> @param box_anchor Coordinates of the left front bottom corner of the
    !> simulation box (in internal length units).
    !> @param box_sides Side lengths of the simulation box (in internal length
    !> units).
    !> @param mapping_type Type of density mapping to use.
    !> @param talk Print output to the terminal window?
    !-
    subroutine cmi_init_periodic_dp(parameter_file, num_thread, &
                                    unit_length_in_SI, unit_mass_in_SI, &
                                    box_anchor, box_sides, mapping_type, &
                                    talk)

      use iso_c_binding
      implicit none

      character (len = *), intent(in) :: parameter_file
      integer, intent(in) :: num_thread
      real*8, intent(in) :: unit_length_in_SI
      real*8, intent(in) :: unit_mass_in_SI
      real*8, intent(in) :: box_anchor(3)
      real*8, intent(in) :: box_sides(3)
      character (len = *), intent(in) :: mapping_type
      integer, intent(in) :: talk

      call cmi_init_periodic_dp_c(trim(parameter_file)//C_NULL_CHAR, &
                                  num_thread, unit_length_in_SI, &
                                  unit_mass_in_SI, box_anchor, box_sides, &
                                  trim(mapping_type)//C_NULL_CHAR, talk)

    end subroutine cmi_init_periodic_dp

    !-
    !> @brief Wrapper for cmi_init_periodic_sp_c().
    !>
    !> We need to convert the Fortran string into a C string before we can pass
    !> it on to the C function.
    !>
    !> @param parameter_file Name of the parameter file to use.
    !> @param num_thread Number of shared memory parallel threads to use.
    !> @param unit_length_in_SI Internal length unit (in m).
    !> @param unit_mass_in_SI Internal mass unit (in kg).
    !> @param box_anchor Coordinates of the left front bottom corner of the
    !> simulation box (in internal length units).
    !> @param box_sides Side lengths of the simulation box (in internal length
    !> units).
    !> @param mapping_type Type of density mapping to use.
    !> @param talk Print output to the terminal window?
    !-
    subroutine cmi_init_periodic_sp(parameter_file, num_thread, &
                                    unit_length_in_SI, unit_mass_in_SI, &
                                    box_anchor, box_sides, mapping_type, &
                                    talk)

      use iso_c_binding
      implicit none

      character (len = *), intent(in) :: parameter_file
      integer, intent(in) :: num_thread
      real*8, intent(in) :: unit_length_in_SI
      real*8, intent(in) :: unit_mass_in_SI
      real*4, intent(in) :: box_anchor(3)
      real*4, intent(in) :: box_sides(3)
      character (len = *), intent(in) :: mapping_type
      integer, intent(in) :: talk

      call cmi_init_periodic_sp_c(trim(parameter_file)//C_NULL_CHAR, &
                                  num_thread, unit_length_in_SI, &
                                  unit_mass_in_SI, box_anchor, box_sides, &
                                  trim(mapping_type)//C_NULL_CHAR, talk)

    end subroutine cmi_init_periodic_sp

end module cmi_fortran_library
