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
!> @file testCMIFortranLibrary.f90
!>
!> @brief Unit test for the CMI Fortran library.
!>
!> @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
!-

!-
!> @brief Unit test for the CMI Fortran library.
!>
!> @return Exit code: 0 on success.
!-
program testCMIFortranLibrary

  use cmi_fortran_library
  implicit none

  real*8 x(1000), y(1000), z(1000), h(1000), m(1000), nH(1000)
  real*8 box_anchor(3), box_sides(3)
  real*8 pc
  integer i, ix, iy, iz

  ! pc in m
  pc = 3.086d16
  ! simulation box dimensions (these are the same as in test_CMI_library.param)
  box_anchor(1) = -5. * pc
  box_anchor(2) = -5. * pc
  box_anchor(3) = -5. * pc
  box_sides(1) = 10. * pc
  box_sides(2) = 10. * pc
  box_sides(3) = 10. * pc
  ! initialize a 10x10x10 grid of particles
  do ix = 1, 10
    do iy = 1, 10
      do iz = 1, 10
        i = (ix-1) * 100 + (iy-1) * 10 + iz
        x(i) = box_anchor(1) + 0.1 * (ix - 0.5) * box_sides(1)
        y(i) = box_anchor(2) + 0.1 * (iy - 0.5) * box_sides(2)
        z(i) = box_anchor(3) + 0.1 * (iz - 0.5) * box_sides(3)
        h(i) = 0.2 * box_sides(1)
        ! 100. cm^-3 * (10.pc)^3 * 1.67*10^{-27} kg / 1000
        m(i) = 4.9d30
      end do
    end do
  end do

  ! initialize the library
  call cmi_init_periodic_dp("test_CMI_library.param", 1, 1.d0, 1.d0, box_anchor, &
                            box_sides, "Petkova", 0)

  ! run the simulation
  call cmi_compute_neutral_fraction_dp(x, y, z, h, m, nH, int8(1000))

  ! write an output file for visual checking
  open(unit = 1, file = "test_CMI_fortran_library.txt")
  do i = 1, 1000
    write(1, *) x(i), y(i), z(i), nH(i)
  end do
  close(1)

  ! clean up the library
  call cmi_destroy()

end program testCMIFortranLibrary
