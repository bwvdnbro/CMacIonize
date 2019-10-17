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
 * @file testCMICLibrary.c
 *
 * @brief Unit test for the CMI C library.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "cmi_c_library.h"
#include <stdio.h>
#include <stdlib.h>

/*! @brief Number of particles in a single coordinate direction (change as you
 *  wish). */
#define TEST_CMICLIBRARY_NPART1D 10

/*! @brief Number of particles in a plane of the particle grid. */
#define TEST_CMICLIBRARY_NPART2D                                               \
  (TEST_CMICLIBRARY_NPART1D * TEST_CMICLIBRARY_NPART1D)

/*! @brief Total number of particles. */
#define TEST_CMICLIBRARY_NPART3D                                               \
  (TEST_CMICLIBRARY_NPART1D * TEST_CMICLIBRARY_NPART2D)

/*! @brief Grid spacing. */
#define TEST_CMICLIBRARY_DX (1. / TEST_CMICLIBRARY_NPART1D)

/**
 * @brief Unit test for the CMI C library.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /* declare variables */
  double *x, *y, *z, *h, *m, *nH;
  double box_anchor[3], box_sides[3];
  unsigned int i, ix, iy, iz;
  FILE *file;
  /* pc in m */
  const double pc = 3.086e16;

  /* allocate buffers */
  x = (double *)calloc(TEST_CMICLIBRARY_NPART3D, sizeof(double));
  y = (double *)calloc(TEST_CMICLIBRARY_NPART3D, sizeof(double));
  z = (double *)calloc(TEST_CMICLIBRARY_NPART3D, sizeof(double));
  h = (double *)calloc(TEST_CMICLIBRARY_NPART3D, sizeof(double));
  m = (double *)calloc(TEST_CMICLIBRARY_NPART3D, sizeof(double));
  nH = (double *)calloc(TEST_CMICLIBRARY_NPART3D, sizeof(double));

  /* set the box dimensions to the same dimensions as the ionization simulation
   * box in test_CMI_library.param */
  box_anchor[0] = -5. * pc;
  box_anchor[1] = -5. * pc;
  box_anchor[2] = -5. * pc;
  box_sides[0] = 10. * pc;
  box_sides[1] = 10. * pc;
  box_sides[2] = 10. * pc;
  /* set up the SPH particles on a regular Cartesian 10x10x10 grid */
  for (ix = 0; ix < TEST_CMICLIBRARY_NPART1D; ++ix) {
    for (iy = 0; iy < TEST_CMICLIBRARY_NPART1D; ++iy) {
      for (iz = 0; iz < TEST_CMICLIBRARY_NPART1D; ++iz) {
        i = ix * TEST_CMICLIBRARY_NPART2D + iy * TEST_CMICLIBRARY_NPART1D + iz;
        x[i] = box_anchor[0] + TEST_CMICLIBRARY_DX * (ix + 0.5) * box_sides[0];
        y[i] = box_anchor[1] + TEST_CMICLIBRARY_DX * (iy + 0.5) * box_sides[1];
        z[i] = box_anchor[2] + TEST_CMICLIBRARY_DX * (iz + 0.5) * box_sides[2];
        h[i] = 2. * TEST_CMICLIBRARY_DX * box_sides[0];
        /* 100. cm^-3 * (10.pc)^3 * 1.67*10^{-27} kg / 1000 */
        m[i] = 4.9e33 / TEST_CMICLIBRARY_NPART3D;
      }
    }
  }

  /* initialize the library */
  cmi_init_periodic_dp("test_CMI_library.param", 1, 1., 1., box_anchor,
                       box_sides, "Petkova", 0);

  /* run the simulation */
  cmi_compute_neutral_fraction_dp(x, y, z, h, m, nH, TEST_CMICLIBRARY_NPART3D);

  /* write an output file for visual checking */
  file = fopen("test_CMI_C_library.txt", "w");
  for (i = 0; i < TEST_CMICLIBRARY_NPART3D; ++i) {
    fprintf(file, "%g %g %g %g\n", x[i], y[i], z[i], nH[i]);
  }
  fclose(file);

  /* clean up the library */
  cmi_destroy();

  free(x);
  free(y);
  free(z);
  free(h);
  free(m);
  free(nH);

  return 0;
}
