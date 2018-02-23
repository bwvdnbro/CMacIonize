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
  x = (double *)calloc(1000, sizeof(double));
  y = (double *)calloc(1000, sizeof(double));
  z = (double *)calloc(1000, sizeof(double));
  h = (double *)calloc(1000, sizeof(double));
  m = (double *)calloc(1000, sizeof(double));
  nH = (double *)calloc(1000, sizeof(double));

  /* set the box dimensions to the same dimensions as the ionization simulation
   * box in test_CMI_library.param */
  box_anchor[0] = -5. * pc;
  box_anchor[1] = -5. * pc;
  box_anchor[2] = -5. * pc;
  box_sides[0] = 10. * pc;
  box_sides[1] = 10. * pc;
  box_sides[2] = 10. * pc;
  /* set up the SPH particles on a regular Cartesian 10x10x10 grid */
  for (ix = 0; ix < 10; ++ix) {
    for (iy = 0; iy < 10; ++iy) {
      for (iz = 0; iz < 10; ++iz) {
        i = ix * 100 + iy * 10 + iz;
        x[i] = box_anchor[0] + 0.1 * (ix + 0.5) * box_sides[0];
        y[i] = box_anchor[1] + 0.1 * (iy + 0.5) * box_sides[1];
        z[i] = box_anchor[2] + 0.1 * (iz + 0.5) * box_sides[2];
        h[i] = 0.2 * box_sides[0];
        /* 100. cm^-3 * (10.pc)^3 * 1.67*10^{-27} kg / 1000 */
        m[i] = 4.9e30;
      }
    }
  }

  /* initialize the library */
  cmi_init_periodic_dp("test_CMI_library.param", 1, 1., 1., box_anchor,
                       box_sides);

  /* run the simulation */
  cmi_compute_neutral_fraction_dp(x, y, z, h, m, nH, 1000);

  /* write an output file for visual checking */
  file = fopen("test_CMI_C_library.txt", "w");
  for (i = 0; i < 1000; ++i) {
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
