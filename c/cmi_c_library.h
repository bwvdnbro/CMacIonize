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
 * @file cmi_c_library.h
 *
 * @brief CMI library C header file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CMI_C_LIBRARY_H
#define CMI_C_LIBRARY_H

#include <stddef.h>

void cmi_init(const char *parameter_file, const int num_thread,
              const double unit_length_in_SI, const double unit_mass_in_SI,
              const char *mapping_type, const int talk);
void cmi_init_periodic_dp(const char *parameter_file, const int num_thread,
                          const double unit_length_in_SI,
                          const double unit_mass_in_SI,
                          const double *box_anchor, const double *box_sides,
                          const char *mapping_type, const int talk);
void cmi_init_periodic_sp(const char *parameter_file, const int num_thread,
                          const double unit_length_in_SI,
                          const double unit_mass_in_SI, const float *box_anchor,
                          const float *box_sides, const char *mapping_type,
                          const int talk);
void cmi_destroy();

void cmi_compute_neutral_fraction_dp(const double *x, const double *y,
                                     const double *z, const double *h,
                                     const double *m, double *nH,
                                     const size_t N);
void cmi_compute_neutral_fraction_mp(const double *x, const double *y,
                                     const double *z, const float *h,
                                     const float *m, float *nH, const size_t N);
void cmi_compute_neutral_fraction_sp(const float *x, const float *y,
                                     const float *z, const float *h,
                                     const float *m, float *nH, const size_t N);

#endif // CMI_C_LIBRARY_H
