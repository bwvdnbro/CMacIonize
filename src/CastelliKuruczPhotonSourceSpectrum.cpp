/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file CastelliKuruczPhotonSourceSpectrum.cpp
 *
 * @brief CastelliKuruczPhotonSourceSpectrum implementation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "CastelliKuruczPhotonSourceSpectrum.hpp"
#include "CastelliKuruczDataLocation.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

/**
 * @brief Interpolate the data block based on the given fractions and bounding
 * indices.
 *
 * @param fZ Metallicity axis fraction.
 * @param fT Temperature axis fraction.
 * @param fg Surface gravity axis fraction.
 * @param fl Wavelength axis fraction.
 * @param iZ Metallicity bounding index.
 * @param iT Temperature bounding index.
 * @param ig Surface gravity bounding index.
 * @param il Wavelength bounding index.
 * @param F Data block.
 * @return Interpolated value.
 */
double CastelliKuruczPhotonSourceSpectrum::interpolate(
    const double fZ, const double fT, const double fg, const double fl,
    const uint_fast32_t iZ, const uint_fast32_t iT, const uint_fast32_t ig,
    const uint_fast32_t il, const HDF5Tools::HDF5DataBlock< double, 4 > &F) {

  const double factor[16] = {
      (1. - fZ) * (1. - fT) * (1. - fg) * (1. - fl),
      (1. - fZ) * (1. - fT) * (1. - fg) * fl,
      (1. - fZ) * (1. - fT) * fg * (1. - fl),
      (1. - fZ) * (1. - fT) * fg * fl,
      (1. - fZ) * fT * (1. - fg) * (1. - fl),
      (1. - fZ) * fT * (1. - fg) * fl,
      (1. - fZ) * fT * fg * (1. - fl),
      (1. - fZ) * fT * fg * fl,
      fZ * (1. - fT) * (1. - fg) * (1. - fl),
      fZ * (1. - fT) * (1. - fg) * fl,
      fZ * (1. - fT) * fg * (1. - fl),
      fZ * (1. - fT) * fg * fl,
      fZ * fT * (1. - fg) * (1. - fl),
      fZ * fT * (1. - fg) * fl,
      fZ * fT * fg * (1. - fl),
      fZ * fT * fg * fl,
  };
  const double Fvals[16] = {
      F[{il, iZ, iT, ig}],
      F[{il + 1, iZ, iT, ig}],
      F[{il, iZ, iT, ig + 1}],
      F[{
          il + 1,
          iZ,
          iT,
          ig + 1,
      }],
      F[{il, iZ, iT + 1, ig}],
      F[{il + 1, iZ, iT + 1, ig}],
      F[{il, iZ, iT + 1, ig + 1}],
      F[{il + 1, iZ, iT + 1, ig + 1}],
      F[{il, iZ + 1, iT, ig}],
      F[{il + 1, iZ + 1, iT, ig}],
      F[{il, iZ + 1, iT, ig + 1}],
      F[{il + 1, iZ + 1, iT, ig + 1}],
      F[{il, iZ + 1, iT + 1, ig}],
      F[{il + 1, iZ + 1, iT + 1, ig}],
      F[{il, iZ + 1, iT + 1, ig + 1}],
      F[{il + 1, iZ + 1, iT + 1, ig + 1}],
  };
  double Fval = 0.;
  for (uint_fast8_t j = 0; j < 16; ++j) {
    cmac_assert(factor[j] >= 0.);
    cmac_assert(factor[j] <= 1.);
    if (factor[j] > 0. && Fvals[j] > 0.) {
      Fval += factor[j] * std::log(Fvals[j]);
    } else {
      if (Fvals[j] <= 0.) {
        return 0.;
      }
    }
  }
  return std::exp(Fval);
}

/**
 * @brief Constructor.
 *
 * Reads in the correct model spectrum and resamples it on the 1000 bin internal
 * frequency grid.
 *
 * @param temperature Effective temperature of the star (in K).
 * @param surface_gravity Surface gravity of the star (in m s^-2).
 * @param metallicity Metallicity of the star.
 * @param log Log to write logging info to.
 */
CastelliKuruczPhotonSourceSpectrum::CastelliKuruczPhotonSourceSpectrum(
    const double temperature, const double surface_gravity,
    const double metallicity, Log *log) {

  // check if the input parameters are valid
  if (!is_valid(temperature, surface_gravity, metallicity)) {
    if (log) {
      log->write_error(
          "Invalid parameters for Castelli-Kurucz spectrum: Teff = ",
          temperature, " K, g = ", surface_gravity,
          " m s^-2, Z = ", metallicity, "!");
    }
    cmac_error("Invalid parameters for Castelli-Kurucz spectrum: Teff = %g K, "
               "g = %g m s^-2, Z = %g!",
               temperature, surface_gravity, metallicity);
  }

  // open the data file and read the data
  HDF5Tools::HDF5File file = HDF5Tools::open_file(CASTELLIKURUCZDATALOCATION,
                                                  HDF5Tools::HDF5FILEMODE_READ);

  // this is the order of the axes in the HDF5 file
  const std::vector< double > lambda =
      HDF5Tools::read_dataset< double >(file, "lambda");
  const std::vector< double > Z = HDF5Tools::read_dataset< double >(file, "Z");
  const std::vector< double > Teff =
      HDF5Tools::read_dataset< double >(file, "Teff");
  const std::vector< double > g = HDF5Tools::read_dataset< double >(file, "g");
  const HDF5Tools::HDF5DataBlock< double, 4 > Flambda =
      HDF5Tools::read_dataset< double, 4 >(file, "Flambda");

  HDF5Tools::close_file(file);

  // locate the parameter values in their respective arrays
  const uint_fast32_t iZ = Utilities::locate(metallicity, &Z[0], Z.size());
  const uint_fast32_t iT =
      Utilities::locate(temperature, &Teff[0], Teff.size());
  const uint_fast32_t ig = Utilities::locate(surface_gravity, &g[0], g.size());
  // compute the (log) fractions of the values in the brackets
  const double Zmin = std::log(Z[iZ]);
  const double Zmax = std::log(Z[iZ + 1]);
  const double Zval = std::log(metallicity);
  cmac_assert(Zval >= Zmin);
  cmac_assert(Zval <= Zmax);
  const double fZ = (Zval - Zmin) / (Zmax - Zmin);
  const double Tmin = std::log(Teff[iT]);
  const double Tmax = std::log(Teff[iT + 1]);
  const double Tval = std::log(temperature);
  cmac_assert(Tval >= Tmin);
  cmac_assert(Tval <= Tmax);
  const double fT = (Tval - Tmin) / (Tmax - Tmin);
  const double gmin = std::log(g[ig]);
  const double gmax = std::log(g[ig + 1]);
  const double gval = std::log(surface_gravity);
  cmac_assert(gval >= gmin);
  cmac_assert(gval <= gmax);
  const double fg = (gval - gmin) / (gmax - gmin);

  if (log) {
    log->write_info("Bracketing intervals:");
    log->write_info("Z: ", Z[iZ], " <= ", metallicity, " <= ", Z[iZ + 1], " (",
                    iZ, " - ", iZ + 1, ")");
    log->write_info("T: ", Teff[iT], " <= ", temperature, " <= ", Teff[iT + 1],
                    " K (", iT, " - ", iT + 1, ")");
    log->write_info("g: ", g[ig], " <= ", surface_gravity, " <= ", g[ig + 1],
                    " m s^-2 (", ig, " - ", ig + 1, ")");
  }

  // allocate memory for the data tables
  _frequencies.resize(CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ, 0.);
  _cumulative_distribution.resize(CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ,
                                  0.);

  // construct the frequency bins
  // 13.6 eV in Hz
  const double min_frequency = 3.289e15;
  // 54.4 eV
  const double max_frequency = 4. * min_frequency;
  for (uint_fast32_t i = 0; i < CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ;
       ++i) {
    _frequencies[i] =
        min_frequency + i * (max_frequency - min_frequency) /
                            (CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ - 1.);
  }

  // create the spectrum in the bin range of interest
  _cumulative_distribution[0] = 0.;
  for (uint_fast32_t i = 1; i < CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ;
       ++i) {
    // get the bracketing frequencies
    const double nu1 = _frequencies[i - 1];
    const double nu2 = _frequencies[i];
    // convert the frequencies to wavelengths
    const double l1 =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED) /
        nu1;
    const double l2 =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_LIGHTSPEED) /
        nu2;

    // locate the wavelengths in the wavelength axis
    const uint_fast32_t il1 = Utilities::locate(l1, &lambda[0], lambda.size());
    const uint_fast32_t il2 = Utilities::locate(l2, &lambda[0], lambda.size());

    // get the corresponding lambda brackets and factors
    const double l1min = std::log(lambda[il1]);
    const double l1max = std::log(lambda[il1 + 1]);
    const double l1val = std::log(l1);
    cmac_assert(l1val >= l1min);
    cmac_assert(l1val <= l1max);
    const double fl1 = (l1val - l1min) / (l1max - l1min);
    const double l2min = std::log(lambda[il2]);
    const double l2max = std::log(lambda[il2 + 1]);
    const double l2val = std::log(l2);
    cmac_assert(l2val >= l2min);
    cmac_assert(l2val <= l2max);
    const double fl2 = (l2val - l2min) / (l2max - l2min);

    // get the corresponding Flambda values and convert from flux per unit
    // wavelength to flux per unit frequency
    const double Flambda1 =
        interpolate(fZ, fT, fg, fl1, iZ, iT, ig, il1, Flambda) * l1 / nu1;
    const double Flambda2 =
        interpolate(fZ, fT, fg, fl2, iZ, iT, ig, il2, Flambda) * l2 / nu2;

    cmac_assert(Flambda1 == Flambda1);
    cmac_assert(Flambda2 == Flambda2);

    // first order quadrature to get the value of the spectrum in this
    // particular bin
    // note that we need to divide by the frequency to keep the units consistent
    // (check against PlanckPhotonSourceSpectrum if unsure)
    _cumulative_distribution[i] =
        0.5 * (Flambda1 / nu1 + Flambda2 / nu2) * (nu2 - nu1);
  }

  // _cumulative_distribution now contains the actual ionizing spectrum
  // make it cumulative (and at the same time get the total ionizing
  // luminosity)
  for (uint_fast32_t i = 1; i < CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ;
       ++i) {
    _cumulative_distribution[i] += _cumulative_distribution[i - 1];
  }

  // the total ionizing luminosity is the last element of
  // _cumulative_distribution (using a zeroth order quadrature)
  // its value is in J Hz^-1 s^-1 m^-2 sr^-1
  // we convert to s^-1 m^-2 sr^-1 by dividing by Planck's constant
  // (in J s = J Hz^-1)
  _total_flux =
      _cumulative_distribution[CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ - 1] /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK);
  // we integrate out over all solid angles
  _total_flux *= 4. * M_PI;
  // normalize the spectrum
  for (uint_fast32_t i = 0; i < CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ;
       ++i) {
    _cumulative_distribution[i] /=
        _cumulative_distribution[CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ -
                                 1];
  }

  if (log) {
    log->write_status(
        "Constructed CastelliKuruczPhotonSourceSpectrum with temperature ",
        temperature, " K, surface gravity ", surface_gravity,
        " m s^-2, metallicity ", metallicity, ", and with total ionizing flux ",
        _total_flux, " m^-2 s^-1.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - temperature: Temperature of the star (default: 4.e4 K)
 *  - surface gravity: Surface gravity of the star (default: 317. m s^-2)
 *  - metallicity: Metallicity of the star (default: 0.02)
 *
 * @param role Role the spectrum will fulfil in the simulation. Parameters are
 * read from the corresponding block in the parameter file.
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
CastelliKuruczPhotonSourceSpectrum::CastelliKuruczPhotonSourceSpectrum(
    std::string role, ParameterFile &params, Log *log)
    : CastelliKuruczPhotonSourceSpectrum(
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              role + ":temperature", "4.e4 K"),
          params.get_physical_value< QUANTITY_ACCELERATION >(
              role + ":surface gravity", "317. m s^-2"),
          params.get_value< double >(role + ":metallicity", 0.02), log) {}

/**
 * @brief Check if the given parameter combination points to a valid spectrum.
 *
 * @param temperature Temperature of the star (in K).
 * @param surface_gravity Surface gravity of the star (in m s^-2).
 * @param metallicity Metallicity of the star.
 * @return True if the given parameter combination corresponds to a valid
 * spectrum.
 */
bool CastelliKuruczPhotonSourceSpectrum::is_valid(const double temperature,
                                                  const double surface_gravity,
                                                  const double metallicity) {

  if (temperature < 3500. || temperature > 50000.) {
    return false;
  }
  if (surface_gravity < 0.01 || surface_gravity > 1.e3) {
    return false;
  }
  if (metallicity < 6.325e-5 || metallicity > 6.325e-2) {
    return false;
  }
  const double logg = std::log10(surface_gravity * 100.);
  if (temperature <= 6000.) {
    return true;
  } else if (temperature <= 7500.) {
    return logg >= 0.5;
  } else if (temperature <= 8250.) {
    return logg >= 1.;
  } else if (temperature <= 9000.) {
    return logg >= 1.5;
  } else if (temperature <= 11750.) {
    return logg >= 2.;
  } else if (temperature <= 19000.) {
    return logg >= 2.5;
  } else if (temperature <= 26000.) {
    return logg >= 3.;
  } else if (temperature <= 31000.) {
    return logg >= 3.5;
  } else if (temperature <= 39000.) {
    return logg >= 4.;
  } else if (temperature <= 49000.) {
    return logg >= 4.5;
  } else if (temperature <= 50000.) {
    return logg == 5.;
  } else {
    return false;
  }
}

/**
 * @brief Get a random frequency from a stellar model spectrum.
 *
 * We draw a random uniform number in the range [0, 1] and find the bin in the
 * tabulated cumulative distribution that contains that number. The sampled
 * frequency is then the linearly interpolated value corresponding to that bin.
 *
 * @param random_generator RandomGenerator to use.
 * @param temperature Not used for this spectrum.
 * @return Random frequency (in Hz).
 */
double CastelliKuruczPhotonSourceSpectrum::get_random_frequency(
    RandomGenerator &random_generator, double temperature) const {

  const double x = random_generator.get_uniform_random_double();
  const uint_fast32_t inu =
      Utilities::locate(x, _cumulative_distribution.data(),
                        CASTELLIKURUCZPHOTONSOURCESPECTRUM_NUMFREQ);
  const double frequency =
      _frequencies[inu] +
      (_frequencies[inu + 1] - _frequencies[inu]) *
          (x - _cumulative_distribution[inu]) /
          (_cumulative_distribution[inu + 1] - _cumulative_distribution[inu]);
  return frequency;
}

/**
 * @brief Get the total ionizing flux of the spectrum.
 *
 * The total ionizing flux depends on the model spectrum that is used and is
 * computed in the constructor, after the spectrum has been regridded.
 *
 * @return Total ionizing flux (in m^-2 s^-1).
 */
double CastelliKuruczPhotonSourceSpectrum::get_total_flux() const {
  return _total_flux;
}
