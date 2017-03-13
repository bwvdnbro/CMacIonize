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
 * @file testPhotonPool.cpp
 *
 * @brief Unit test for the PhotonPool classes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "PhotonBatch.hpp"
#include "PhotonBucket.hpp"
#include "PhotonPool.hpp"

/**
 * @brief Unit test for the PhotonPool classes.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  /// test PhotonBatch
  {
    PhotonBatch batch(100);

    // check if we can add photons to the batch until it is full
    Photon photon(CoordinateVector<>(), CoordinateVector<>(), 0.);
    for (unsigned int i = 0; i < 100; ++i) {
      assert_condition(batch.add_photon(photon) == true);
    }
    assert_condition(batch.add_photon(photon) == false);

    // check if the iterator visits all photons (and set some values)
    int numphoton = 0;
    for (auto it = batch.begin(); it != batch.end(); ++it) {
      it.set_new_index(numphoton);
      ++numphoton;
    }
    assert_condition(numphoton == 100);

    // check if we can read values (and if the order is always the same)
    numphoton = 0;
    for (auto it = batch.begin(); it != batch.end(); ++it) {
      assert_condition(numphoton == it.get_new_index());
      ++numphoton;
    }
  }

  /// test PhotonBucket
  {
    PhotonBucket bucket(100);

    Photon photon(CoordinateVector<>(), CoordinateVector<>(), 0.);
    for (unsigned int i = 0; i < 150; ++i) {
      bucket.add_photon(photon);
    }

    // check if the bucket contains exactly two batches, with respectively
    // 100 and 50 photons
    unsigned int numbucket = 0;
    std::vector< unsigned int > numphoton;
    PhotonBatch *batch = bucket.get_batch();
    while (batch != nullptr) {
      numphoton.push_back(0);
      for (auto it = batch->begin(); it != batch->end(); ++it) {
        ++numphoton[numbucket];
      }
      ++numbucket;
      batch = bucket.get_batch();
    }
    assert_condition(numbucket == 2);
    assert_condition(numphoton[0] == 100);
    assert_condition(numphoton[1] == 50);
  }

  /// test PhotonPool
  {
    PhotonPool pool(1, 100);

    Photon photon(CoordinateVector<>(), CoordinateVector<>(), 0.);
    for (unsigned int i = 0; i < 150; ++i) {
      pool.add_photon(0, photon);
    }

    // check if the pool contains exactly two batches, with respectively
    // 100 and 50 photons
    unsigned int numbucket = 0;
    std::vector< unsigned int > numphoton;
    PhotonBatch *batch = pool.get_batch(0);
    while (batch != nullptr) {
      numphoton.push_back(0);
      for (auto it = batch->begin(); it != batch->end(); ++it) {
        ++numphoton[numbucket];
      }
      ++numbucket;
      batch = pool.get_batch(0);
    }
    assert_condition(numbucket == 2);
    assert_condition(numphoton[0] == 100);
    assert_condition(numphoton[1] == 50);
  }

  return 0;
}
