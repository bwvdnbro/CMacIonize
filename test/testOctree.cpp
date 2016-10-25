/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testOctree.cpp
 *
 * @brief Unit test for the Octree class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "Octree.hpp"
#include "Utilities.hpp"
#include <vector>

/**
 * @brief Unit test for the Octree class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // first test the Box-CoordinateVector distance function
  Box box(CoordinateVector<>(), CoordinateVector<>(1.));
  {
    // point outside and smaller coordinate
    CoordinateVector<> caseA(-0.5, 0., 0.);
    assert_condition(box.get_distance(caseA) == 0.5);
  }
  {
    // point inside
    CoordinateVector<> caseB(0., 0.5, 0.);
    assert_condition(box.get_distance(caseB) == 0.);
  }
  {
    // point outside and larger coordinate
    CoordinateVector<> caseC(0., 0., 1.5);
    assert_condition(box.get_distance(caseC) == 0.5);
  }

  unsigned int numpos = 100;
  std::vector< CoordinateVector<> > positions(numpos);
  std::vector< double > hs(numpos);
  for (unsigned int i = 0; i < numpos; ++i) {
    positions[i][0] = Utilities::random_double();
    positions[i][1] = Utilities::random_double();
    positions[i][2] = Utilities::random_double();
    hs[i] = 0.5 * Utilities::random_double();
  }

  Octree tree(positions, box, false);
  tree.set_auxiliaries(hs, Octree::max< double >);

  CoordinateVector<> centre(0.5);
  std::vector< unsigned int > ngbs_brute_force;
  for (unsigned int i = 0; i < numpos; ++i) {
    double r = (positions[i] - centre).norm();
    if (r < hs[i]) {
      ngbs_brute_force.push_back(i);
    }
  }
  status("Number of ngbs (brute force): %lu.", ngbs_brute_force.size());

  std::vector< unsigned int > ngbs_tree = tree.get_ngbs(centre);
  status("Number of ngbs (tree): %lu.", ngbs_tree.size());

  assert_condition(ngbs_brute_force.size() == ngbs_tree.size());
  std::sort(ngbs_brute_force.begin(), ngbs_brute_force.end());
  std::sort(ngbs_tree.begin(), ngbs_tree.end());
  for (unsigned int i = 0; i < ngbs_tree.size(); ++i) {
    assert_condition(ngbs_brute_force[i] == ngbs_tree[i]);
  }

  return 0;
}
