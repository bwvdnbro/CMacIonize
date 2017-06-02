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
 * @file NewVoronoiCell.cpp
 *
 * @brief NewVoronoiCell implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "NewVoronoiCell.hpp"
#include "Error.hpp"

/**
 * @brief Constructor.
 *
 * @param generator Index of the generator of the cell.
 */
NewVoronoiCell::NewVoronoiCell(unsigned int generator) {
  _ngbs.resize(7);
  _ngbs[0] = generator;
  _ngbs[1] = NEWVORONOICELL_BOX_LEFT;
  _ngbs[2] = NEWVORONOICELL_BOX_RIGHT;
  _ngbs[3] = NEWVORONOICELL_BOX_FRONT;
  _ngbs[4] = NEWVORONOICELL_BOX_BACK;
  _ngbs[5] = NEWVORONOICELL_BOX_BOTTOM;
  _ngbs[6] = NEWVORONOICELL_BOX_TOP;

  _tetrahedra.resize(8);
  _tetrahedra[0] = VoronoiTetrahedron(0, 1, 3, 6);
  _tetrahedra[1] = VoronoiTetrahedron(0, 4, 1, 6);
  _tetrahedra[2] = VoronoiTetrahedron(0, 2, 4, 6);
  _tetrahedra[3] = VoronoiTetrahedron(0, 3, 2, 6);
  _tetrahedra[4] = VoronoiTetrahedron(0, 1, 4, 5);
  _tetrahedra[5] = VoronoiTetrahedron(0, 4, 2, 5);
  _tetrahedra[6] = VoronoiTetrahedron(0, 2, 3, 5);
  _tetrahedra[7] = VoronoiTetrahedron(0, 3, 1, 5);

  _connections.resize(7);
  _connections[0].resize(8);
  _connections[0][0] = 0;
  _connections[0][1] = 1;
  _connections[0][2] = 2;
  _connections[0][3] = 3;
  _connections[0][4] = 4;
  _connections[0][5] = 5;
  _connections[0][6] = 6;
  _connections[0][7] = 7;

  // we order the tetrahedra counterclockwise around the axis, when looking
  // along the axis from outside the cell towards the cell generator
  _connections[1].resize(4);
  _connections[1][0] = 0;
  _connections[1][1] = 7;
  _connections[1][2] = 4;
  _connections[1][3] = 1;

  _connections[2].resize(4);
  _connections[2][0] = 2;
  _connections[2][1] = 5;
  _connections[2][2] = 6;
  _connections[2][3] = 3;

  _connections[3].resize(4);
  _connections[3][0] = 0;
  _connections[3][1] = 3;
  _connections[3][2] = 6;
  _connections[3][3] = 7;

  _connections[4].resize(4);
  _connections[4][0] = 1;
  _connections[4][1] = 4;
  _connections[4][2] = 5;
  _connections[4][3] = 2;

  _connections[5].resize(4);
  _connections[5][0] = 4;
  _connections[5][1] = 7;
  _connections[5][2] = 6;
  _connections[5][3] = 5;

  _connections[6].resize(4);
  _connections[6][0] = 0;
  _connections[6][1] = 1;
  _connections[6][2] = 2;
  _connections[6][3] = 3;
}

/**
 * @brief Get the volume of the cell.
 *
 * @return Volume of the cell (in m^3).
 */
double NewVoronoiCell::get_volume() const { return _volume; }

/**
 * @brief Get the centroid of the cell.
 *
 * @return Centroid of the cell (in m).
 */
const CoordinateVector<> &NewVoronoiCell::get_centroid() const {
  return _centroid;
}

/**
 * @brief Get the faces of the cell.
 *
 * @return Faces of the cell.
 */
const std::vector< VoronoiFace > &NewVoronoiCell::get_faces() const {
  return _faces;
}

/**
 * @brief Update the cell structure by interacting it with the generator with
 * the given index.
 *
 * @param ngb Index of the neighbouring generator.
 * @param positions Generator positions (integer representation).
 */
void NewVoronoiCell::intersect(
    unsigned int ngb,
    const std::vector< CoordinateVector< unsigned long > > &positions) {
  cmac_error("If you *just* implement this method, the algorithm will work!");
}

/**
 * @brief Get the maximum distance (squared) between the cell generator and an
 * arbitrary other generator that still could change the cell structure.
 *
 * @return Maximum influence radius squared (in m^2).
 */
double NewVoronoiCell::get_max_radius_squared() const {
  cmac_error("This method has not been implemented!");
  return 0;
}

/**
 * @brief Compute geometrical properties of the cell and clean up the cell
 * structure information.
 *
 * @param box Bounding box of the grid (in m).
 * @param positions Positions of the generators (in m).
 */
void NewVoronoiCell::finalize(
    const Box &box, const std::vector< CoordinateVector<> > &positions) {
  std::vector< CoordinateVector<> > real_positions(_ngbs.size());
  for (unsigned int i = 0; i < _ngbs.size(); ++i) {
    if (_ngbs[i] < NEWVORONOICELL_MAX_INDEX) {
      real_positions[i] = positions[_ngbs[i]];
    } else {
      real_positions[i] = positions[_ngbs[0]];
      if (_ngbs[i] == NEWVORONOICELL_BOX_LEFT) {
        real_positions[i][0] = 2. * box.get_anchor().x() - real_positions[i][0];
      } else if (_ngbs[i] == NEWVORONOICELL_BOX_RIGHT) {
        real_positions[i][0] =
            2. * (box.get_anchor().x() + box.get_sides().x()) -
            real_positions[i][0];
      } else if (_ngbs[i] == NEWVORONOICELL_BOX_FRONT) {
        real_positions[i][1] = 2. * box.get_anchor().y() - real_positions[i][1];
      } else if (_ngbs[i] == NEWVORONOICELL_BOX_BACK) {
        real_positions[i][1] =
            2. * (box.get_anchor().y() + box.get_sides().y()) -
            real_positions[i][1];
      } else if (_ngbs[i] == NEWVORONOICELL_BOX_BOTTOM) {
        real_positions[i][2] = 2. * box.get_anchor().z() - real_positions[i][2];
      } else if (_ngbs[i] == NEWVORONOICELL_BOX_TOP) {
        real_positions[i][2] =
            2. * (box.get_anchor().z() + box.get_sides().z()) -
            real_positions[i][2];
      } else {
        cmac_error("Invalid generator index!");
      }
    }
  }

  std::vector< CoordinateVector<> > cell_vertices(_tetrahedra.size() + 1);
  cell_vertices[0] = real_positions[0];
  for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
    cell_vertices[i + 1] =
        _tetrahedra[i].get_midpoint_circumsphere(real_positions);
  }

  // due to the ordering of the connections, this constructs the faces in a
  // counterclockwise direction when looking from outside the cell towards the
  // cell generator, through the face
  _volume = 0.;
  for (unsigned int i = 1; i < _connections.size(); ++i) {
    double area = 0.;
    CoordinateVector<> midpoint;
    for (unsigned int j = 2; j < _connections[i].size(); ++j) {
      const VoronoiTetrahedron tetrahedron(0, _connections[i][0] + 1,
                                           _connections[i][j] + 1,
                                           _connections[i][j - 1] + 1);
      const double tvol = tetrahedron.get_volume(cell_vertices);
      const CoordinateVector<> tcentroid =
          tetrahedron.get_centroid(cell_vertices);
      _volume += tvol;
      _centroid += tvol * tcentroid;

      const CoordinateVector<> r1 = cell_vertices[_connections[i][j] + 1] -
                                    cell_vertices[_connections[i][0] + 1];
      const CoordinateVector<> r2 = cell_vertices[_connections[i][j - 1] + 1] -
                                    cell_vertices[_connections[i][0] + 1];
      const CoordinateVector<> w = CoordinateVector<>::cross_product(r1, r2);
      const double tarea = 0.5 * w.norm();
      const CoordinateVector<> tmidpoint =
          (cell_vertices[_connections[i][0] + 1] +
           cell_vertices[_connections[i][j - 1] + 1] +
           cell_vertices[_connections[i][j] + 1]) /
          3.;
      area += tarea;
      midpoint += tarea * tmidpoint;
    }
    midpoint /= area;
    _faces.push_back(VoronoiFace(area, midpoint, _ngbs[i]));
  }
  _centroid /= _volume;
}
