#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CMacIonize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMacIonize is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file write_FLASHtest.py
#
# @brief Write "FLASHtest.hdf5", a minimal FLASH snapshot that can be used to
# test the FLASHSnapshotDensityFunction.
#
# We set up a simple AMR grid of 64x32x32 cells, with a refinement level of
# 128x64x64 in the central x slab.
#
# The density in the entire box is given by 1.+x+y+z.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import numpy as np
import h5py

##
# @brief Recursive AMR cell representation.
##
class AMRGridCell:
  ##
  # @brief Constructor.
  #
  # @param self Reference to this object.
  # @param box Geometrical extents of the cell.
  # @param level Refinement level of the cell.
  # @param density Density in the cell.
  ##
  def __init__(self, box, level, density):
    self._box = box
    self._level = level
    self._density = density
    self._node_type = 1
    self._children = np.array([None, None, None, None, None, None, None, None])
    self._next = None

  ##
  # @brief Add a child to this cell.
  #
  # @param self Reference to this object.
  # @param cell Cell to add.
  # @param key Child key where the cell needs to be added.
  ##
  def add_child(self, cell, key):
    self._children[key] = cell

  ##
  # @brief Get the node type for this cell.
  #
  # The node type is equal to the number of hierarchical generations below the
  # cell: if a cell has children and grand-children, its node type is 3. If it
  # has no children, its node type is 1.
  #
  # @param self Reference to this object.
  # @return Node type of this cell.
  ##
  def get_node_type(self):
    maxchild = 1
    for i in range(8):
      if self._children[i]:
        maxchild = max(maxchild, self._children[i].get_node_type() + 1)
    self._node_type = maxchild
    return self._node_type

  ##
  # @brief Get the requested child.
  #
  # @param self Reference to this object.
  # @param key Key of a child.
  # @return Reference to that child.
  ##
  def get_child(self, key):
    return self._children[key]

  ##
  # @brief Get the number of cells in this cell (including the self itself).
  #
  # @param self Reference to this object.
  # @return Number of cells in this cell (including self).
  ##
  def get_num_cells(self):
    num = 1
    for i in range(8):
      if self._children[i]:
        num += self._children[i].get_num_cells()
    return num

  ##
  # @brief Link the cell together with its children and siblings in a long
  # chain.
  #
  # The goal is to be able to traverse all cells by just calling cell._next on
  # each cell, until cell._next == None.
  #
  # @param self Reference to this object.
  # @param next Reference to the next cell in the chain (or None if there is
  # none).
  ##
  def set_order(self, next):
    # find first child
    i = 0
    while i < 8 and not self._children[i]:
      i += 1
    if i == 8:
      # this cell has no children. The next cell in the chain is the argument.
      self._next = next
    else:
      # this cell has children. The next cell in the chain is its first child.
      child = self._children[i]
      self._next = child
      i += 1
      # find next child
      while i < 8:
        while i < 8 and not self._children[i]:
          i += 1
        if i < 8:
          # this child has a sibling, which is its next cell in the chain
          child.set_order(self._children[i])
          child = self._children[i]
          i += 1
      if i == 8:
        # this was the last child, the next cell in the chain is this cell's
        # sibling
        child.set_order(next)

##
# @brief Density function.
#
# @param x X coordinate of a position.
# @param y Y coordinate of a position.
# @param z Z coordinate of a position.
# @return Density at that position.
##
def density(x, y, z):
  return 1. + x + y + z

# open the file for writing
file = h5py.File("FLASHtest.hdf5", 'w')

# write real runtime parameters
comp_type = np.dtype([("name", np.str_, 80), ("value", 'd')])
dataset = file.create_dataset("real runtime parameters", (6,), comp_type)
data = np.array([("xmin", 0.), ("xmax", 2.), ("ymin", 0.), ("ymax", 1.),
                 ("zmin", 0.), ("zmax", 1.)], dtype = comp_type)
dataset[...] = data

# write integer runtime parameters
comp_type = np.dtype([("name", np.str_, 80), ("value", 'i')])
dataset = file.create_dataset("integer runtime parameters", (3,), comp_type)
data = np.array([("nblockx", 2), ("nblocky", 1), ("nblockz", 1)],
                dtype = comp_type)
dataset[...] = data

# make grid
grid = np.array([AMRGridCell(np.array([[0., 0., 0.], [1., 1., 1.]]),
                             1, density(0.5, 0.5, 0.5)),
                 AMRGridCell(np.array([[1., 0., 0.], [1., 1., 1.]]),
                             1, density(1.5, 0.5, 0.5))])

for iblock in range(2):
  block = grid[iblock]
  for ix in range(2):
    for iy in range(2):
      for iz in range(2):
        box = np.array(block._box)
        box[1] *= 0.5
        box[0][0] += ix*box[1][0]
        box[0][1] += iy*box[1][1]
        box[0][2] += iz*box[1][2]
        x = box[0][0] + 0.5*box[1][0]
        y = box[0][1] + 0.5*box[1][1]
        z = box[0][2] + 0.5*box[1][2]
        block.add_child(AMRGridCell(box, 2, density(x, y, z)), 4*ix+2*iy+iz)
        # we make a refined region near the centre of the box
        if ix != iblock:
          child = block.get_child(4*ix+2*iy+iz)
          for iix in range(2):
            for iiy in range(2):
              for iiz in range(2):
                box = np.array(child._box)
                box[1] *= 0.5
                box[0][0] += iix*box[1][0]
                box[0][1] += iiy*box[1][1]
                box[0][2] += iiz*box[1][2]
                x = box[0][0] + 0.5*box[1][0]
                y = box[0][1] + 0.5*box[1][1]
                z = box[0][2] + 0.5*box[1][2]
                child.add_child(AMRGridCell(box, 3, density(x, y, z)),
                                4*iix+2*iiy+iiz)
# set up the chain for easy grid traversal
grid[0].set_order(None)
grid[1].set_order(None)

# gather grid data
numblocks = grid[0].get_num_cells() + grid[1].get_num_cells()
bbox = np.zeros((numblocks, 3, 2), dtype='f')
rho = np.zeros((numblocks, 8, 8, 8), dtype='f')
levels = np.zeros((numblocks, 1), dtype='i')
ntypes = np.zeros((numblocks, 1), dtype='i')
index = 0
for iblock in range(2):
  child = grid[iblock]
  while child:
    bbox[index][0][0] = child._box[0][0]
    bbox[index][1][0] = child._box[0][1]
    bbox[index][2][0] = child._box[0][2]
    bbox[index][0][1] = child._box[0][0] + child._box[1][0]
    bbox[index][1][1] = child._box[0][1] + child._box[1][1]
    bbox[index][2][1] = child._box[0][2] + child._box[1][2]
    for ix in range(8):
      for iy in range(8):
        for iz in range(8):
          x = child._box[0][0] + 0.125*(ix+0.5)*child._box[1][0]
          y = child._box[0][1] + 0.125*(iy+0.5)*child._box[1][1]
          z = child._box[0][2] + 0.125*(iz+0.5)*child._box[1][2]
          rho[index][iz][iy][ix] = density(x, y, z)
    levels[index] = child._level
    ntypes[index] = child.get_node_type()
    child = child._next
    index += 1

# write datasets to file
file.create_dataset("bounding box", data = bbox, dtype = 'f')
file.create_dataset("dens", data = rho, dtype = 'f')
file.create_dataset("refine level", data = levels, dtype = 'i')
file.create_dataset("node type", data = ntypes, dtype = 'i')
