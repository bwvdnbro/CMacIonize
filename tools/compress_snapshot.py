#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file compress_snapshot.py
#
# @brief Script that checks if the given snapshot file has been compressed and,
# if not, replaces it with a compressed version.
#
# The script takes the name of the snapshot file as an argument (using standard
# GNU parameter syntax). It first checks compression on the snapshot, if the
# snapshot was not compressed, a new snapshot file is created with the same
# contents but maximally gzip-compressed datasets. Once the data are
# successfully copied, the original file is replaced by this new file.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules:
#  - h5py for Python HDF5 manipulation
#  - argparse for GNU style command line argument parsing
#  - os for file operations
import h5py
import argparse
import os

# parse the command line argument
argparser = argparse.ArgumentParser("Compress the given HDF5 snapshot.")
argparser.add_argument("-f", "--file",
                       help = "name of the snapshot file",
                       action = "store", required = True)
args = argparser.parse_args()

# check compression by looping over the datasets in the snapshot
isCompressed = True
file = h5py.File(args.file, "r")
for dataset in file["/PartType0"].keys():
  if not file["/PartType0/" + dataset].compression == "gzip" or \
     not file["/PartType0/" + dataset].compression_opts == 9:
    isCompressed = False

# make a compressed copy of the file if not yet compressed
if not isCompressed:
  # we create a new file with a temporary name
  new_file = h5py.File(args.file + ".copy", "w")
  # we loop over the groups in the original file
  for group in file.keys():
    # create a new group with the same name
    new_group = new_file.create_group(group)
    # if the group contains particle data, we have to copy and compress the
    # datasets. If not, we simply copy the group attributes.
    if group == "PartType0":
      # loop over the datasets
      for ds in file[group].keys():
        # read the original data
        data = file[group + "/" + ds][:]
        # create a new compressed dataset with the same data
        new_group.create_dataset(ds, shape = data.shape, dtype = data.dtype,
                                 data = data, fletcher32 = True, shuffle = True,
                                 compression = "gzip", compression_opts = 9)
    else:
      # loop over the attributes and copy them
      for attr_key in file[group].attrs.keys():
        new_group.attrs[attr_key] = file[group].attrs[attr_key]

  # close the original and new file
  file.close()
  new_file.close()

  # overwrite the old file with the new file
  os.rename(args.file + ".copy", args.file)

  print "Done compressing", args.file

else:
  print args.file, "already compressed!"
