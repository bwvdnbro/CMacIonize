#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file write_parameterfile.py
#
# @brief Script that can be used to read and write CMacIonize parameter files.
#
# This script is meant to be imported from another script and makes available
# the functions `read_parameterfile()` and `write_parameterfile` that can be
# used to respectively parse the contents of a parameter file into a Python
# dictionary, or write out the contents of a Python dictionary to a parameter
# file. When used on its own, the script takes two required arguments, `input`
# and `output`. In this case, the parameter file corresponding to `input` is
# read in and then written out to `output`. This is only useful as a test.
#
# @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
##

import sys
import re
import time
import datetime


# read the parameter file with the given name and return a dictionary containing
# its contents
def read_parameterfile(filename):
    file = open(filename, "r")

    params = {}
    indent = 0
    groupname = []
    lines = 0
    for line in file.readlines():
        if line[0] == "#":
            continue
        m = re.match(r"(?P<indent>[ ]*)(?P<groupname>.*?):([ ]?#.*?)?$", line)
        if m:
            lines += 1
            current_indent = len(m.group("indent"))
            if current_indent > indent:
                indent = current_indent
            if current_indent < indent:
                indent = current_indent
                groupname.pop()
            groupname.append(m.group("groupname"))
        m = re.match(
            r"(?P<indent>[ ]*)(?P<key>.*?):[ ]*" "(?P<value>.+?)[#\n]", line
        )
        if m:
            lines += 1
            current_indent = len(m.group("indent"))
            if current_indent > indent:
                indent = current_indent
            if current_indent < indent:
                indent = current_indent
                groupname.pop()
            groupname.append(m.group("key"))
            full_name = ":".join(groupname)
            groupname.pop()
            params[full_name] = m.group("value").strip()

    return params


# write the given parameter dictionary to a new parameter file with the given
# name
def write_parameterfile(filename, params):
    file = open(filename, "w")

    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime(
        "%d/%m/%Y, %H:%M:%S"
    )
    file.write(
        "# file automatically written by write_parameterfile.py on "
        "{timestamp}\n".format(timestamp=timestamp)
    )

    groupname = []
    stream = ""
    for key in sorted(params):
        keygroups = key.split(":")
        keyname = keygroups[-1]
        keygroups.pop()
        indent = ""
        if len(keygroups) > len(groupname):
            i = 0
            while i < len(groupname) and groupname[i] == keygroups[i]:
                i += 1
            for j in range(i):
                indent += "  "
            for j in range(i, len(groupname)):
                groupname.pop()
            for j in range(i, len(keygroups)):
                groupname.append(keygroups[j])
                stream += indent + keygroups[j] + ":\n"
                indent += "  "
        else:
            while len(keygroups) < len(groupname):
                groupname.pop()
            i = 0
            while i < len(keygroups) and groupname[i] == keygroups[i]:
                i += 1
            for j in range(i):
                indent += "  "
            for j in range(i, len(keygroups)):
                groupname.pop()
            for j in range(i, len(keygroups)):
                groupname.append(keygroups[j])
                stream += indent + keygroups[j] + ":\n"
                indent += "  "

        stream += indent + keyname + ": " + params[key] + "\n"

    file.write(stream)


# unit test the read and write functions
if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser(
        "Copy the given input parameter file into the given output by reading it into a Python dictionary and writing out that dictionary again."
    )
    argparser.add_argument("--input", "-i", action="store", required=True)
    argparser.add_argument("--output", "-o", action="store", required=True)
    args = argparser.parse_args()

    params = read_parameterfile(args.input)
    write_parameterfile(args.output, params)
