#! /usr/bin/python

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file hilbertkeygenerator_spacefilling_curve.py
#
# @brief Script that generates a TikZ TeX file that shows a two level 3D space-
# filling Hilbert curve with the corresponding level keys.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

##
# @brief Get the three bit components of the given one level key.
#
# @param key 3-bit key.
# @return String representation of the three bits of the key.
##
def get_bits(key):
    highbit = key >> 2
    midbit = (key - (highbit<<2)) >> 1
    lowbit = key - (highbit<<2) - (midbit<<1)
    return "{h}{m}{l}".format(h = highbit, m = midbit, l = lowbit)

##
# @brief Split the given 6-bit key into two 3-bit keys.
#
# @param key 6-bit key.
# @return Two 3-bit key components.
##
def split_key(key):
    nodekeyint = key >> 3;
    smallkeyint = key - (nodekeyint<<3)
    
    nodekey = get_bits(nodekeyint)
    smallkey = get_bits(smallkeyint)
    
    return nodekey, smallkey

##
# @brief Compute the 2D (x, y) shift in projection for the given 3D position.
#
# @param position 3D position.
# @return x and y shift in the 2D projection space.
##
def get_shifts(position):
    x = position["x"]
    y = position["y"]
    z = position["z"]
    
    xs = z*4.5
    ys = y*4.5
    
    xs -= 1.2*(3-x)
    ys -= 1.8*(3-x)
    
    xshift = "{xs}cm".format(xs = xs)
    yshift = "{ys}cm".format(ys = ys)
    
    return xshift, yshift

# read positions from the data file and link keys to them (the positions in the
# data file are in key order)
positions = {}
file = open("hilbertkeygenerator_spacefilling_curve.txt", "r")
key = 0
for line in file.readlines():
    cols = line.split()
    positions[key] = {}
    positions[key]["x"] = int(cols[0])
    positions[key]["y"] = int(cols[1])
    positions[key]["z"] = int(cols[2])
    key += 1

# we do our own z-sort by ordering the coordinates in 4 depth planes
back = []
front1 = []
front2 = []
front3 = []
# traverse the keys and make labels for them
# add them to the right depth plane
for key in positions:
    nodekey, smallkey = split_key(key)
    label = "{{\\color{{gray}}{nodekey}}}{smallkey}".format(nodekey = nodekey,
                                                            smallkey = smallkey)
    nodename = "{nodekey}{smallkey}".format(nodekey = nodekey,
                                            smallkey = smallkey)
    positions[key]["label"] = label
    positions[key]["name"] = nodename
    xshift, yshift = get_shifts(positions[key])
    positions[key]["xshift"] = xshift
    positions[key]["yshift"] = yshift
    if positions[key]["x"] == 3:
        back.append(key)
    if positions[key]["x"] == 2:
        front1.append(key)
    if positions[key]["x"] == 1:
        front2.append(key)
    if positions[key]["x"] == 0:
        front3.append(key)

# sort the depth planes on key
back = sorted(back)
front1 = sorted(front1)
front2 = sorted(front2)
front3 = sorted(front3)

# start writing the TikZ TeX file
file = open("hilbertkeygenerator_spacefilling_curve.tex", "w")
file.write(r"""\documentclass[convert={density=72,outext=.png}]{standalone}
\usepackage{tikz}
\usetikzlibrary{shapes.geometric, arrows}
\usepackage{amsmath}

\tikzstyle{node} = [rectangle, text centered, align=center]
\tikzstyle{stealthnode} = [rectangle, text centered, align=center, color=gray]
\tikzstyle{arrow} = [thick,>=stealth]
\tikzstyle{stealtharrow} = [thick,>=stealth, color=gray, dashed]

\tikzstyle{back1} = [xshift=0.8cm, yshift=1.2cm]
\tikzstyle{back2} = [xshift=1.6cm, yshift=2.4cm]
\tikzstyle{back3} = [xshift=2.4cm, yshift=3.6cm]
\tikzstyle{front1} = [xshift=-1.2cm, yshift=-1.8cm]

\newcommand{\drawLinewithBG}[3]
{
    \draw [white,line width=3pt, opacity=1.0]  (#1) -- (#2);
    \draw [#3] (#1) -- (#2);
}

\begin{document}

\begin{tikzpicture}[node distance=4.5cm]
""")

file.write("% background\n")
for key in back:
    file.write("\\node ({nodename}) ".format(nodename = positions[key]["name"]))
    file.write("[node, xshift={xshift}".format(xshift = positions[key]["xshift"]))
    file.write(", yshift={yshift}] ".format(yshift = positions[key]["yshift"]))
    file.write("{{{label}}};\n".format(label = positions[key]["label"]))

for ikey in range(len(back)-1):
    if back[ikey+1]-back[ikey] == 1:
        file.write("\\draw [arrow] ({n1})".format(n1 = positions[back[ikey]]["name"]))
        file.write(" -- ({n2});\n".format(n2 = positions[back[ikey+1]]["name"]))

file.write("% front1\n")
for key in front1:
    file.write("\\node ({nodename}) ".format(nodename = positions[key]["name"]))
    file.write("[node, xshift={xshift}".format(xshift = positions[key]["xshift"]))
    file.write(", yshift={yshift}] ".format(yshift = positions[key]["yshift"]))
    file.write("{{{label}}};\n".format(label = positions[key]["label"]))

file.write(r"""\path (010000) -- (010101) node[stealthnode, midway] (010) {010};
\path (011000) -- (011101) node[stealthnode, midway] (011) {011};
\path (100000) -- (100101) node[stealthnode, midway] (100) {100};
\path (101000) -- (101101) node[stealthnode, midway] (101) {101};
\drawLinewithBG{010}{011}{stealtharrow};
\drawLinewithBG{011}{100}{stealtharrow};
\drawLinewithBG{100}{101}{stealtharrow};
""")

# write the depth planes to the file
for ikey in range(len(front1)-1):
    if front1[ikey+1]-front1[ikey] == 1:
        file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front1[ikey]]["name"]))
        file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front1[ikey+1]]["name"]))

for ikey in range(len(front1)):
    if ikey < len(front1)-1 and not front1[ikey]+1 in front1:
        file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front1[ikey]]["name"]))
        file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front1[ikey]+1]["name"]))
    if ikey > 0 and not front1[ikey]-1 in front1:
        file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front1[ikey]]["name"]))
        file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front1[ikey]-1]["name"]))

file.write("% front2\n")
for key in front2:
    file.write("\\node ({nodename}) ".format(nodename = positions[key]["name"]))
    file.write("[node, xshift={xshift}".format(xshift = positions[key]["xshift"]))
    file.write(", yshift={yshift}] ".format(yshift = positions[key]["yshift"]))
    file.write("{{{label}}};\n".format(label = positions[key]["label"]))

file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front1[0]]["name"]))
file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front1[0]-1]["name"]))
file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front1[-1]]["name"]))
file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front1[-1]+1]["name"]))

for ikey in range(len(front2)-1):
    if front2[ikey+1]-front2[ikey] == 1:
        file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front2[ikey]]["name"]))
        file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front2[ikey+1]]["name"]))

file.write("% front3\n")
for key in front3:
    file.write("\\node ({nodename}) ".format(nodename = positions[key]["name"]))
    file.write("[node, xshift={xshift}".format(xshift = positions[key]["xshift"]))
    file.write(", yshift={yshift}] ".format(yshift = positions[key]["yshift"]))
    file.write("{{{label}}};\n".format(label = positions[key]["label"]))

file.write(r"""\path (000000) -- (000101) node[stealthnode, midway] (000) {000};
\path (001000) -- (001101) node[stealthnode, midway] (001) {001};
\path (110000) -- (110101) node[stealthnode, midway] (110) {110};
\path (111000) -- (111101) node[stealthnode, midway] (111) {111};
\drawLinewithBG{001}{010}{stealtharrow};
\drawLinewithBG{101}{110}{stealtharrow};
\drawLinewithBG{000}{001}{stealtharrow};
\drawLinewithBG{110}{111}{stealtharrow};
\drawLinewithBG{001100}{001101}{arrow};
\drawLinewithBG{110010}{110011}{arrow};
""")

for ikey in range(len(front3)-1):
    if front3[ikey+1]-front3[ikey] == 1:
        file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front3[ikey]]["name"]))
        file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front3[ikey+1]]["name"]))

for ikey in range(len(front3)):
    if ikey < len(front3)-1 and not front3[ikey]+1 in front3:
        file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front3[ikey]]["name"]))
        file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front3[ikey]+1]["name"]))
    if ikey > 0 and not front3[ikey]-1 in front3:
        file.write("\\drawLinewithBG{{{n1}}}".format(n1 = positions[front3[ikey]]["name"]))
        file.write("{{{n2}}}{{arrow}};\n".format(n2 = positions[front3[ikey]-1]["name"]))

file.write(r"""\end{tikzpicture}

\end{document}""")
