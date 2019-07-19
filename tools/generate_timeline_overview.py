################################################################################
# This file is part of CMacIonize
# Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
# @file generate_timeline_overview.py
#
# @brief Script that generates an interactive web page that allows for a
# detailed analysis of a simulation time line.
#
# The script takes two required command line arguments: the name of the file to
# plot, and the name of the desired output folder (this folder is erased if
# it already exists).
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules:
#  - numpy for file reading and array operations
#  - argparse for smart command line argument parsing
#  - os for file system checks
#  - shutil for folder manipulation
import numpy as np
import argparse
import os
import shutil

chart_index = 0

class TimelineChart:
  def __init__(self, starts, ends, labels, levels):
    global chart_index
    self._starts = starts
    self._ends = ends
    self._labels = labels
    self._levels = levels
    self._index = chart_index
    chart_index += 1

  def write_callback(self, ofile):
    ofile.write("google.charts.setOnLoadCallback(drawChart{0});".format(
      self._index))

  def write_function(self, ofile):
    ofile.write("function drawChart{0}() {{".format(self._index))
    ofile.write("var data = new google.visualization.DataTable();")
    ofile.write("data.addColumn({type: 'string', id: 'Level'});")
    ofile.write("data.addColumn({type: 'string', id: 'Label'});")
    ofile.write("data.addColumn({type: 'number', id: 'Start'});")
    ofile.write("data.addColumn({type: 'number', id: 'End'});")
    ofile.write("data.addRows([")
    for i in range(len(self._labels)):
      ofile.write("['Level {0}', '{1}', {2}, {3}],".format(
        self._levels[i], self._labels[i], 1000. * self._starts[i],
        1000. * self._ends[i]))
    ofile.write("]);")
    ofile.write("var chart = new google.visualization.Timeline("
                " document.getElementById('timelinechart{0}'));".format(
      self._index))
    ofile.write("chart.draw(data)}")

  def write_div(self, ofile):
    ofile.write("<div id=\"timelinechart{0}\"".format(self._index))
    ofile.write(" style=\"width: 100%; height: 250px;\"></div>")

class PieChart:
  def __init__(self, wedges, labels, level):
    global chart_index
    self._wedges = wedges
    self._labels = labels
    self._index = chart_index
    self._level = level
    chart_index += 1

  def write_callback(self, ofile):
    ofile.write("google.charts.setOnLoadCallback(drawChart{0});".format(
      self._index))

  def write_function(self, ofile):
    ofile.write("function drawChart{0}() {{".format(self._index))
    ofile.write("var data = google.visualization.arrayToDataTable([")
    ofile.write("['Task', 'Total time (s)'],")
    for i in range(len(self._wedges)):
      ofile.write("['{0}', {1}],".format(self._labels[i], self._wedges[i]))
    ofile.write("]);")
    ofile.write("var options = {{title: 'Level {0}'}};".format(self._level))
    ofile.write("var chart = new google.visualization.PieChart("
                " document.getElementById('piechart{0}'));".format(self._index))
    ofile.write("chart.draw(data, options)}")

  def write_div(self, ofile):
    ofile.write("<div id=\"piechart{0}\"".format(self._index))
    ofile.write(" style=\"width: 100%; height: 500px;\"></div>")

class HTMLFile:
  def __init__(self):
    self._charts = []
    self._widths = []

  def add_chart(self, chart, full_width = False):
    self._charts.append(chart)
    if full_width:
      self._widths.append(12)
    else:
      self._widths.append(6)

  def write(self, fname):
    ofile = open(fname, "w")
    ofile.write("""<html><head><script type="text/javascript"
src="https://www.gstatic.com/charts/loader.js"></script>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js"></script>
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js"></script>
<script type="text/javascript">
google.charts.load('current', {'packages':['corechart', 'timeline']});""")

    for chart in self._charts:
      chart.write_callback(ofile)
    for chart in self._charts:
      chart.write_function(ofile)
    ofile.write("</script></head><body>")
    ofile.write("<div class=\"jumbotron text-center\">")
    ofile.write("<h1>Simulation time line</h1></div>")
    ofile.write("<div class=\"container-fluid\"><div class=\"row\">")
    for i in range(len(self._charts)):
      ofile.write("<div class=\"col-sm-{0}\">".format(self._widths[i]))
      self._charts[i].write_div(ofile)
      ofile.write("</div>")
    ofile.write("</div></div>")
    ofile.write("</body></html>")

# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("--file", "-f", action="store", required=True)
argparser.add_argument("--output-folder", "-o", action="store", required=True)
args = argparser.parse_args()

# create the output folder
if os.path.exists(args.output_folder):
  shutil.rmtree(args.output_folder)
os.mkdir(args.output_folder)

# load the time log
data = np.loadtxt(
    args.file,
    delimiter="\t",
    dtype={
        "names": ("id", "pid", "depth", "tic", "toc", "start", "end", "label"),
        "formats": ("u4", "u4", "u4", "u8", "u8", "f8", "f8", ("U", 100)),
    },
)

total_time = data["end"].max()

fields = np.unique(data["label"])
times = {}
for field in fields:
  subdata = data[data["label"] == field]
  times[field] = {}
  times[field]["total_time"] = (subdata["end"] - subdata["start"]).sum()
  times[field]["level"] = subdata["depth"][0]

minlevel = data["depth"].min()
maxlevel = data["depth"].max()
hfile = HTMLFile()
hfile.add_chart(TimelineChart(data["start"], data["end"], data["label"],
                              data["depth"]), True)
for level in range(minlevel, maxlevel + 1):
  wedges = []
  labels = []
  dist = []
  for field in fields:
    if times[field]["level"] == level:
      wedges.append(times[field]["total_time"])
      labels.append(field)
      dist.append(0.1)
  hfile.add_chart(PieChart(wedges, labels, level))

hfile.write("{0}/index.html".format(args.output_folder))
