#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script reads a *.bin file and outputs a *.vtr file to be viewed with Paraview.
#
# usage: ./bin_to_vtk.py <fibers.bin input file> <output_file (without .vtr)>

import sys, os
import numpy as np
import struct
import stl
from stl import mesh
import datetime
import pickle
from pyevtk.hl import gridToVTK     # pip3 install pyevtk

input_filename = "../../examples/electrophysiology/input/13x13fibers.bin"
if len(sys.argv) >= 2:
  input_filename = sys.argv[1]

output_filename = "{}".format(input_filename)
if len(sys.argv) >= 3:
  output_filename = sys.argv[2]
  
print("input filename: {}".format(input_filename))
print("output filename: {}".format(output_filename))

# parse input file and save all streamlines
with open(input_filename, "rb") as infile:
  
  # parse header
  bytes_raw = infile.read(32)
  header_str = struct.unpack('32s', bytes_raw)[0]
  
  header_length_raw = infile.read(4)
  header_length = struct.unpack('i', header_length_raw)[0]
  #header_length = 32+8
  parameters = []
  for i in range(int(header_length/4) - 1):
    int_raw = infile.read(4)
    value = struct.unpack('i', int_raw)[0]
    parameters.append(value)
    
  n_fibers_total = parameters[0]
  n_points_whole_fiber = parameters[1]
  n_fibers_x = (int)(np.sqrt(parameters[0]))
  n_fibers_y = n_fibers_x
  
  if "version 2" in header_str:   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_x = parameters[2]
    n_fibers_y = parameters[3]
  
  print("header: {}".format(header_str))
  print("nFibersTotal:      {n_fibers} = {n_fibers_x} x {n_fibers_y}".format(n_fibers=parameters[0], n_fibers_x=n_fibers_x, n_fibers_y=n_fibers_y))
  print("nPointsWholeFiber: {}".format(parameters[1]))
  if "version 2" not in header_str:
    print("nBoundaryPointsXNew: {}".format(parameters[2]))
    print("nBoundaryPointsZNew: {}".format(parameters[3]))
  print("nFineGridFibers_:  {}".format(parameters[4]))
  print("nRanks:            {}".format(parameters[5]))
  print("nRanksZ:           {}".format(parameters[6]))
  print("nFibersPerRank:    {}".format(parameters[7]))
  print("date:              {:%d.%m.%Y %H:%M:%S}".format(datetime.datetime.fromtimestamp(parameters[8])))
  
  streamlines = []
  n_streamlines_valid = 0
  n_streamlines_invalid = 0
  
  # loop over fibers
  for streamline_no in range(n_fibers_total):
    streamline = []
    streamline_valid = True
    
    # loop over points of fiber
    for point_no in range(n_points_whole_fiber):
      point = []
      
      # parse point
      for i in range(3):
        double_raw = infile.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
        
      # check if point is valid
      if point[0] == 0.0 and point[1] == 0.0 and point[2] == 0.0:
        if streamline_valid:
          coordinate_x = streamline_no % n_fibers_x
          coordinate_y = (int)(streamline_no / n_fibers_x)
          print("Error: streamline {}, ({},{})/({},{}) is invalid ({}. point)".format(streamline_no, coordinate_x, coordinate_y, n_fibers_x, n_fibers_y, point_no))
          print("streamline so far: ",streamline[0:10])
        streamline_valid = False
      streamline.append(point)
      
    if streamline_valid:
      n_streamlines_valid += 1
    else:
      n_streamlines_invalid += 1
      streamline = []
    streamlines.append(streamline)
  
  print("n valid: {}, n invalid: {}".format(n_streamlines_valid, n_streamlines_invalid))
  
# create geometry for vtk structured grid
n_points_x = n_fibers_x
n_points_y = n_fibers_y
n_points_z = n_points_whole_fiber
n_points = n_points_x * n_points_y * n_points_z

print("create geometry with {} x {} x {} = {} points".format(n_points_x, n_points_y, n_points_z, n_points))

positions_x = np.zeros((n_points_x, n_points_y, n_points_z))
positions_y = np.zeros((n_points_x, n_points_y, n_points_z))
positions_z = np.zeros((n_points_x, n_points_y, n_points_z))

field_fiber_no = np.zeros((n_points_x, n_points_y, n_points_z))

# loop over points of geometry and set positions
for k in range(n_points_z):
  for j in range(n_points_y):
    for i in range(n_points_x):
      
      point = streamlines[j*n_fibers_x + i][k]
      positions_x[i,j,k] = point[0]
      positions_y[i,j,k] = point[1]
      positions_z[i,j,k] = point[2]
      
      fiber_no = j*n_points_x + i
      field_fiber_no[i,j,k] = fiber_no

# write vtk file
gridToVTK(output_filename, positions_x, positions_y, positions_z, pointData = {"fiber_no" : field_fiber_no})
