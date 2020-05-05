#!python3
# -*- coding: utf-8 -*-
#
# Utility to extract the bottom or top points of a muscle fiber file and set them as seed points to do fiber tracing for the tendons.
#
# usage: set_seed_points.py <muscle_fiber_filename> <tracing_mesh_filename> <output filename> <is_bottom>

import sys
import numpy as np
import os
import pickle
import struct
import datetime

muscle_fiber_filename = ""
tracing_mesh_filename = ""
output_filename = ""
is_bottom = True

n_fibers_x_extract = 1

if len(sys.argv) >= 3:
  muscle_fiber_filename = sys.argv[1]
  tracing_mesh_filename = sys.argv[2]
  output_filename = sys.argv[3]
  is_bottom = False if (int)(sys.argv[4])==0 else True
else:
  print("usage: set_seed_points.py <muscle_fiber_filename> <tracing_mesh_filename> <output filename> <is_bottom>")
  quit()

print("input muscle fiber file: \"{}\"".format(muscle_fiber_filename))
print("tracing_mesh_filename: \"{}\"".format(tracing_mesh_filename))
print("is_bottom: {}".format(is_bottom))

# parse fiber file
with open(muscle_fiber_filename, "rb") as infile:
  
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
  
  if "version 2" in str(header_str):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_x = parameters[2]
    n_fibers_y = parameters[3]
  
  print("file contains {} x {} fibers\n".format(n_fibers_x, n_fibers_y))
  print("date: {:%d.%m.%Y %H:%M:%S}".format(datetime.datetime.fromtimestamp(parameters[8])))
  
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
    
# extract bottom or top points
seed_points = []
for j in range(n_fibers_y):
  for i in range(n_fibers_x):
    if is_bottom:
      k = 0
      mesh_index = k*n_fibers_x*n_fibers_y + j*n_fibers_x + i
      point = streamlines[j*n_fibers_x + i][0]
      #point[2] -= 1e-10
    else:
      k = n_points_whole_fiber-1
      mesh_index = k*n_fibers_x*n_fibers_y + j*n_fibers_x + i
      point = streamlines[j*n_fibers_x + i][-1]
      #point[2] += 1e-10
      
    print("  add point ({},{},{}), index {}, p: {}".format(i,j,k,mesh_index,point))
    seed_points.append(point)
  
# set seed points

# read in data from pickle file
with open(tracing_mesh_filename, 'rb') as f:
  data = pickle.load(f, encoding='latin1')

  # assign seed points
  data["seed_points"] = seed_points
  
  # write out data
  print("Write seed points to file: \"{}\"".format(output_filename))
  with open(output_filename, 'wb') as f:
    pickle.dump(data, f)
  
