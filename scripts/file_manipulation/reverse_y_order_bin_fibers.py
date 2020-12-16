#!/usr/bin/env ../../../dependencies/python/install/bin/python3
# -*- coding: utf-8 -*-
#
# This script reverses the order of the numbering in y direction.
# For .bin files.
#
# usage: reverse_y_order_bin_fibers.py <input filename> [<output filename>]

import sys, os
import numpy as np
import struct
import stl
from stl import mesh
import datetime
import pickle
import time

input_filename = "fibers.bin"
output_filename = "{}.out.bin".format(input_filename)

if len(sys.argv) >= 2:
  input_filename = sys.argv[1]

output_filename = "{}.reversed".format(input_filename)
if len(sys.argv) >= 3:
  output_filename = sys.argv[2]
  
print("input file:  {}".format(input_filename))
print("output file: {}".format(output_filename))

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
  
  if "version 2" in header_str.decode('utf-8'):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_x = parameters[2]
    n_fibers_y = parameters[3]
  
  print("reverse y numbering in file with {} x {} fibers\n".format(n_fibers_x,n_fibers_y))

  print("header: {}".format(header_str))
  print("nFibersTotal:      {n_fibers} = {n_fibers_x} x {n_fibers_y}".format(n_fibers=parameters[0], n_fibers_x=n_fibers_x, n_fibers_y=n_fibers_y))
  print("nPointsWholeFiber: {}".format(parameters[1]))
  if "version 2" not in str(header_str):
    print("nBoundaryPointsXNew: {}".format(parameters[2]))
    print("nBoundaryPointsZNew: {}".format(parameters[3]))
  print("nFineGridFibers_:  {}".format(parameters[4]))
  print("nRanks:            {}".format(parameters[5]))
  print("nRanksZ:           {}".format(parameters[6]))
  print("nFibersPerRank:    {}".format(parameters[7]))
  print("date:              {:%d.%m.%Y %H:%M:%S}".format(datetime.datetime.fromtimestamp(parameters[8])))
  
  streamlines = []
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
      streamline.append(point)
    streamlines.append(streamline)
    
  # create output file
  with open(output_filename,"wb") as outfile:
    
    infile.seek(0)
    header_str = str.encode("opendihu binary fibers version 2".format(header_str))
    outfile.write(struct.pack('32s', header_str))
    outfile.write(header_length_raw)
     
    # write parameter[0]: n_fibers_total
    n_fibers = n_fibers_x * n_fibers_y
    outfile.seek(32+4)
    outfile.write(struct.pack('i', n_fibers))
    
    # write parameter[1]: n_points_whole_fiber
    outfile.seek(32+2*4)
    outfile.write(struct.pack('i', n_points_whole_fiber))
    
    # write parameter[2]: n_fibers_x
    outfile.seek(32+3*4)
    outfile.write(struct.pack('i', n_fibers_x))
    
    # write parameter[3]: n_fibers_y
    outfile.seek(32+4*4)
    outfile.write(struct.pack('i', n_fibers_y))
    
    # write timestamp
    outfile.seek(32+9*4)
    outfile.write(struct.pack('i', (int)(time.time())))
    
    # write fiber
    for y in reversed(range(n_fibers_y)):
      for x in range(n_fibers_x):
        fiber = streamlines[y*n_fibers_x + x]
        # loop over points of fiber
        for point_no in range(n_points_whole_fiber):
          point = fiber[point_no]
          
          # parse point
          for component_no in range(3):
            double_raw = struct.pack('d', point[component_no])
            outfile.write(double_raw)
            
    print("File {} written.".format(output_filename))
    
