#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script outputs the bounding box of the file. 
# usage: get_bounding_box_fibers.py <filename>

import sys, os
import numpy as np
import struct
import stl
from stl import mesh
import datetime
import pickle
import time

input_filename = "fibers.bin"

if len(sys.argv) == 2:
  input_filename = sys.argv[1]
  
else:
  print("usage: get_bounding_box_fibers.py <filename>\n The output will be: xmin xmax ymin ymax zmin zmax")
  sys.exit(0)

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
  
  bounding_box_min = [None,None,None]
  bounding_box_max = [None,None,None]
  
  # loop over fibers
  for streamline_no in range(n_fibers_total):
    
    # loop over points of fiber
    for point_no in range(n_points_whole_fiber):
      point = []
      
      # parse point
      for i in range(3):
        double_raw = infile.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
        
      # update bounding box
      for i in range(3):
        if bounding_box_min[i] is None:
          bounding_box_min[i] = point[i]
          bounding_box_max[i] = point[i]
        bounding_box_min[i] = min(bounding_box_min[i],point[i])
        bounding_box_max[i] = max(bounding_box_max[i],point[i])
      
  print("{} {} {} {} {} {}".format(bounding_box_min[0],bounding_box_max[0],bounding_box_min[1],bounding_box_max[1],
    bounding_box_min[2],bounding_box_max[2]))
    
