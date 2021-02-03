#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script scales the points in the fiber file by given scaling factors [scale_x,scale_y,scale_z]. 
# The origin of the scaling is at [x,y,z] or (0,0,0) if not specified.
# For .bin files.
#
# usage: scale_bin_fibers.py <input filename> <output filename> <scale_x> <scale_y> <scale_z> [<x>, <y>, <z>]

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

origin_x = 0
origin_y = 0
origin_z = 0

if len(sys.argv) >= 6:
  input_filename = sys.argv[1]
  output_filename = sys.argv[2]
  scale_x = float(sys.argv[3])
  scale_y = float(sys.argv[4])
  scale_z = float(sys.argv[5])
  
  if len(sys.argv) >= 9:
    origin_x = float(sys.argv[6])
    origin_y = float(sys.argv[7])
    origin_z = float(sys.argv[8])
else:
  print("usage: scale_bin_fibers.py <input filename> <output filename> <scale_x> <scale_y> <scale_z> [<x>, <y>, <z>]")
  sys.exit(0)
  
print("input file:         {}".format(input_filename))
print("output file:        {}".format(output_filename))
print("scaling factors:    [{}, {}, {}]".format(scale_x, scale_y, scale_z))
print("origin point that will not change: [{}, {}, {}]".format(origin_x, origin_y, origin_z))

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
  
  print("header: {}".format(header_str))
  print("nFibersTotal:      {n_fibers} = {n_fibers_x} x {n_fibers_y}".format(n_fibers=parameters[0], n_fibers_x=n_fibers_x, n_fibers_y=n_fibers_y))
  print("nPointsWholeFiber: {}".format(parameters[1]))
  if "version 2" not in header_str.decode('utf-8'):
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
  
  bounding_box_min = [None,None,None]
  bounding_box_max = [None,None,None]
  
  # create output file
  with open(output_filename,"wb") as outfile:
    
    infile.seek(0)
    header_str = str.encode("opendihu binary fibers version 2".format(header_str))
    outfile.write(struct.pack('32s', header_str))
    outfile.write(header_length_raw)
     
    # write parameter[0]: n_fibers_total
    outfile.seek(32+4)
    outfile.write(struct.pack('i', n_fibers_total))
    
    # write parameter[1]: n_fibers_total
    outfile.seek(32+2*4)
    outfile.write(struct.pack('i', parameters[1]))
    
    # write parameter[2]: n_fibers_x
    outfile.seek(32+3*4)
    outfile.write(struct.pack('i', n_fibers_x))
    
    # write parameter[3]: n_fibers_y
    outfile.seek(32+4*4)
    outfile.write(struct.pack('i', n_fibers_y))
    
    for i in range(5,9):
      outfile.seek(32+i*4)
      outfile.write(struct.pack('i', parameters[i-1]))
    
    # write timestamp
    outfile.seek(32+9*4)
    outfile.write(struct.pack('i', (int)(time.time())))
    
    # write fiber
  
    for streamline_no in range(n_fibers_total):
      fiber = streamlines[streamline_no]
      
      # loop over points of fiber
      for point_no in range(n_points_whole_fiber):
        point = fiber[point_no]
        
        # transform point
        point[0] = (point[0] - origin_x) * scale_x + origin_x
        point[1] = (point[1] - origin_y) * scale_y + origin_y
        point[2] = (point[2] - origin_z) * scale_z + origin_z
        
        # update bounding box
        for i in range(3):
          if bounding_box_min[i] is None:
            bounding_box_min[i] = point[i]
            bounding_box_max[i] = point[i]
          bounding_box_min[i] = min(bounding_box_min[i],point[i])
          bounding_box_max[i] = max(bounding_box_max[i],point[i])
        
        # write point
        for component_no in range(3):
          double_raw = struct.pack('d', point[component_no])
          outfile.write(double_raw)
          
    print("New bounding box: [{},{}] x [{},{}] x [{},{}],\n size: {} x {} x {}".
      format(bounding_box_min[0],bounding_box_max[0],bounding_box_min[1],bounding_box_max[1],bounding_box_min[2],bounding_box_max[2],
             bounding_box_max[0]-bounding_box_min[0],bounding_box_max[1]-bounding_box_min[1],bounding_box_max[2]-bounding_box_min[2]))
    print("File {} written.".format(output_filename))
    
