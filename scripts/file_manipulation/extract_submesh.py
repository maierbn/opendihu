#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script reads a fiber file (.bin) and extracts a cuboid of the mesh.
# For .bin files.
#
# usage: extract_submesh.py <input filename> <output filename> <x_begin> <x_end> <y_begin> <y_end> <z_begin> <z_end>

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

if len(sys.argv) >= 9:
  input_filename = sys.argv[1]
  output_filename = sys.argv[2]
  x_begin = int(sys.argv[3])
  x_end = int(sys.argv[4])
  y_begin = int(sys.argv[5])
  y_end = int(sys.argv[6])
  z_begin = int(sys.argv[7])
  z_end = int(sys.argv[8])
else:
  print("usage: extract_submesh.py <input filename> <output filename> <x_begin> <x_end> <y_begin> <y_end> <z_begin> <z_end>")
  quit()

offset = -1  # -1 means compute
if len(sys.argv) >= 5:
  offset = (int)(sys.argv[4])
  
print("input filename: {}\noutput filename: {}".format(input_filename, output_filename))
print("x range: [{},{}]".format(x_begin, x_end))
print("y range: [{},{}]".format(y_begin, y_end))
print("z range: [{},{}]".format(z_begin, z_end))

n_points_extracted_x = x_end - x_begin
n_points_extracted_y = y_end - y_begin
n_points_extracted_z = z_end - z_begin

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
  
  if "version 2" in str(header_str):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_x = parameters[2]
    n_fibers_y = parameters[3]
  
  print("input file has a mesh {} x {} x {}".format(n_fibers_x,n_fibers_y,n_points_whole_fiber))
  print("output file will have a mesh {} x {} x {}\n".format(n_points_extracted_x, n_points_extracted_y, n_points_extracted_z))

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
  
  fibers = []
  # loop over fibers
  for fiber_no in range(n_fibers_total):
    fiber = []
    
    # loop over points of fiber
    for point_no in range(n_points_whole_fiber):
      point = []
      
      # parse point
      for i in range(3):
        double_raw = infile.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
        
      fiber.append(point)
    fibers.append(fiber)
  
  if x_begin >= n_fibers_x or x_end > n_fibers_x \
    or y_begin >= n_fibers_y or y_end > n_fibers_y \
    or z_begin >= n_points_whole_fiber or z_end > n_points_whole_fiber:
    print("Error, selected range does not match mesh in input file.")
    quit()
  
  # create output file
  with open(output_filename,"wb") as outfile:
    
    infile.seek(0)
    header_str = "opendihu binary fibers version 2".format(header_str)
    outfile.write(struct.pack('32s', str(header_str).encode('utf-8')))
    outfile.write(header_length_raw)
     
    # write parameter[0]: n_fibers_total
    n_fibers = n_points_extracted_x * n_points_extracted_y
    outfile.write(struct.pack('i', n_fibers))
    
    # write parameter[1]: n_points_extracted_z
    outfile.write(struct.pack('i', n_points_extracted_z))
    
    # write parameter[2]: n_points_extracted_x
    outfile.write(struct.pack('i', n_points_extracted_x))
    
    # write parameter[3]: n_points_extracted_y
    outfile.write(struct.pack('i', n_points_extracted_y))
    
    # write parameter[4]:
    outfile.write(struct.pack('i', 0))
    
    # write parameter[5]: n_ranks
    outfile.write(struct.pack('i', 1))
    
    # write parameter[6]: n_ranks_z
    outfile.write(struct.pack('i', 1))
    
    # write parameter[7]: n_fibers_per_rank
    outfile.write(struct.pack('i', 1))
    
    # write parameters[8]: timestamp
    outfile.write(struct.pack('i', (int)(time.time())))
    
    # write fiber
    for y in range(y_begin, y_end):
      for x in range(x_begin, x_end):
        fiber = fibers[y*n_fibers_x + x]
        
        # loop over points of fiber
        for point_no in range(z_begin, z_end):
          point = fiber[point_no]
          
          # parse point
          for component_no in range(3):
            double_raw = struct.pack('d', point[component_no])
            outfile.write(double_raw)
            
    print("File {} written.".format(output_filename))
    
