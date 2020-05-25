#!python3
# -*- coding: utf-8 -*-
#
# Utility to change the top or bottom nodes of one mesh to be a copy of another mesh.
# For *.bin files.
#
# usage: set_nodes_to_match_other_mesh.py <input_filename mesh where to change nodes> <input_filename mesh to which to adapt> <output filename> <is_bottom if change bottom layer of mesh else top layer>
#
#
#


import sys
import numpy as np
import os
import pickle
import struct
import datetime
import time

input_filename_1 = ""
input_filename_2 = ""
output_filename = ""
is_bottom = True

if len(sys.argv) >= 4:
  input_filename_1 = sys.argv[1]
  input_filename_2 = sys.argv[2]
  output_filename = sys.argv[3]
  is_bottom = False if (int)(sys.argv[4])==0 else True
else:
  print("usage: set_nodes_to_match_other_mesh.py <input_filename mesh where to change nodes> <input_filename mesh to which to adapt> <output filename> <is_bottom if change bottom layer of mesh else top layer>")
  quit()

print("input mesh to change: \"{}\"".format(input_filename_1))
print("input mesh to match:  \"{}\"".format(input_filename_2))
print("output_filename:      \"{}\"".format(output_filename))
print("is_bottom: {}".format(is_bottom))

# parse fiber file 1
with open(input_filename_1, "rb") as infile:
  
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
    
  n_fibers_1_total = parameters[0]
  n_points_1_whole_fiber = parameters[1]
  n_fibers_1_x = (int)(np.sqrt(parameters[0]))
  n_fibers_1_y = n_fibers_1_x
  
  if "version 2" in str(header_str):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_1_x = parameters[2]
    n_fibers_1_y = parameters[3]
  
  print("file \"{}\" contains {} x {} fibers with {} points each\n".format(input_filename_1, n_fibers_1_x, n_fibers_1_y, n_points_1_whole_fiber))
  print("created: {:%d.%m.%Y %H:%M:%S}".format(datetime.datetime.fromtimestamp(parameters[8])))
  
  fibers_mesh_1 = []
  
  # loop over fibers
  for fiber_no in range(n_fibers_1_total):
    fiber = []
    
    # loop over points of fiber
    for point_no in range(n_points_1_whole_fiber):
      point = []
      
      # parse point
      for i in range(3):
        double_raw = infile.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
        
      fiber.append(point)
      
    fibers_mesh_1.append(fiber)
    
# parse fiber file 2
with open(input_filename_2, "rb") as infile:
  
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
    
  n_fibers_2_total = parameters[0]
  n_points_2_whole_fiber = parameters[1]
  n_fibers_2_x = (int)(np.sqrt(parameters[0]))
  n_fibers_2_y = n_fibers_2_x
  
  if "version 2" in str(header_str):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_2_x = parameters[2]
    n_fibers_2_y = parameters[3]
  
  print("file \"{}\" contains {} x {} fibers with {} points each\n".format(input_filename_2, n_fibers_2_x, n_fibers_2_y, n_points_2_whole_fiber))
  print("created: {:%d.%m.%Y %H:%M:%S}".format(datetime.datetime.fromtimestamp(parameters[8])))
  
  fibers_mesh_2 = []
  
  # loop over fibers
  for fiber_no in range(n_fibers_2_total):
    fiber = []
    
    # loop over points of fiber
    for point_no in range(n_points_2_whole_fiber):
      point = []
      
      # parse point
      for i in range(3):
        double_raw = infile.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
        
      fiber.append(point)
      
    fibers_mesh_2.append(fiber)
      
  if n_fibers_1_x != n_fibers_2_x or n_fibers_1_y != n_fibers_2_y:
    print("Number of fibers / mesh dimension mismatches!")

# write out result
# write bin file
with open(output_filename, "wb") as outfile:
  
  header_str = "opendihu binary fibers version 2"
  header_length_raw = 40
  outfile.write(struct.pack('32s', str(header_str).encode('utf-8')))
  outfile.write(struct.pack('i',header_length_raw))
   
  # write parameter[0]: n_fibers_total
  outfile.seek(32+4)
  outfile.write(struct.pack('i', n_fibers_1_total))
  
  # write parameter[1]: n_points_whole_fiber
  outfile.seek(32+2*4)
  outfile.write(struct.pack('i', n_points_1_whole_fiber))
  
  # write parameter[2]: n_fibers_x
  outfile.seek(32+3*4)
  outfile.write(struct.pack('i', n_fibers_1_x))
  
  # write parameter[3]: n_fibers_y
  outfile.seek(32+4*4)
  outfile.write(struct.pack('i', n_fibers_1_y))
  
  # write timestamp
  outfile.seek(32+9*4)
  outfile.write(struct.pack('i', (int)(time.time())))
  
  # write fibers
  for j in range(n_fibers_1_y):
    for i in range(n_fibers_1_x):
      for k in range(n_points_1_whole_fiber):
        point = fibers_mesh_1[j*n_fibers_1_x + i][k]
        
        if is_bottom and k == 0:
          point = fibers_mesh_2[j*n_fibers_1_x + i][n_points_2_whole_fiber-1]
          
        elif not is_bottom and k == n_points_1_whole_fiber-1:
          point = fibers_mesh_2[j*n_fibers_1_x + i][0]
          
        # store point
        for component_no in range(3):
          double_raw = struct.pack('d', point[component_no])
          outfile.write(double_raw)
          
  print("Saved {}x{}={} fibers with {} points each to \"{}\".".format(n_fibers_1_x, n_fibers_1_y, n_fibers_1_x*n_fibers_1_y, n_points_1_whole_fiber, output_filename))
  
