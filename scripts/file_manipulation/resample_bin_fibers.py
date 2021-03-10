#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This script reads a fiber file (.bin) and creates a new one with different mesh width.
# It interpolates in the given mesh. The four "corner" fibers will remain the same, all others are interpolated
# For .bin files.
#
# usage: resample_bin_fibers.py <input filename> <output filename> <n_fibers_x> [<n_fibers_y> [<n_points_z>]]

import sys, os
import numpy as np
import datetime
import struct
import time

input_filename = "fibers.bin"
output_filename = "{}.out.bin".format(input_filename)

n_points_z = None
if len(sys.argv) >= 4:
  input_filename = sys.argv[1]
  output_filename = sys.argv[2]
  n_fibers_x = (int)(sys.argv[3])
  n_fibers_y = n_fibers_x
else:
  print("usage: resample_bin_fibers.py <input filename> <output filename> <n_fibers_x> [<n_fibers_y> [<n_points_z>]]")
  quit()

if len(sys.argv) >= 5:
  n_fibers_y = (int)(sys.argv[4])
  
if len(sys.argv) >= 6:
  n_points_z = (int)(sys.argv[5])

print("input filename:  {}".format(input_filename))
print("output filename: {}".format(output_filename))

if n_fibers_x < 2 or n_fibers_y < 2:
  print("Number of fibers is {} x {}, but it has to be at least 2 x 2!".format(n_fibers_x, n_fibers_y))
  quit()

# source https://stackoverflow.com/a/1094933/10290071
def format_filesize(num, suffix='B'):
  for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
    if abs(num) < 1024.0:
      return "%3.1f %s%s" % (num, unit, suffix)
    num /= 1024.0
  return "%.1f %s%s" % (num, 'Yi', suffix)

# parse input file
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
    
  n_fibers_total_input = parameters[0]
  n_points_whole_fiber_input = parameters[1]
  n_fibers_x_input = (int)(np.sqrt(parameters[0]))
  n_fibers_y_input = n_fibers_x_input
  
  # set n_points_z if not specified
  if n_points_z is None:
    n_points_z = n_points_whole_fiber_input
  
  n_fibers_total = n_fibers_x*n_fibers_y
  
  if "version 2" in str(header_str):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_x_input = parameters[2]
    n_fibers_y_input = parameters[3]

  filesize_input = (32+header_length) + (n_points_whole_fiber_input * 3 * 8) * n_fibers_total_input
  filesize_output = (32+header_length) + (n_points_z * 3 * 8) * n_fibers_total
  print("date of input file: {:%d.%m.%Y %H:%M:%S}\n".format(datetime.datetime.fromtimestamp(parameters[8])))
  print("input file has a mesh        {} x {} = {} fibers à {} points, in total {} points, {}".format(
    n_fibers_x_input, n_fibers_y_input, n_fibers_total_input, n_points_whole_fiber_input, n_fibers_total_input*n_points_whole_fiber_input, 
    format_filesize(filesize_input)))
  print("output file will have a mesh {} x {} = {} fibers à {} points, in total {} points, {}".format(
    n_fibers_x, n_fibers_y, n_fibers_total, n_points_z, n_fibers_total*n_points_z,
    format_filesize(filesize_output)))
  print("\nNote that you can also refine meshes using the `refine` program in `parallel_fiber_estimation`, ")
  print("which is less flexible, but might be faster for large meshes as not all data is loaded to RAM first.\n")
  
  fibers = []
  # loop over fibers
  for fiber_no in range(n_fibers_total_input):
    fiber = []
    
    # loop over points of fiber
    for point_no in range(n_points_whole_fiber_input):
      point = []
      
      # parse point
      for i in range(3):
        double_raw = infile.read(8)
        value = struct.unpack('d', double_raw)[0]
        point.append(value)
        
      fiber.append(np.array(point))
    fibers.append(fiber)
    
  # create output file
  with open(output_filename,"wb") as outfile:
    
    infile.seek(0)
    header_str = "opendihu binary fibers version 2".format(header_str)
    outfile.write(struct.pack('32s', str(header_str).encode('utf-8')))
    outfile.write(header_length_raw)
     
    # write parameter[0]: n_fibers_total
    outfile.write(struct.pack('i', n_fibers_total))
    
    # write parameter[1]: n_points_z
    outfile.write(struct.pack('i', n_points_z))
    
    # write parameter[2]: n_fibers_x
    outfile.write(struct.pack('i', n_fibers_x))
    
    # write parameter[3]: n_fibers_y
    outfile.write(struct.pack('i', n_fibers_y))
    
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
    
    # write fibers, iterate over new fibers
    for j in range(n_fibers_y):
      
      # determine coordinates in (x,y) plane in [0,1]^2 
      coordinate_y = j / (n_fibers_y-1)
      
      # determine coordinates of surrounding rectangle of fibers, 
      # the current fiber (x,y)=(i,j) is contained in [i_input_previous,i_input_next] x [j_input_previous,j_input_next] of the input fiber mesh
      j_input_previous = (int)(coordinate_y * (n_fibers_y_input-1))
      j_input_next = j_input_previous + 1
      if j_input_next >= n_fibers_y_input:
        j_input_next = n_fibers_y_input-1
      
      factor_y = coordinate_y * (n_fibers_y_input-1) - j_input_previous
      
      for i in range(n_fibers_x):
        
        # determine coordinates in (x,y) plane in [0,1]^2 
        coordinate_x = i / (n_fibers_x-1)
        
        # determine coordinates of surrounding rectangle of fibers, 
        # the current fiber (x,y)=(i,j) is contained in [i_input_previous,i_input_next] x [j_input_previous,j_input_next] of the input fiber mesh
        i_input_previous = (int)(coordinate_x * (n_fibers_x_input-1))
        i_input_next = i_input_previous + 1
        if i_input_next >= n_fibers_x_input:
          i_input_next = n_fibers_x_input-1
        
        factor_x = coordinate_x * (n_fibers_x_input-1) - i_input_previous
        
        # 2 3
        # 0 1
        fiber0 = fibers[j_input_previous*n_fibers_x_input + i_input_previous]
        fiber1 = fibers[j_input_previous*n_fibers_x_input + i_input_next]
        fiber2 = fibers[j_input_next    *n_fibers_x_input + i_input_previous]
        fiber3 = fibers[j_input_next    *n_fibers_x_input + i_input_next]
        
        # loop over points of fiber
        for k in range(n_points_z):
          
          coordinate_z = k / (n_points_z-1)
          k_input_previous = (int)(coordinate_z * (n_points_whole_fiber_input-1))
          k_input_next = k_input_previous + 1
          if k_input_next >= n_points_whole_fiber_input:
            k_input_next = n_points_whole_fiber_input-1
            
          factor_z = coordinate_z * (n_points_whole_fiber_input-1) - k_input_previous
          
          # interpolate new point
          point = \
              (1-factor_z)*(1-factor_y)*(1-factor_x)*fiber0[k_input_previous] \
            + (1-factor_z)*(1-factor_y)*factor_x    *fiber1[k_input_previous] \
            + (1-factor_z)*factor_y    *(1-factor_x)*fiber2[k_input_previous] \
            + (1-factor_z)*factor_y    *factor_x    *fiber3[k_input_previous] \
            + factor_z    *(1-factor_y)*(1-factor_x)*fiber0[k_input_next] \
            + factor_z    *(1-factor_y)*factor_x    *fiber1[k_input_next] \
            + factor_z    *factor_y    *(1-factor_x)*fiber2[k_input_next] \
            + factor_z    *factor_y    *factor_x    *fiber3[k_input_next]
          
          # write point to file
          for component_no in range(3):
            double_raw = struct.pack('d', point[component_no])
            outfile.write(double_raw)
            
    print("File {} written.".format(output_filename))
    
