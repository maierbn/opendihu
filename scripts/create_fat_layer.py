#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This script creates a fat and skin layer on top of a muscle mesh.
#
# usage: ./create_fat_layer.py [<fibers input file> [<output file of fat tissue> [<thickness in cm> [<y_size>]]]]

import sys, os
import numpy as np
import struct
import stl
from stl import mesh
import datetime
import pickle
import time

input_filename = "../examples/electrophysiology/input/13x13fibers.bin"

if len(sys.argv) == 1:
  print("usage: create_fat_layer.py [<fibers input file> [<output file of fat tissue> [<thickness in cm> [<y_size>]]]]\n default values: thickness = 1cm, y_size = 5\n")

if len(sys.argv) >= 2:
  input_filename = sys.argv[1]

output_filename = "{}_fat.bin".format(input_filename)
y_size = 5
thickness = 1.0  #[cm]
  
if len(sys.argv) >= 3:
  output_filename = sys.argv[2]

if len(sys.argv) >= 4:
  thickness = (float)(sys.argv[3])

if len(sys.argv) >= 5:
  y_size = (int)(sys.argv[4])

print("input filename:  {}".format(input_filename))
print("output filename: {}".format(output_filename))
print("thickness:       {} cm".format(thickness))
print("y_size:          {}".format(y_size))

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
  
  if "version 2" in str(header_str):   # the version 2 has number of fibers explicitly stored and thus also allows non-square dimension of fibers
    n_fibers_x = parameters[2]
    n_fibers_y = parameters[3]
  
  print("header: {}".format(str(header_str)))
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
  
  print("file offset: {}".format(infile.tell()))
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
      if point[0] == 0.0 and point[1] == 0.0 and point[2] == 0.0 and False:
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
  
  # create mesh, oriented in x and z direction along fibers, in y direction normal to muscle belly
  # set the sampling stride to 1, the mesh will be sampled by the script that does the partitioning in the simulation
  x_stride = 1
  z_stride = 1
  
  result_n_points_x = len(range(0, n_fibers_x+n_fibers_y-1, x_stride))
  result_n_points_y = y_size
  
  # for z range, ensure that last value is contained
  z_range = list(range(0, n_points_whole_fiber, z_stride))+([n_points_whole_fiber-1] if (n_points_whole_fiber-1)%z_stride !=0 else [])
  result_n_points_z = len(z_range)
  result_n_points = result_n_points_x * result_n_points_y * result_n_points_z
  
  print("create mesh of size {}x{}x{} = {}".format(result_n_points_x, result_n_points_y, result_n_points_z, result_n_points))
  
  result_mesh = []
  for k in z_range:
    for j in range(0, y_size):
      for i in range(0, n_fibers_x+n_fibers_y-1, x_stride):
        
        #print("i: {}, j: {}, k: {}".format(i,j,k))
        
        # compute cog at z-plane
        cog = np.zeros(3)
        for j2 in range(n_fibers_y):
          for i2 in range(n_fibers_x):
            cog += np.array(streamlines[j2*n_fibers_x + i2][k])
        
        cog /= (n_fibers_x * n_fibers_y)
        
        # compute outward normal
        if i < n_fibers_x:
          p0 = np.array(streamlines[(n_fibers_y-1)*n_fibers_x + i][k])
        else:
          p0 = np.array(streamlines[(n_fibers_y-2-(i-n_fibers_x))*n_fibers_x + (n_fibers_x-1)][k])
        n = -cog + p0
        n /= np.linalg.norm(n, ord=1)
        
        # compute current point
        point = p0 + j*thickness/y_size*n
        
        result_mesh.append(point)
  
  print("done")
  print("n points: {} = {}".format(result_n_points, len(result_mesh)))
  
  # write output file
  with open(output_filename,"wb") as outfile:
    
    infile.seek(0)
    header_str = "opendihu binary fibers version 2".format(header_str)
    outfile.write(struct.pack('32s', str(header_str).encode('utf-8')))
    outfile.write(header_length_raw)
    
    # write header
    outfile.seek(32+4)
    
    # write parameter[0]: n_fibers_total
    n_points = result_n_points_x * result_n_points_y
    outfile.write(struct.pack('i', n_points))
    
    # write parameter[1]: n_points_whole_fiber
    outfile.write(struct.pack('i', result_n_points_z))
    
    # write parameter[2]: result_n_points_x
    outfile.write(struct.pack('i', result_n_points_x))
    
    # write parameter[3]: result_n_points_y
    outfile.write(struct.pack('i', result_n_points_y))
    
    # write parameters[8]: timestamp
    outfile.seek(32+9*4)
    outfile.write(struct.pack('i', (int)(time.time())))
    
    outfile.seek(32+10*4)
    
    # write result mesh in the same file format as fibers
    for j in range(result_n_points_y):
      for i in range(result_n_points_x):
        
        # loop over points of "fiber", i.e. line of points in z-direction
        for k in range(result_n_points_z):
          point = result_mesh[k*result_n_points_x*result_n_points_y + j*result_n_points_x + i]
          
          # parse point
          for component_no in range(3):
            double_raw = struct.pack('d', point[component_no])
            outfile.write(double_raw)
            
    print("File {} written.".format(output_filename))
    
