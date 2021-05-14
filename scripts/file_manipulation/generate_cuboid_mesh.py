#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script creates a cuboid mesh file in .bin format.
#
# usage: generate_cuboid_mesh.py <output filename> <n points x> <n points y> <n points z> <mesh width x> <mesh width y> <mesh width z> [<offset x> <offset y> <offset z>]

import sys, os
import numpy as np
import struct
import datetime
import pickle
import time

output_filename = "cuboid.bin"

offset_x = 0
offset_y = 0
offset_z = 0

if len(sys.argv) >= 8:
  output_filename = sys.argv[1]
  n_points_x = int(sys.argv[2])
  n_points_y = int(sys.argv[3])
  n_points_z = int(sys.argv[4])
  mesh_width_x = float(sys.argv[5])
  mesh_width_y = float(sys.argv[6])
  mesh_width_z = float(sys.argv[7])
else:
  print("usage: generate_cuboid_mesh.py <output filename> <n points x> <n points y> <n points z> <mesh width x> <mesh width y> <mesh width z> [<offset x> <offset y> <offset z>]")
  quit()

if len(sys.argv) >= 11:
  offset_x = float(sys.argv[8])
  offset_y = float(sys.argv[9])
  offset_z = float(sys.argv[10])

print("output filename: \"{}\"".format(output_filename))
print("size:       {} x {} x {}".format(n_points_x, n_points_y, n_points_z))
print("mesh width: {} x {} x {}".format(mesh_width_x, mesh_width_y, mesh_width_z))
print("offset:     {} x {} x {}".format(offset_x, offset_y, offset_z))
print("bounding box will be: [{},{}] x [{},{}] x [{},{}]".format(offset_x,offset_x + (n_points_x-1)*mesh_width_x, offset_y,offset_y + (n_points_y-1)*mesh_width_y, offset_z,offset_z + (n_points_z-1)*mesh_width_z))

# create output file
with open(output_filename,"wb") as outfile:
  
  header_str = "opendihu binary fibers version 2"
  outfile.write(struct.pack('32s', str(header_str).encode('utf-8')))
  outfile.write(struct.pack('i', 40))  # header length
   
  # write parameter[0]: n_fibers_total
  n_fibers = n_points_x * n_points_y
  outfile.write(struct.pack('i', n_fibers))
  
  # write parameter[1]: n_points_z
  outfile.write(struct.pack('i', n_points_z))
  
  # write parameter[2]: n_points_x
  outfile.write(struct.pack('i', n_points_x))
  
  # write parameter[3]: n_points_y
  outfile.write(struct.pack('i', n_points_y))
  
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
  for y in range(n_points_y):
    for x in range(n_points_x):
      
      # loop over points of fiber
      for z in range(n_points_z):
        point = [offset_x + x*mesh_width_x, offset_y + y*mesh_width_y, offset_z + z*mesh_width_z]
        
        # parse point
        for component_no in range(3):
          double_raw = struct.pack('d', point[component_no])
          outfile.write(double_raw)
          
  print("File {} written.".format(output_filename))
  
