#!python3
# -*- coding: utf-8 -*-
#
# This script reads a pickle file with streamlines and outputs a .bin file.
#
# usage: ./convert_pickle_to_bin.py <input pickle file> <output file>

import datetime
now = datetime.datetime.now()
print(" ======= convert_pickle_to_bin.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

import sys, os
import numpy as np
import pickle
import time

input_filename = "streamlines.csv"
output_filename = "out"
  
if len(sys.argv) < 3:
  print("usage: usage: ./convert_pickle_to_bin.py <input pickle file> <output file>")
  sys.exit(0)
else:
  if os.path.isfile(sys.argv[1]):
    input_filename = sys.argv[1]
  else:
    print("File \"{}\" does not exist.".format(sys.argv[1]))
    sys.exit(0)
    
  output_filename = sys.argv[2]
    
print("input file:  \"{}\"".format(input_filename))
print("output file: \"{}\"".format(output_filename))

# read in data from pickle file
with open(input_filename, 'rb') as f:
  streamlines = pickle.load(f, encoding='latin1')

if len(streamlines) == 0:
  print("Error, file contains no streamlines.")
  quit()

n_streamlines = len(streamlines)
n_fibers_x = (int)(np.sqrt(n_streamlines))
n_fibers_y = (int)(n_streamlines / n_fibers_x)
n_fibers = n_fibers_x * n_fibers_y

if n_fibers_x * n_fibers_y != n_streamlines:
  print("Error number of fibers {} is not equal to {} x {} = {} ".format(n_streamlinse, n_fibers_x, n_fibers_y, n_fibers))
  
n_points_whole_fiber = len(streamlines[0])

# write bin file
with open(output_filename, "wb") as outfile:
  
  header_str = "opendihu binary fibers version 2"
  header_length_raw = 40
  outfile.write(struct.pack('32s', str(header_str).encode('utf-8')))
  outfile.write(struct.pack('i',header_length_raw))
   
  # write parameter[0]: n_fibers_total
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
  
  # write fibers
  for j in range(n_fibers_y):
    for i in range(n_fibers_x):
      for k in range(n_points_whole_fiber):
        point = node_positions[k*n_fibers_x*n_fibers_y + j*n_fibers_x + i]
        
        # store point
        for component_no in range(3):
          double_raw = struct.pack('d', point[component_no])
          outfile.write(double_raw)
          
  print("Saved {}x{}={} fibers with {} points each to \"{}\".".format(n_fibers_x, n_fibers_y, n_fibers_x*n_fibers_y, n_points_whole_fiber, output_filename))
  
