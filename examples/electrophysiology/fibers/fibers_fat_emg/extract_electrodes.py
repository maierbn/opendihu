#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to extract a number of electrodes from the given data.
# arguments: <n_extracted_points_xy> <n_extracted_points_z> <path>
#

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

input_filename = "build_release/out/biceps/electrodes.csv"
output_filename = "build_release/out/biceps/electrodes_out.csv"
n_extracted_points_xy = 2
n_extracted_points_z = 2
offset_xy = 0
offset_z = 0


# parse file
if len(sys.argv) == 1:
  print("usage: {} <n_extracted_points_xy> <n_extracted_points_z> <offset_xy> <offset_z> <input file> <output file>".format(sys.argv[0]))
  quit()
if len(sys.argv) > 1:
  n_extracted_points_xy = (int)(sys.argv[1])
if len(sys.argv) > 2:
  n_extracted_points_z = (int)(sys.argv[2])
if len(sys.argv) > 3:
  offset_xy = (int)(sys.argv[3])
if len(sys.argv) > 4:
  offset_z = (int)(sys.argv[4])
if len(sys.argv) > 5:
  input_filename = sys.argv[5]
if len(sys.argv) > 6:
  output_filename = sys.argv[6]
  
print("input file: {}".format(input_filename))
print("output file: {}".format(output_filename))
print("n electrodes to extract: {} x {}".format(n_extracted_points_xy, n_extracted_points_z))
print("offset: {} x {}".format(offset_xy, offset_z))

# parse csv file
data = np.genfromtxt(input_filename, delimiter=';')
t_list = data[:,1]
n_points = (int)(data[0][2])

# parse timestamps which are in the first column and are not parsed by genfromtxt
timestamps = []
with open(input_filename) as file:
  for line in file:
    timestamp = line[0:line.find(";")]
    timestamps.append(timestamp)
timestamps = timestamps[1:]
  

# determine number of electrodes in x direction
# extract all z positions (along muscle)
position_data = data[0][3:3+3*n_points]
z_positions = np.array([position_data[3*i+2] for i in range(n_points)])

# compute differences in z position between neighbour points, this will have a jump (high increment) where a new line starts
differences = z_positions[1:] - z_positions[0:-1]

# determine the indices where the increment is higher than 3 times the mean 
jumps = np.where(differences > np.mean(differences)*3)[0]

# number of points in x/y-direction (across muscle) is position of first jump + 1
n_points_xy = jumps[0] + 1
n_points_z = (int)(n_points / n_points_xy)

n_extracted_points = n_extracted_points_xy*n_extracted_points_z

print("Input file contains {} x {} = {} points".format(n_points_xy, n_points_z, n_points))
print("Extract {} x {} = {} points, between ({},{}) and ({},{})".format(n_extracted_points_xy,n_extracted_points_z,n_extracted_points,offset_xy,offset_z,offset_xy+n_extracted_points_xy,offset_z+n_extracted_points_z))

if offset_xy + n_extracted_points_xy > n_points_xy:
  print("Error, too many points in xy direction ({}, maximum: {})".format(offset_xy+n_extracted_points_xy, n_points_xy))
  quit()
  
if offset_z + n_extracted_points_z > n_points_z:
  print("Error, too many points in z direction ({}, maximum: {})".format(offset_z+n_extracted_points_z, n_points_z))
  quit()

# write out file
with open(output_filename, "w") as output_file:
  
  output_file.write("#timestamp;t;n_points")
  
  # write header
  for i in range(n_extracted_points):
    output_file.write(";p{}_x;p{}_y;p{}_z".format(i,i,i))

  for i in range(n_extracted_points):
    output_file.write(";p{}_value".format(i))
  output_file.write("\n")

  # write data
  for timestep_no in range(len(t_list)):
    
    # write timestamp
    output_file.write(timestamps[timestep_no])
    
    # write t
    output_file.write(";")
    output_file.write(str(data[timestep_no,1]))
    
    # write n_points
    output_file.write(";")
    output_file.write(str(n_extracted_points))
    
    # write point positions
    # loop over extracted points
    for j in range(n_extracted_points_z):
      for i in range(n_extracted_points_xy):
        point_index = (j + offset_z) * n_points_xy + i + offset_xy
        
        for k in range(3):
          output_file.write(";")
          output_file.write(str(data[timestep_no,3 + 3*point_index + k]))
  
    # write value
    for j in range(n_extracted_points_z):
      for i in range(n_extracted_points_xy):
        point_index = (j + offset_z) * n_points_xy + i + offset_xy
        
        output_file.write(";")
        output_file.write(str(data[timestep_no,3 + 3*n_points + point_index]))

    output_file.write("\n")

print("Created file \"{}\".".format(output_filename))

