#!/usr/bin/env ../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This scripts reads a fibers.bin file and evaluates the quality of the mesh, how the spacing between the fibers is.
#
# usage: ./mesh_evaluate_quality.py [<filename>]

import sys, os
import numpy as np
import struct
import stl
from stl import mesh
import datetime
import pickle

input_filename = "fibers.bin"
include_boundary = True

if len(sys.argv) >= 2:
  input_filename = sys.argv[1]

print("\nEvaluate the quality of the generated mesh by computing variances of relative element lengths and angles.")
print("Input filename: {}".format(input_filename))

with open(input_filename, "rb") as infile:
  
  # parse header
  bytes_raw = infile.read(32)
  scenario_name = struct.unpack('32s', bytes_raw)[0].decode('utf-8')
  print(scenario_name)
  
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
  
  if False:
    print("nFibersTotal:      {n_fibers} = {n_fibers_x} x {n_fibers_x}".format(n_fibers=parameters[0], n_fibers_x=n_fibers_x))
    print("nPointsWholeFiber: {}".format(parameters[1]))
    print("nBoundaryPointsXNew: {}".format(parameters[2]))
    print("nBoundaryPointsZNew: {}".format(parameters[3]))
    print("nFineGridFibers_:  {}".format(parameters[4]))
    print("nRanks:            {}".format(parameters[5]))
    print("nRanksZ:           {}".format(parameters[6]))
    print("nFibersPerRank:    {}".format(parameters[7]))
    print("date:              {:%d.%m.%Y %H:%M:%S}".format(datetime.datetime.fromtimestamp(parameters[8])))
    
    input("Press any key to continue.")
  
  # parse data from file
  n_points_x = n_fibers_x
  n_points_y = n_points_x
  n_points_z = parameters[1]
  points = np.zeros((n_points_z, n_points_y, n_points_x, 3))
  
  n_fibers_valid = 0
  n_fibers_invalid = 0
  
  # loop over fibers
  for fiber_no in range(n_fibers_total):
    fiber_value = True
    
    # loop over points of fiber
    for z in range(n_points_whole_fiber):
      point = np.zeros((3))
      
      # parse point
      for i in range(3):
        double_raw = infile.read(8)
        value = struct.unpack('d', double_raw)[0]
        point[i] = value
        
      # check if point is valid
      if point[0] == 0.0 and point[1] == 0.0 and point[2] == 0.0:
        if fiber_value:
          coordinate_x = fiber_no % n_fibers_x
          coordinate_y = (int)(fiber_no / n_fibers_x)
          print("Error: streamline {}, ({},{})/({},{}) is invalid ({}. point)".format(fiber_no, coordinate_x, coordinate_y, n_fibers_x, n_fibers_y, z))
          #print("streamline so far: ",streamline[0:10])
        fiber_value = False
      
      x = fiber_no % n_fibers_x
      y = (int)(fiber_no / n_fibers_x)
      points[z,y,x,:] = point
      
    if fiber_value:
      n_fibers_valid += 1
    else:
      n_fibers_invalid += 1
  
  print("Fibers n valid: {}, n invalid: {}, {} points for fiber. Computing all angles for the mesh can take a bit ...".format(n_fibers_valid, n_fibers_invalid, n_points_whole_fiber))
  
  # compute variance
  # loop over z levels
  variance_sum = 0
  if include_boundary:
    x_start = 0
    x_end = n_points_x
    y_start = 0
    y_end = n_points_y
  else:
    x_start = 1
    x_end = n_points_x-1
    
    y_start = 1
    y_end = n_points_y-1
  
  angles = []
  for z in range(n_points_whole_fiber):
    edge_lengths = []
    for y in range(y_start,y_end):
      for x in range(x_start,x_end):
        if x < x_end-1:
          edge_lengths.append(np.linalg.norm(points[z,y,x,:] - points[z,y,x+1,:]))
        elif y < y_end-1:
          edge_lengths.append(np.linalg.norm(points[z,y,x,:] - points[z,y+1,x,:]))

        #     u1
        #  u2 p  u0
        #     u3
        # angle between u0 and u1
        if x < x_end-1 and y < y_end-1:
          u0 = -points[z,y,x,:] + points[z,y,x+1,:]
          if np.linalg.norm(u0) != 0: 
            u0 /= np.linalg.norm(u0)
          u1 = -points[z,y,x,:] + points[z,y+1,x,:]
          if np.linalg.norm(u1) != 0: 
            u1 /= np.linalg.norm(u1)
          angle = np.arccos(u0.dot(u1))
          angles.append(angle)
          
        # angle between u1 and u2
        if x > 0 and y < y_end-1:
          u1 = -points[z,y,x,:] + points[z,y+1,x,:]
          if np.linalg.norm(u1) != 0: 
            u1 /= np.linalg.norm(u1)
          u2 = -points[z,y,x,:] + points[z,y,x-1,:]
          if np.linalg.norm(u2) != 0: 
            u2 /= np.linalg.norm(u2)
          angle = np.arccos(u1.dot(u2))
          angles.append(angle)
          
        # angle between u2 and u3
        if x > 0 and y > 0:
          u2 = -points[z,y,x,:] + points[z,y,x-1,:]
          if np.linalg.norm(u2) != 0: 
            u2 /= np.linalg.norm(u2)
          u3 = -points[z,y,x,:] + points[z,y-1,x,:]
          if np.linalg.norm(u3) != 0: 
            u3 /= np.linalg.norm(u3)
          angle = np.arccos(u2.dot(u3))
          angles.append(angle)
          
        # angle between u3 and u0
        if x < x_end-1 and y > 0:
          u3 = -points[z,y,x,:] + points[z,y-1,x,:]
          if np.linalg.norm(u3) != 0: 
            u3 /= np.linalg.norm(u3)
          u0 = -points[z,y,x,:] + points[z,y,x+1,:]
          if np.linalg.norm(u0) != 0: 
            u0 /= np.linalg.norm(u0)
          angle = np.arccos(u0.dot(u3))
          angles.append(angle)

    mean = np.mean(edge_lengths)
    variance = np.var(edge_lengths)
    #print("z: {}, mean: {}, var: {}, var/mean^2: {}".format(z, mean, variance, variance / (mean*mean)))
    variance_sum += variance / (mean*mean)

  variance = variance_sum / n_points_whole_fiber
  variance_angles = np.var(angles)
  mean_angle = np.mean(angles)
  
  print("lengths variance: {}, angles variance: {}, angles mean: {}".format(variance, variance_angles, mean_angle/np.pi*180))
  
  # write to log file
  with open("variances.csv","a") as f:
    f.write("{};{};{};{};{}\n".format(scenario_name,input_filename,n_fibers_valid,variance,variance_angles))
  print("appended to variances.csv")
