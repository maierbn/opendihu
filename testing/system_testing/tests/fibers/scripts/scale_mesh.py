#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This script scales the node positions of a mesh and writes it to another file.
# This is necessary to adjust the mesh in same way the fibers were scaled.
#
# usage: ./scale_mesh.py <input filename> <output filename> <scaling_factor>

import datetime
now = datetime.datetime.now()
print(" ======= scale_mesh.py =======")
print(now.strftime("%d/%m/%Y %H:%M:%S"))

import sys, os
import numpy as np
import pickle
import timeit

scaling_factor = 0.0627878
  
if len(sys.argv) < 4:
  print("usage: ./scale_mesh.py <input filename> <output filename> <scaling_factor>")
  sys.exit(0)

input_filename = sys.argv[1]
output_filename = sys.argv[2]
scaling_factor = (float)(sys.argv[3])
  
print("input_filename: {}".format(input_filename))
print("output_filename: {}".format(output_filename))
print("scaling_factor: {}".format(scaling_factor))

with open(input_filename, 'rb') as f:
  data = pickle.load(f, encoding='latin1')
#
#  "node_positions": node_positions, 
#  "linear_elements": linear_elements, 
#  "quadratic_elements": quadratic_elements, 
#  "seed_points": seed_points,
#  "bottom_nodes": bottom_node_indices,
#  "top_nodes": top_node_indices,
#  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
#  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,

# scale node_positions
for i,node_position in enumerate(data["node_positions"]):
  data["node_positions"][i] = list(np.array(data["node_positions"][i]) * scaling_factor)
  
# write data to output file 

with open(output_filename, 'wb') as f:
  pickle.dump(data, f)
