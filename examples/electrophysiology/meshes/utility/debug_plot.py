#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to visualize pickle file with seed points in it
#

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle

input_filename = "processed_meshes/left_biceps_brachii_05_tendon1.pickle"

if len(sys.argv) >= 2:
  input_filename = sys.argv[1]

print("input filename: {}".format(input_filename))

# read in data from pickle file
with open(input_filename, 'rb') as f:
  data = pickle.load(f, encoding='latin1')

# load node positions and seed points
node_positions = data["node_positions"]
seed_points = data["seed_points"]
nx,ny,nz = data["n_linear_elements_per_coordinate_direction"]
nx += 1
ny += 1
nz += 1

print("node_positions:{}".format(node_positions))
print("seed_ponits: {}".format(seed_points))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

points = node_positions

x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]

# points
ax.scatter(x_list, y_list, z_list, c='b',s=10)

# lines
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

for j in range(ny):
  for i in range(nx):
    for k in range(nz-1):
      p0 = node_positions[k*nx*ny + j*nx + i]
      p1 = node_positions[(k+1)*nx*ny + j*nx + i]
      ax.plot([p0[0],p1[0]], [p0[1],p1[1]], [p0[2],p1[2]], '-', color=colors[i%len(colors)])

points = seed_points

x_list = [x for [x,y,z] in points]
y_list = [y for [x,y,z] in points]
z_list = [z for [x,y,z] in points]

# points
ax.scatter(x_list, y_list, z_list, c='r',s=20)



plt.show()
