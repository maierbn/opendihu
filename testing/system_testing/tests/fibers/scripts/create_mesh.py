#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This script creates a hexaeder mesh out of the rings given by create_rings.py. 
# It reads in a pickle file `rings` and creates a pickle file `mesh` containing 
#  {
#  "node_positions": nodes, 
#  "linear_elements": linear_elements, 
#  "quadratic_elements": quadratic_elements, 
#  "seed_points": seed_points,
#  "bottom_nodes": bottom_node_indices,
#  "top_nodes": top_node_indices,
#  "n_linear_elements_per_coordinate_direction": n_linear_elements_per_coordinate_direction,
#  "n_quadratic_elements_per_coordinate_direction": n_quadratic_elements_per_coordinate_direction,
#  }.
#
#  The core functionality is implemented in opendihu/scripts/geometry_manipulation/stl_create_mesh.py. The main functions are 
#  stl_create_mesh.standardize_loop, stl_create_mesh.create_planar_mesh, stl_create_mesh.create_3d_mesh
#
# usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x>]]]]"

import datetime
now = datetime.datetime.now()
print(" ======= create_mesh.py =======")
print(now.strftime("%d/%m/%Y %H:%M:%S"))

import sys, os
import numpy as np
import matplotlib

havedisplay = False
if not havedisplay:
  print("use Agg backend")
  matplotlib.use('Agg')
else:
  print("use Tk backend")

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import collections, patches
import struct
import stl
from stl import mesh
import pickle
import timeit
import stl_create_mesh

duration = 0

# constant parameters
triangulation_type = 1  # 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
parametric_space_shape = 3   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
max_area_factor = 2.    # only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
show_plot = False
debug = False  
debugging_stl_output = True

n_points_x = 4    # number of points per dimension for the initial triangulation, assuming a cube reference area (4*n_points_x on circumference for unit circle)

n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
n_grid_points_y = n_points_x+1

if parametric_space_shape == 0:  # for unit circle 
  n_grid_points_y = 20
  
  
if len(sys.argv) < 2:
  print("usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x>]]]]")
  sys.exit(0)

if len(sys.argv) >= 2:
  triangulation_type = int(sys.argv[1])
  
if len(sys.argv) >= 3:
  parametric_space_shape = int(sys.argv[2])
  
if len(sys.argv) >= 4:
  n_points_x = int(sys.argv[3])    
  n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
  n_grid_points_y = n_points_x+1
  
if len(sys.argv) >= 5:
  n_grid_points_x = int(sys.argv[4])
  n_grid_points_y = n_grid_points_x
  
print("triangulation_type: {}".format(triangulation_type))
print("parametric_space_shape: {}".format(parametric_space_shape))
print("n_points_x: {}".format(n_points_x))
print("n_grid_points_x: {}".format(n_grid_points_x))
print("n_grid_points_y: {}".format(n_grid_points_y))
  
# read in loops
with open('rings', 'rb') as f:
  loops = pickle.load(f)
  
n_loops = len(loops)
print("{} loops".format(n_loops))

# sample loop with 4*n_points_x equidistant points
n_points = 4*n_points_x
border_point_loops,lengths = stl_create_mesh.rings_to_border_points(loops, n_points)

# triangle lists for debugging output to stl files
out_triangulation_world_space = []
markers_border_points_world_space = []
out_triangulation_parametric_space = []
grid_triangles_world_space = []
grid_triangles_parametric_space = []
markers_grid_points_parametric_space = []
markers_grid_points_world_space = []
out_3d_mesh_triangles = []
 
loop_grid_points = []  # list of grid point, for every slice, only contains loops that are not empty  
distances_between_world_mesh_nodes_std = []   # list of distances between neighboring grid points
relative_distances_between_world_mesh_nodes_std = []   # list of relative distances between neighbouring grid points, relative per slice

# loop over all loops of border points
for loop_no,(border_points,length) in enumerate(zip(border_point_loops,lengths)):
    
  print("")
  print("Loop {}/{} with {} border points, length: {}".format(loop_no, n_loops, len(border_points), length))
  
  # create 2D mesh with border_points
  show_plot = False
  grid_points_world_space,duration_1d = stl_create_mesh.create_planar_mesh(border_points, loop_no, n_points, \
    n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, show_plot, debugging_stl_output,\
    [out_triangulation_world_space, markers_border_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
      markers_grid_points_parametric_space, markers_grid_points_world_space])

  duration += duration_1d

  # store grid points in world space of current loop
  loop_grid_points.append(grid_points_world_space)

  distances_current_loop,relative_distances_current_loop = stl_create_mesh.compute_mean_distances(grid_points_world_space,n_grid_points_x,n_grid_points_y)
  
  #print("loop {} has relative_distances_current_loop std {} from values {}".format(loop_no, np.std(relative_distances_current_loop), relative_distances_current_loop)
  
  relative_distances_between_world_mesh_nodes_std.append(np.std(relative_distances_current_loop))
  distances_between_world_mesh_nodes_std.append(np.std(distances_current_loop))
  
standard_deviation_distance_between_world_mesh_nodes = np.mean(distances_between_world_mesh_nodes_std)
standard_deviation_relative_distance_between_world_mesh_nodes = np.mean(relative_distances_between_world_mesh_nodes_std)

# save mean distance and duration
if not os.path.isfile("mesh_quality.csv"):
  with open("mesh_quality.csv", "w") as f:
    f.write("# triangulation_type; parametric_space_shape; n_grid_points_x; n_grid_points_y; number of rings; standard deviation of distance; standard deviation of relative distances (distance/mean distance on every slice); duration\n")
with open("mesh_quality.csv", "a") as f:
  f.write("{};{};{};{};{};{};{};{}\n".\
  format(triangulation_type, parametric_space_shape, n_grid_points_x, n_grid_points_y, len(loops),\
    standard_deviation_distance_between_world_mesh_nodes,\
    standard_deviation_relative_distance_between_world_mesh_nodes,duration))

data = stl_create_mesh.create_3d_mesh(loop_grid_points, n_grid_points_x, n_grid_points_y, debugging_stl_output, out_3d_mesh_triangles)

with open('mesh', 'wb') as f:
  pickle.dump(data, f)
  
# write debugging output stl meshes
def write_stl(triangles, outfile, description):
  # create output mesh
  n_triangles = len(triangles)

  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
    #for j in range(3):
      #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
      
  out_mesh.update_normals()

  out_mesh.save(outfile) #, mode=stl.Mode.ASCI
  print("saved {} triangles to \"{}\" ({})".format(n_triangles,outfile,description))

if debugging_stl_output:
  write_stl(markers_border_points_world_space,   "out/mesh_02_border_points_w.stl", "border points")
  write_stl(out_triangulation_world_space,       "out/mesh_03_triangulation_w.stl", "triangulation world space")
  write_stl(out_triangulation_parametric_space,  "out/mesh_04_triangulation_p.stl", "triangulation parametric space")
  write_stl(grid_triangles_parametric_space,     "out/mesh_05_grid_triangles_p.stl","grid parametric space")
  write_stl(markers_grid_points_parametric_space,"out/mesh_06_grid_points_p.stl",   "grid points parametric space")
  write_stl(grid_triangles_world_space,          "out/mesh_07_grid_triangles_w.stl","grid world space")
  write_stl(markers_grid_points_world_space,     "out/mesh_08_grid_points_w.stl",   "grid points world space")
  write_stl(out_3d_mesh_triangles,               "out/mesh_09_3d_mesh_w.stl",       "3d mesh world space")
