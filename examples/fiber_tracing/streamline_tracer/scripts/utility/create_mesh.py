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
#  Note, that you need to create the pickle file `rings` first. This can be done with `scripts/utility/create_rings.py`.
#  Instead of calling this python script, use $OPENDIHU_HOME/examples/fiber_tracing/streamline_tracer/scripts/run_evaluation.sh
#
# usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x> [<improve_mesh> [<pickle output filename> <bin output filename>]]]]]]"

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
import datetime
import time

# add to python path in order to load stl_create_mesh
if os.environ.get("OPENDIHU_HOME"):
  sys.path.append(os.path.join(os.environ.get("OPENDIHU_HOME"), "scripts/geometry_manipulation"))
import stl_create_mesh

duration = 0

# constant parameters
triangulation_type = 1          # 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
parametric_space_shape = 3      # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid, 4 = like 3 but move grid points such that mean distance between points in world space gets optimal
max_area_factor = 2.            # only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
show_plot = False
debug = False  
debugging_stl_output = True    # enable this to get the png files which give useful information about the process of generating the mesh
improve_mesh = False           # if the spacing of the nodes should be improved
pickle_output_filename = "mesh"
bin_output_filename = "mesh.bin"

n_points_x = 4    # number of points per dimension for the initial triangulation, assuming a cube reference area (4*n_points_x on circumference for unit circle)

n_grid_points_x = n_points_x+1   # grid width of generated 2d mesh
n_grid_points_y = n_points_x+1

if parametric_space_shape == 0:  # for unit circle 
  n_grid_points_y = 20
  
  
if len(sys.argv) < 2:
  print("usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x (-1=auto)> [<improve_mesh> [<pickle output filename> <bin output filename>]]]]]]")
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
  argument_value = int(sys.argv[4])
  # if set to a negative value, do not use this value
  if argument_value > 0:
    n_grid_points_x = argument_value
    n_grid_points_y = n_grid_points_x
  
if len(sys.argv) >= 6:
  improve_mesh = False if int(sys.argv[5]) == 0 else True
  
if len(sys.argv) >= 8:
  pickle_output_filename = sys.argv[6]
  bin_output_filename = sys.argv[7]
  
print("triangulation_type: {} (0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance)".format(triangulation_type))
print("parametric_space_shape: {} (0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid, 4 = like 3 but move grid points such that mean distance between points in world space gets optimal)".format(parametric_space_shape))
print("n_points_x: {}".format(n_points_x))
print("n_grid_points_x: {}".format(n_grid_points_x))
print("n_grid_points_y: {}".format(n_grid_points_y))
print("improve_mesh: {}".format(improve_mesh))
print("pickle_output_filename: \"{}\"".format(pickle_output_filename))
print("bin_output_filename:    \"{}\"".format(bin_output_filename))
  
# read in loops
with open('rings', 'rb') as f:
  loops = pickle.load(f)
  
n_loops = len(loops)
print("{} loops".format(n_loops))

# sample loop with 4*n_points_x equidistant points
n_points = 4*n_points_x
boundary_point_loops,lengths = stl_create_mesh.rings_to_boundary_points(loops, n_points)

# triangle lists for debugging output to stl files
out_triangulation_world_space = []
markers_boundary_points_world_space = []
out_triangulation_parametric_space = []
grid_triangles_world_space = []
grid_triangles_parametric_space = []
markers_grid_points_parametric_space = []
markers_grid_points_world_space = []
out_3d_mesh_triangles = []
 
loop_grid_points = []  # list of grid point, for every slice, only contains loops that are not empty  
distances_between_world_mesh_nodes_std = []   # list of distances between neighboring grid points
relative_distances_between_world_mesh_nodes_std = []   # list of relative distances between neighbouring grid points, relative per slice

# loop over all loops of boundary points
for loop_no,(boundary_points,length) in enumerate(zip(boundary_point_loops,lengths)):
    
  print("")
  print("Loop {}/{} with {} boundary points, length: {}".format(loop_no, n_loops, len(boundary_points), length))
  
  # create 2D mesh with boundary_points
  show_plot = False
  grid_points_world_space,duration_1d = stl_create_mesh.create_planar_mesh(boundary_points, loop_no, n_points, \
    n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output,\
    [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,\
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
    f.write("# triangulation_type; parametric_space_shape; improve_mesh; debugging_stl_output; n_grid_points_x; n_grid_points_y; number of rings; standard deviation of distance; standard deviation of relative distances (distance/mean distance on every slice); duration\n")
with open("mesh_quality.csv", "a") as f:
  f.write("{};{};{};{};{};{};{};{};{};{}\n".\
  format(triangulation_type, parametric_space_shape, improve_mesh, debugging_stl_output, n_grid_points_x, n_grid_points_y, len(loops),\
    standard_deviation_distance_between_world_mesh_nodes,\
    standard_deviation_relative_distance_between_world_mesh_nodes,duration))

data = stl_create_mesh.create_3d_mesh(loop_grid_points, n_grid_points_x, n_grid_points_y, debugging_stl_output, out_3d_mesh_triangles)

# write result to pickle file
try:
  with open(pickle_output_filename, 'wb') as f:
    pickle.dump(data, f)
    print("Saved pickle file \"{}\"".format(pickle_output_filename))
except:
  print("Could not save pickle file \"{}\".".format(pickle_output_filename))

# write result to .bin file
node_positions = data["node_positions"]
n_elements = data["n_linear_elements_per_coordinate_direction"]
n_fibers_x = n_elements[0]+1
n_fibers_y = n_elements[1]+1
n_points_whole_fiber = n_elements[2]+1

# create output file
with open(bin_output_filename, "wb") as outfile:
  
  header_str = "opendihu binary fibers version 2"
  header_length_raw = 40
  outfile.write(struct.pack('32s', str(header_str).encode('utf-8')))
  outfile.write(struct.pack('i',header_length_raw))
   
  # write parameter[0]: n_fibers_total
  n_fibers = n_fibers_x * n_fibers_y
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
        for ii in range(3):
          double_raw = struct.pack('d', point[ii])
          outfile.write(double_raw)
          
  print("Saved {}x{}={} fibers with {} points each to \"{}\".".format(n_fibers_x, n_fibers_y, n_fibers_x*n_fibers_y, n_points_whole_fiber, bin_output_filename))
  
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
  print("current working directory: {}".format(os.getcwd()))
  write_stl(markers_boundary_points_world_space,   "out/mesh_02_boundary_points_w.stl", "boundary points")
  write_stl(out_triangulation_world_space,       "out/mesh_03_triangulation_w.stl", "triangulation world space")
  write_stl(out_triangulation_parametric_space,  "out/mesh_04_triangulation_p.stl", "triangulation parametric space")
  write_stl(grid_triangles_parametric_space,     "out/mesh_05_grid_triangles_p.stl","grid parametric space")
  write_stl(markers_grid_points_parametric_space,"out/mesh_06_grid_points_p.stl",   "grid points parametric space")
  write_stl(grid_triangles_world_space,          "out/mesh_07_grid_triangles_w.stl","grid world space")
  write_stl(markers_grid_points_world_space,     "out/mesh_08_grid_points_w.stl",   "grid points world space")
  write_stl(out_3d_mesh_triangles,               "out/mesh_09_3d_mesh_w.stl",       "3d mesh world space")
