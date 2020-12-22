#!/usr/bin/env ../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This script calls the stl_create_mesh.create_planar_mesh function from geometry_manipulation/stl_create_mesh.py

import sys,os
import pickle
import stl_create_mesh

import datetime
now = datetime.datetime.now()
print(" ======= debug_create_planar_mesh.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

filename = "dump_mesh_2d_28813.py"
if len(sys.argv) > 1:
  filename = sys.argv[1]

print("use file \"{}\"".format(filename))

# read in data
with open(filename, 'rb') as f:
  data = pickle.load(f)
  print("data:")
  print(data)
  print("len: ",len(data))

boundary_points_faces = data
stl_create_mesh.create_3d_mesh_from_boundary_points_faces(boundary_points_faces)

sys.exit(0)











#-------------------- create_planar_mesh
debugging_stl_output = True
if debugging_stl_output:
  out_triangulation_world_space = []
  markers_boundary_points_world_space = []
  out_triangulation_parametric_space = []
  grid_triangles_world_space = []
  grid_triangles_parametric_space = []
  markers_grid_points_parametric_space = []
  markers_grid_points_world_space = []
  debugging_output_lists = [out_triangulation_world_space, markers_boundary_points_world_space, out_triangulation_parametric_space, grid_triangles_world_space, grid_triangles_parametric_space,markers_grid_points_parametric_space, markers_grid_points_world_space]
else:
  debugging_output_lists = []
    
# create_planar_mesh(boundary_points, loop_no, n_points, n_grid_points_x, n_grid_points_y, triangulation_type, parametric_space_shape, max_area_factor, show_plot, debugging_stl_output, debugging_output_lists)
triangulation_type = data[5]
parametric_space_shape = data[6]
max_area_factor = data[7]
show_plot = data[8]

# constant parameters
triangulation_type = 1  # 0 = scipy, 1 = triangle, 2 = center pie (2 is best), 3 = minimized distance
parametric_space_shape = 3   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid
max_area_factor = 80.    # only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
show_plot = False
improve_mesh = True

grid_points_world_space,duration_1d = stl_create_mesh.create_planar_mesh(data[0], data[1], data[2], data[3], data[4], triangulation_type, parametric_space_shape, max_area_factor, improve_mesh, show_plot, debugging_stl_output, debugging_output_lists)
  
print("duration: {}".format(duration_1d))
