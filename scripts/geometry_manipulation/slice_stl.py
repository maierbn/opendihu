#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Extracts a subset of a given volume bounded by an STL mesh. The extracted part is between two planes that are parallel to the z direction.
# Input is a tubular surface as STL file, output is an enclosed surface of the extracted part (also the bottom and top planes are triangulated, such that the output is fully enclosed.)
#
# usage:  ./slice_stl.py [<input filename> [<output filename> [<bottom clip> [<top clip> [<n loops>]]]]]

import sys, os
import numpy as np
import stl
from stl import mesh
import pickle

import stl_create_rings
import stl_create_mesh

input_filename = "biceps_full.stl"
output_filename = "sliced.stl"
bottom_clip = 20
top_clip = 100
n_loops = 5

if len(sys.argv) < 2:
  print("usage:  ./slice_stl.py [<input filename> [<output filename> [<bottom clip> [<top clip> [<n loops>]]]]]")
  sys.exit(0)

if len(sys.argv) >= 2:
  input_filename = sys.argv[1]

if len(sys.argv) >= 3:
  output_filename = sys.argv[2]

if len(sys.argv) >= 4:
  bottom_clip = (float)(sys.argv[3])

if len(sys.argv) >= 5:
  top_clip = (float)(sys.argv[4])

if len(sys.argv) >= 6:
  n_loops = max(2,(int)(sys.argv[5]))

print("input_filename: {}".format(input_filename))
print("output_filename: {}".format(output_filename))
print("bottom_clip: {}".format(bottom_clip))
print("top_clip: {}".format(top_clip))

write_output_mesh = False
loops = stl_create_rings.create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh)
n_points = 20

boundary_point_loops,lengths = stl_create_mesh.rings_to_boundary_points(loops, n_points)

triangles = []

n_loops = len(boundary_point_loops)

# mantle triangles
for i in range(n_loops-1):
  for j in range(n_points):
    p0 = boundary_point_loops[i][j]
    p1 = boundary_point_loops[i][(j+1)%n_points]
    p2 = boundary_point_loops[i+1][j]
    p3 = boundary_point_loops[i+1][(j+1)%n_points]
    triangles.append([p0,p1,p3])
    triangles.append([p0,p3,p2])

# bottom and top triangles
is_bottom = True
for z_value,points in zip([bottom_clip, top_clip],[boundary_point_loops[0], boundary_point_loops[-1]]):

  print("z : {}".format(z_value))

  n_points = len(points)

  # project points on xy=z_value plane
  projected_points = []
  for point in points:
    projected_points.append(np.array([point[0], point[1]]))

  projected_points = np.reshape(projected_points, (-1,2))

  # delaunay triangulation of triangle package, adds new points, is constrained (works for concave domains)
  
  import triangle   # sudo pip triangle

  # create delaunay triangulation of points
  segments = np.reshape([[i,i+1] for i in range(n_points)], (n_points,2))
  segments[n_points-1] = np.array([n_points-1,0])
  
  data = {"vertices": projected_points, "segments": segments}

  triangulation = triangle.triangulate(data, 'pq')
  triangulated_projected_points = np.array(triangulation['vertices'])
  
  # transform projected points back to 3D points
  points = []
  for projected_point in triangulated_projected_points:
    points.append(np.array([projected_point[0], projected_point[1], z_value]))
  
  # update n_points
  n_points = len(points)
  points = np.reshape(points, (-1,3))
    
  for triangle in triangulation["triangles"]:
    if is_bottom:
      triangles.append([points[triangle[0]], points[triangle[2]], points[triangle[1]]])
    else:
      triangles.append([points[triangle[0]], points[triangle[1]], points[triangle[2]]])
  is_bottom = False

#---------------------------------------
# Create the mesh
out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
for i, f in enumerate(triangles):
  out_mesh.vectors[i] = f
#out_mesh.update_normals()

out_mesh.save(output_filename)
print("saved {} triangles to \"{}\"".format(len(triangles),output_filename))
