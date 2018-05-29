#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import collections
import copy
import scipy.spatial
import os
import pymp

import stl
from stl import mesh
from sets import Set
from svg.path import parse_path
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier

# source: https://github.com/kliment/Printrun/blob/master/printrun/stltool.py
def ray_triangle_intersection(ray_origin, ray_direction, (v1, v2, v3)):
  """
  Möller–Trumbore intersection algorithm in pure python
  Based on http://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
  """
  eps = 1e-6
  edge1 = v2 - v1
  edge2 = v3 - v1
  #print "edges: ",edge1,edge2
  h = np.cross(ray_direction, edge2)
  a = edge1.dot(h)
    
  #print "pvec=",pvec,"det=",det
  if abs(a) < eps:
    return False
    
  f = 1. / a
  s = ray_origin - v1
  u = f * s.dot(h)
  #print "u=",u
  
  if u < 0. or u > 1.:
    return False
    
  q = np.cross(s, edge1)
  v = f * ray_direction.dot(q)
  
  if v < 0. or u + v > 1.:
    return False

  t = f * edge2.dot(q)
  #print "t=",t
  if t > eps:
    return True
  else:
    return False

# count the number of intersections of a ray origin + t*direction with the stl_mesh
def get_n_intersections(origin, direction, stl_mesh):
  
  n_intersections_array = pymp.shared.array((1,), dtype='uint8')
  n_intersections_array[0] = 0
  
  # loop over all triangles in mesh
  n_triangles = len(stl_mesh.points)
  with pymp.Parallel(4) as par:
    for i in par.range(0,n_triangles):
      p = stl_mesh.points[i]
      # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

      p1 = np.array(p[0:3])
      p2 = np.array(p[3:6])
      p3 = np.array(p[6:9])
      
      center = (p1+p2+p3)/3.
      if np.linalg.norm(origin-center) < 1e-13:
        continue

      # check if ray intersects triangle
      intersects = ray_triangle_intersection(origin, direction, (p1, p2, p3))

      if intersects:
        with par.lock:
          n_intersections_array[0] += 1
      
  return n_intersections_array[0]

# parse command line arguments
outfile = "out.stl"
infile = "in.stl"
direction1 = np.array([0.0,1.0,0.3])
direction2 = np.array([1.0,0.0,0.3])

# normalize direction
direction1 = direction1/np.linalg.norm(direction1)
direction2 = direction2/np.linalg.norm(direction2)

if len(sys.argv) < 2:
  print "usage: remove_inside_triangles <input file> [<output file>]"
  sys.exit(0)

if len(sys.argv) >= 2:
  if os.path.isfile(sys.argv[1]):
    infile = sys.argv[1]
  else:
    print "File \"{}\" does not exists".format(sys.argv[1])
    sys.exit(0)
  
if len(sys.argv) >= 3:
  outfile = sys.argv[2]
else:
  outfile = os.path.splitext(infile)[0]+"_out.stl"

print "Input file: \"{}\"".format(infile)
print "Output file: \"{}\"".format(outfile)

stl_mesh = mesh.Mesh.from_file(infile)

out_triangles = []
origin_index = 0

n_is_inside_1 = 0
n_is_inside_2 = 0

min = [None,None,None]
max = [None,None,None]

n_triangles = len(stl_mesh.points)

# loop over triangles in mesh
for p in stl_mesh.points:
  # p contains the 9 entries [p1x p1y p1z p2x p2y p2z p3x p3y p3z] of the triangle with corner points (p1,p2,p3)

  origin_index += 1

  p1 = np.array(p[0:3])
  p2 = np.array(p[3:6])
  p3 = np.array(p[6:9])
  
  center = (p1+p2+p3)/3.
  
  # update boundary box
  for i in range(3):
    if min[i] == None or min[i] > center[i]:
      min[i] = center[i]
    if max[i] == None or max[i] < center[i]:
      max[i] = center[i]
  
  is_inside = False
  
  # check how often a ray intersects the surface
  n_intersections = get_n_intersections(center, direction1, stl_mesh)    
  #is_inside_1 = (n_intersections%2 == 1)
  is_inside_1 = (n_intersections > 0)
  
  n_intersections = get_n_intersections(center, -direction1, stl_mesh)    
  is_inside_2 = (n_intersections > 0)

  n_intersections = get_n_intersections(center, direction2, stl_mesh)    
  is_inside_3 = (n_intersections > 0)
  
  n_intersections = get_n_intersections(center, -direction2, stl_mesh)
  is_inside_4 = (n_intersections > 0)

  print "\b\b\b\b\b{:.2f}%".format(float(origin_index) / n_triangles * 100.0)

  if is_inside_1:
    n_is_inside_1 += 1
  if is_inside_2:
    n_is_inside_2 += 1
    
  is_inside = (is_inside_1 and is_inside_2) and (is_inside_3 and is_inside_4)
  
  if not is_inside:
    out_triangles += [[p1, p2, p3]]
    
print "n_is_inside_1:",n_is_inside_1
print "n_is_inside_2:",n_is_inside_2
print "ray directions: ", direction1, direction2

print "bounding box: x[",min[0],",",max[0],"], y[",min[1],",",max[1],"], z[",min[2],",",max[2],"], "

# create output mesh
triangles = out_triangles
n_triangles = len(out_triangles)
print "n_triangles: ",n_triangles

# Create the mesh
out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
for i, f in enumerate(triangles):
  out_mesh.vectors[i] = f
  #for j in range(3):
    #print "set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
    
    
out_mesh.update_normals()

out_mesh.save(outfile, mode=stl.Mode.ASCII)
