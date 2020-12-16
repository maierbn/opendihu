#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Extrude a svg file to create a 3D volume.

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import collections
import copy
import scipy.spatial
import os

from stl import mesh
from sets import Set
from svg.path import parse_path
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier

def get_n_intersections(center, n_points, svg_path):

  n_intersections = 0
  avg_dist_sum = 0.0
  avg_dist_count = 0
  last_p = None
  for t in np.linspace(0,1,n_points):
    point_complex = svg_path.point(t)
    p = np.array([point_complex.real, point_complex.imag])
    
    if last_p is not None:
      
      current_dist = np.linalg.norm(last_p - p)
      if avg_dist_count != 0:
        avg_dist = avg_dist_sum / avg_dist_count
      
        if current_dist / avg_dist > 3.:
          last_p = p
          continue
      
      avg_dist_sum += current_dist
      avg_dist_count += 1
    
      if p[1] - last_p[1] < 1e-10:
        if last_p[0] < center[0] and p[0] < center[0] and (last_p[1] <= center[1] <= p[1] or p[1] <= center[1] <= last_p[1]):
          n_intersections += 1
      else:
        a = (center[1] - last_p[1]) / (p[1] - last_p[1])
        x = last_p[0] + (p[0] - last_p[0])*a
        
        if p[1] != last_p[1] and (0 <= a <= 1) and x <= center[0]:
          n_intersections += 1
      
    last_p = p
  return n_intersections

def svg_path_to_mesh(svg_path, n_points, output_no, thickness, invert_normals=False):
  
  svg_path = parse_path(path)

  #print "path ",svg_path.d()
  
  points = np.zeros((n_points, 2))
  for (i,t) in enumerate(np.linspace(0,1.0,n_points)):
    point_complex = svg_path.point(t)
    x = point_complex.real
    y = point_complex.imag
    
    points[i,0] = x
    points[i,1] = y
    
  # create triangulation
  if True:
    triangulation = scipy.spatial.Delaunay(points)
    

    # remove outside triangles
    interior_triangles = []
    interior_triangles_indices = []
    for [i1, i2, i3] in triangulation.simplices:
      p1 = np.array(points[i1])
      p2 = np.array(points[i2])
      p3 = np.array(points[i3])
      
      center = (p1+p2+p3)/3.
      
      # check how often horizontal line to center intersects path
      n_intersections = get_n_intersections(center, 200, svg_path)
      
      if n_intersections%2 == 1:
        interior_triangles.append([p1, p2, p3])
        if invert_normals:
          interior_triangles_indices.append([i1, i2, i3])
        else:
          interior_triangles_indices.append([i1, i3, i2])

  else:

    # do triangulation with triangle package which does not create outside triangles
    import triangle   # sudo easy_install triangle

    # create delaunay triangulation of points
    segments = np.reshape([[i,i+1] for i in range(n_points)], (n_points,2))
    segments[n_points-1] = np.array([n_points-1,0])
    
    data = {"vertices": points, "segments": segments}

    triangulation = triangle.triangulate(data, 'pq')
    
    #print triangulation

    print "triangulation vertices: ",len(triangulation["vertices"])
    print "points: ",len(points)

    points = np.array(triangulation["vertices"])
    n_points = len(triangulation["vertices"])

    # process triangulation to extract triangles
    interior_triangles = []
    interior_triangles_indices = []
    for [i1, i2, i3] in triangulation["triangles"]:
      p1 = np.array(points[i1])
      p2 = np.array(points[i2])
      p3 = np.array(points[i3])
      
      interior_triangles.append([p1, p2, p3])
      if invert_normals:
        interior_triangles_indices.append([i1, i2, i3])
      else:
        interior_triangles_indices.append([i1, i3, i2])

  # create extruded object
  mesh_points = np.zeros((2*n_points, 3))
  
  # bottom points
  for i in range(n_points):
    mesh_points[i,:] = list(points[i])+[0.0]
  
  # top points
  for i in range(n_points):
    mesh_points[i+n_points,:] = list(points[i])+[thickness]
  
  # faces
  mesh_faces = []
  
  # bottom faces
  mesh_faces = copy.copy(interior_triangles_indices)
  
  # top faces
  mesh_faces += [[i1+n_points, i3+n_points, i2+n_points] for [i1, i2, i3] in interior_triangles_indices]
  
  # side faces
  avg_dist_sum = 0.0
  avg_dist_count = 0
  for i in range(1, n_points):

    p = np.array([points[i,0], points[i,1]])
    last_p = p
    if i > 0:
      last_p = np.array([points[i-1,0], points[i-1,1]])
      current_dist = np.linalg.norm(last_p - p)
      if avg_dist_count != 0:
        avg_dist = avg_dist_sum / avg_dist_count
      
        if current_dist / avg_dist > 3.:
          last_p = p
          continue
      
      avg_dist_sum += current_dist
      avg_dist_count += 1
    
    v = (last_p - p)*3.0
    c = 0.5*(p + last_p) + np.array([v[1], -v[0]])
    
    # check how often horizontal line to center intersects path
    n_intersections = get_n_intersections(c, 200, svg_path)
        
    if (invert_normals) != (n_intersections%2 == 0):
      mesh_faces += [[i, i-1, (i-1)+n_points]]
      mesh_faces += [[i, (i-1)+n_points, i+n_points]]
    else:
      mesh_faces += [[i-1, i, (i-1)+n_points]]
      mesh_faces += [[(i-1)+n_points, i, i+n_points]]
  
  #print "faces :",len(mesh_faces)
  #print "points: ",len(mesh_points), n_points
  
  # create stl mesh
  mesh_object = mesh.Mesh(np.zeros(len(mesh_faces), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(mesh_faces):
    for j in range(3):
      mesh_object.vectors[i][j] = mesh_points[f[j],:]
        
  if output_no is not None:
    mesh_object.save("mesh_{}.stl".format(output_no))
        
    plt.triplot(points[:,0], points[:,1], interior_triangles_indices)
    plt.plot(points[:,0], points[:,1], 'o')
    plt.savefig("pic_{}.png".format(output_no))
  
  return mesh_object
  
outfile = "out100.stl"
infile = "text.svg"
thickness = 1.0
resolution = 1000
invert_normals = False

if len(sys.argv) < 2:
  print "usage: svg2stl <input file> [<output file> [<thickness (default 1)> [<resolution (default 1000)> [<invert normals (default 0)>]]]]"
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
  outfile = os.path.splitext(infile)[0]+".stl"

if len(sys.argv) >= 4:
  thickness = float(sys.argv[3])
  
if len(sys.argv) >= 5:
  resolution = int(sys.argv[4])
  
if len(sys.argv) >= 6:
  invert_normals = int(sys.argv[5])

print "Input file: \"{}\"".format(infile)
print "Output file: \"{}\"".format(outfile)
print "Thickness: {}".format(thickness)
print "Resolution: {}".format(resolution)
print "invert normals: {}".format(invert_normals)

# extract paths starting with [d="]
paths = []
with open(infile) as file:
  contents = file.read()
  
  while " d=\"" in contents:
    
    pos = contents.find(" d=\"")
    contents = contents[pos+4:]
    
    pos2 = contents.find("\"")
    path = contents[:pos2]
    contents = contents[pos2+1:]
    
    paths.append(path)

print "Found {} paths".format(len(paths))

meshes = []
concatenate_mesh = None
for (path_no, path) in enumerate(paths):

  print "{}. Path".format(path_no+1)
  n_points = resolution
  svg_mesh = svg_path_to_mesh(path, n_points, None, thickness, invert_normals)
  
  meshes.append(svg_mesh)

  if concatenate_mesh is None:
    concatenate_mesh = svg_mesh
  else:
    concatenate_mesh = mesh.Mesh(np.concatenate([
        svg_mesh.data.copy(),
        concatenate_mesh.data.copy(),
    ]))

concatenate_mesh.save(outfile)
