#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import stl
from stl import mesh

import stl_create_rings

def output_points(filename, rankNo, points, size):
  
  triangles = []
  
  factor = 1.0
  for p in points:
    point = np.array([p[0], p[1], p[2]])
    stl_create_rings.create_point_marker(point, triangles, size*factor)
    #factor *= 1.1
    #if factor > 3:
    #  factor = 3.0

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  outfile = "{}.{}.{}.stl".format(filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

def output_streamline(filename, rankNo, points, size):
  
  triangles = []
  
  factor = 1.0

  previous_point = None
  for p in points:
    point = np.array([p[0], p[1], p[2]])
    if np.linalg.norm(point) < 1e-3:
      continue
    if previous_point is not None:
      triangles.append([previous_point, point, 0.5*(previous_point+point)])
    previous_point = point
    
    stl_create_rings.create_point_marker(point, triangles, size*factor)
    #factor *= 1.1
    #if factor > 3:
    #  factor = 3.0

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  outfile = "{}.{}.{}.stl".format(filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

def output_streamlines(filename, rankNo, streamlines, size):
  
  triangles = []
  
  factor = 1.0
  
  for points in streamlines:
    previous_point = None
    
    #print("output_streamlines, streamline: {}".format(points))
    for p in points:
      point = np.array([p[0], p[1], p[2]])
      if np.linalg.norm(point) < 1e-3:
        continue
      if previous_point is not None:
        triangles.append([previous_point, point, 0.5*(previous_point+point)])
      previous_point = point
      
      stl_create_rings.create_point_marker(point, triangles, size*factor)
      #factor *= 1.1
      #if factor > 3:
      #  factor = 3.0

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  outfile = "{}.{}.{}.stl".format(filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

def output_border_points(filename, rankNo, points, size):
  
  triangles = []
  print("output_border_points(filename={}, rankNo={}, points: {}, size={})".format(filename, rankNo, len(points), size))
  
  factor = 1.0
  for p1 in points:
    for p2 in p1:
      previous_point = None
      first_point = None
      for p3 in p2:
        point = np.array([p3[0], p3[1], p3[2]])
        if np.linalg.norm(point) < 1e-3:
          continue
        if previous_point is None:
          first_point = point
        else:
          triangles.append([previous_point, point, 0.5*(previous_point+point)])
        previous_point = point
        
        stl_create_rings.create_point_marker(point, triangles, size*factor)
        #factor *= 1.1
        #if factor > 3:
        #  factor = 3.0        
      # close loop
      if previous_point is not None:
        triangles.append([previous_point, first_point, 0.5*(previous_point+first_point)])

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
    #for j in range(3):
      #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
      
  #out_mesh.update_normals()

  outfile = "{}.{}.{}.stl".format(filename[0:2], rankNo, filename[2:])
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

def output_triangles(filename, triangles):
  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  outfile = "{}.stl".format(filename)
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

def output_ghost_elements(filename, rankNo, point_values, n_elements, size):
  
  triangles = []
  
  n_points = (int)(len(point_values)/3)
  
  n_nodes = [n_elements[0]+1, n_elements[1]+1, n_elements[2]+1]
  
  #print("output_ghost_elements, filename {}, n_elements: {}, n points: {}, n_nodes: {}, n_points: {}".format(filename, n_elements, len(point_values), n_nodes, n_points))
  #print("point_values: {}".format(point_values))
  #print("n_points: {}".format(n_points))
  #print("n_elements: {}".format(n_elements))
  
  factor = 1.0
  for z in range(n_elements[2]):
    for y in range(n_elements[1]):
      for x in range(n_elements[0]):
        p = list()
        
        for k in range(2):
          for j in range(2):
            for i in range(2):
              index = (z+k) * n_nodes[0]*n_nodes[1] + (y+j) * n_nodes[0] + (x+i)
              #print("index: {}".format(index))
              p0 = np.array([point_values[index], point_values[n_points+index], point_values[2*n_points+index]])
              p.append(p0)
        #print("x: {}, y: {}, z: {}, points: {}".format(x,y,z,p))
              
        triangles += [
          [p[0],p[3],p[1]],[p[0],p[2],p[3]],  # bottom
          [p[4],p[5],p[7]],[p[4],p[7],p[6]],  # top
          [p[0],p[1],p[5]],[p[0],p[5],p[4]],  # front
          [p[2],p[7],p[3]],[p[2],p[6],p[7]],  # back
          [p[2],p[0],p[4]],[p[2],p[4],p[6]],  # left
          [p[1],p[3],p[7]],[p[1],p[7],p[5]]  # right
        ]
  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  outfile = "{}.{}.{}.stl".format(filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))



