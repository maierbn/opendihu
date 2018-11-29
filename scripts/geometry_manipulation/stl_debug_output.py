#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import stl
from stl import mesh

import stl_create_rings

def output_points(filename, rankNo, points, size):
  
  triangles = []
  
  for p in points:
    point = np.array([p[0], p[1], p[2]])
    stl_create_rings.create_point_marker(point, triangles, size)

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  outfile = "{}_{}.stl".format(filename, rankNo)
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))



def output_border_points(filename, rankNo, points, size):
  
  triangles = []
  
  for p1 in points:
    for p2 in p1:
      for p3 in p2:
        point = np.array([p3[0], p3[1], p3[2]])
        stl_create_rings.create_point_marker(point, triangles, size)

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
    #for j in range(3):
      #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
      
  #out_mesh.update_normals()

  outfile = "{}.{}.stl".format(filename, rankNo)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

