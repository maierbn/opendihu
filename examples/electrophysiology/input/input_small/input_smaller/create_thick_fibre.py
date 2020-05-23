#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script reads in a set of fibres and create a thick cylindrical fibre with given radius
#
# usage:  ./create_thick_fibre.py <fibre_no> [<radius>]

import sys, os
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import collections, patches
import struct
import stl
from stl import mesh
import pickle

fibre_filename = "../build_debug/laplace3d_structured_linear"

with open(fibre_filename, 'rb') as f:
  data = pickle.load(f)


radius = 0.1
fibre_no = 50

if len(sys.argv) < 2:
  print "usage: ./create_thick_fibre.py <fibre_no> [<radius>]"
  sys.exit(0)

if len(sys.argv) >= 2:
  fibre_no = int(sys.argv[1])
  
if len(sys.argv) >= 3:
  radius = float(sys.argv[2])

  
print "radius: {}".format(radius)
print "n fibres: {}".format(len(data))
print "select fibre no {}".format(fibre_no)


fibre_data = data[fibre_no]

triangles = []
previous_point = None
for point in fibre_data:
  if previous_point is not None:
    
    n_steps = 6
    delta_phi = 2*np.pi / n_steps
    for phi in np.linspace(0, 2*np.pi, n_steps+1):
      
      p0 = np.array(point) + np.array([np.cos(phi)*radius, np.sin(phi)*radius, 0.0])
      p1 = np.array(point) + np.array([np.cos(phi+delta_phi)*radius, np.sin(phi+delta_phi)*radius, 0.0])
      p2 = np.array(previous_point) + np.array([np.cos(phi)*radius, np.sin(phi)*radius, 0.0])
      p3 = np.array(previous_point) + np.array([np.cos(phi+delta_phi)*radius, np.sin(phi+delta_phi)*radius, 0.0])
      
      triangles += [[p0, p1, p3], [p0, p3, p2]]
      
  previous_point = point

  
# write output stl meshes
def write_stl(triangles, outfile, description):
  # create output mesh
  n_triangles = len(triangles)

  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
    #for j in range(3):
      #print "set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
      
  out_mesh.update_normals()

  out_mesh.save(outfile) #, mode=stl.Mode.ASCI
  print "saved {} triangles to \"{}\" ({})".format(n_triangles,outfile,description)

write_stl(triangles, "thick_fibre_{:02}.stl".format(fibre_no), "thick fibre")
