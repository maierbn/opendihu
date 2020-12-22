#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This script needs a STL file of a tubular muscle surface as input and then creates a NURBS surface approximation.
# Sample the input STL mesh at 12 z levels and create a spline surface from these points, save it as stl and pickle file.

import sys
import numpy as np
import csv
import collections
import copy
import scipy.spatial
import os
import pickle
import stl_create_rings # for create_rings
import stl_create_mesh   # for rings_to_boundary_points

import stl
from stl import mesh

import geomdl.operations
import geomdl.fitting
import geomdl.exchange

try:
  from geomdl.visualization import VisPlotly
  from geomdl.visualization import VisMPL
except:
  pass

# load stl mesh and evaluate 
if __name__ == "__main__":

  # parse command arguments
  if len(sys.argv) < 6:
    print("usage: {} <input filename> <output stl filename> <output pickle filename> [<bottom clip> <top clip>]".format(sys.argv[0]))
    sys.exit(0)

  input_filename = sys.argv[1]
  output_stl_filename = sys.argv[2]
  output_pickle_filename = sys.argv[3]
  
  bottom_clip = 70
  top_clip = 250
  if len(sys.argv) == 6:
    bottom_clip = float(sys.argv[4])
    top_clip = float(sys.argv[5])
  n_loops = 12
  
  # define number of control points
  n_points_u = 10          # x-y direction (along rings)
  n_points_v = n_loops    # z direction

  # output parameters
  print("input_filename:         {}".format(input_filename))
  print("output_stl_filename:    {}".format(output_stl_filename))
  print("output_pickle_filename: {}".format(output_pickle_filename))
  print("bottom_clip:            {}".format(bottom_clip))
  print("top_clip:               {}".format(top_clip))
  
  # try to load stored surface as pickle file
  try:
    f = open(output_pickle_filename,"rb")
    surface = pickle.load(f)
  except Exception as e:
    
    # if the pickle file does not exist, create it
    print("Create rings because pickle file {} does not yet exist.".format(output_pickle_filename))
    
    # sample surface at equidistant z intervals using the stl_create_rings script
    write_output_mesh = False
    rings = stl_create_rings.create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh)
    [loops, lengths] = stl_create_mesh.rings_to_boundary_points(rings, n_points_u-1)

    loops = stl_create_mesh.boundary_point_loops_to_list(loops)

    # duplicate the rings: 3x given ring + first point
    new_loops = []
    for loop in loops:
      print("  loop has {} points".format(len(loop)))
      #new_loop = loop + [loop[0]]
      new_loop = loop + loop + loop + [loop[0]]
      n_points_u = len(new_loop)
      new_loops.append(new_loop)
    loops = new_loops

    # collect points on every ring
    evaluated_points = []
    for i in range(n_points_u):
      for j in range(n_points_v):
        evaluated_points.append(loops[j][i])

    print("n_points_u: ",n_points_u)
    print("n_points_v: ",n_points_v)
    print("{} points".format(len(evaluated_points)))
    #print("evaluated_points: ",evaluated_points)
    
    # call NURBS surface fitting algorithm of the geomdl.fitting library
    surface_full = geomdl.fitting.approximate_surface(evaluated_points, n_points_u, n_points_v, 3, 2, centripetal=True, ctrlpts_size_u=n_points_u-3, ctrlpts_size_v=n_points_v-3)

    # cut out the inner part that represents the actual surface
    [surface0, surface1] = geomdl.operations.split_surface_u(surface_full, 0.4)
    [surface, surface2] = geomdl.operations.split_surface_u(surface1, 0.5555)

    # save surface to a pickle file    
    pickle_filename = output_pickle_filename
    print("Write pickle file \"{}\"".format(pickle_filename))

    with open(pickle_filename,"wb") as f:
      pickle.dump(surface,f)
  
  # here we have the surface, either loaded from an existing pickle file or newly created

  # triangulate the surface for output
  for surf in [surface]:
    surf.sample_size = 100
    surf.evaluate()
  
  # plot the result
  if False:
    surface_full.vis = VisPlotly.VisSurface()
    surface_full.render()
    surface.vis = VisPlotly.VisSurface()
    #surface.vis = VisMPL.VisSurfWireframe()
    surface.render()
    
  # export as STL file
  print("write stl file \"{}\"".format(output_stl_filename))
  geomdl.exchange.export_stl(surface, output_stl_filename)
  print("File \"{}\" written.".format(output_stl_filename))
  
