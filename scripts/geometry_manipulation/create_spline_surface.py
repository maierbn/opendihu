#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
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

import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import stl
from stl import mesh

import geomdl
from geomdl import NURBS
from geomdl import BSpline
from geomdl import utilities
from geomdl import exchange
from geomdl import construct
from geomdl import linalg
from geomdl import operations
from geomdl import fitting
from geomdl import exchange

try:
  from geomdl.visualization import VisPlotly
  from geomdl.visualization import VisMPL
except:
  pass

# source: https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_axes_radius(ax, origin, radius):
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    set_axes_radius(ax, origin, radius)

# load stl mesh and evaluate 
if __name__ == "__main__":

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
  
  n_points_u = 10          # x-y direction (along rings)
  n_points_v = n_loops    # z direction

  
  print("input_filename:         {}".format(input_filename))
  print("output_stl_filename:    {}".format(output_stl_filename))
  print("output_pickle_filename: {}".format(output_pickle_filename))
  print("bottom_clip:            {}".format(bottom_clip))
  print("top_clip:               {}".format(top_clip))
  
  debug = False
  
  # try to load stored surface
  try:
    f = open(output_pickle_filename,"rb")
    surface = pickle.load(f)
  except Exception as e:
    
    print("Create rings because pickle file {} does not yet exist.".format(output_pickle_filename))
    # create surface

    write_output_mesh = False
    rings = stl_create_rings.create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh)
    [loops, lengths] = stl_create_mesh.rings_to_boundary_points(rings, n_points_u-1)

    loops = stl_create_mesh.boundary_point_loops_to_list(loops)

    #print("loops: ",loops)

    # duplicate first point
    new_loops = []
    for loop in loops:
      print("  loop has {} points".format(len(loop)))
      #new_loop = loop + [loop[0]]
      new_loop = loop + loop + loop + [loop[0]]
      n_points_u = len(new_loop)
      new_loops.append(new_loop)
    loops = new_loops

    # collect points
    evaluated_points = []
    for i in range(n_points_u):
      for j in range(n_points_v):
        evaluated_points.append(loops[j][i])


    print("n_points_u: ",n_points_u)
    print("n_points_v: ",n_points_v)
    print("{} points".format(len(evaluated_points)))
    #print("evaluated_points: ",evaluated_points)
    surface_full = fitting.approximate_surface(evaluated_points, n_points_u, n_points_v, 3, 2, centripetal=True, ctrlpts_size_u=n_points_u-3, ctrlpts_size_v=n_points_v-3)

    [surface0, surface1] = operations.split_surface_u(surface_full, 0.4)
    [surface, surface2] = operations.split_surface_u(surface1, 0.5555)
    
    pickle_filename = output_pickle_filename
    print("Write pickle file \"{}\"".format(pickle_filename))

    # save surface
    with open(pickle_filename,"wb") as f:
      pickle.dump(surface,f)
  
  # here we have surface
  
  for surf in [surface]:
    surf.sample_size = 100
    surf.evaluate()
  
  #surface.vis = VisMPL.VisSurfWireframe()
  if False:
    surface_full.vis = VisPlotly.VisSurface()
    surface_full.render()
    surface.vis = VisPlotly.VisSurface()
    surface.render()
  print("write stl file \"{}\"".format(output_stl_filename))
  
  exchange.export_stl(surface, output_stl_filename)
  print("File \"{}\" written.".format(output_stl_filename))
  
