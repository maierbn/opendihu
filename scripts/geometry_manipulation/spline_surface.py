#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import csv
import collections
import copy
import scipy.spatial
import os
import pickle
import stl_create_mesh   # for standardize_loop
import stl_create_rings

import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import stl
from stl import mesh

from geomdl import NURBS
from geomdl import BSpline
from geomdl import utilities
from geomdl import exchange
from geomdl import operations
from geomdl import fitting
from geomdl import exchange
from geomdl.visualization import VisPlotly
from geomdl.visualization import VisMPL

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

  if len(sys.argv) < 3:
    print("usage: {} <input filename> <output filename>".format(sys.argv[0]))
    sys.exit(0)

  input_filename = sys.argv[1]
  output_filename = sys.argv[2]
  bottom_clip = 70
  top_clip = 250
  n_loops = 12
  
  n_points_u = 10          # x-y direction (along rings)
  n_points_v = n_loops    # z direction

  
  print("input_filename: {}".format(input_filename))
  print("output_filename: {}".format(output_filename))
  print("bottom_clip: {}".format(bottom_clip))
  print("top_clip: {}".format(top_clip))
  print("input_filename: {}".format(input_filename))
  

  write_output_mesh = False
  rings = stl_create_rings.create_rings(input_filename, bottom_clip, top_clip, n_loops, write_output_mesh)
  [loops, lengths] = stl_create_mesh.rings_to_border_points(rings, n_points_u-1)

  loops = stl_create_mesh.border_point_loops_to_list(loops)

  #print("loops: ",loops)

  # duplicate first point
  new_loops = []
  for loop in loops:
    #new_loop = loop + loop + loop + [loop[0]]
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
  [surface, surface2] = operations.split_surface_u(surface1, 0.5585)

  #surface.evaluate_single()
  # evaluate surfaces
  #for surf in [surface, surface0, surface1, surface2, surface_full]:
  for surf in [surface]:
    surf.sample_size = 200
    surf.evaluate()
  
  if False:
    surf = surface
   
    # Prepare points for plotting
    surfpts = np.array(surf.evalpts)

    # Start plotting of the surface and the control points grid
    fig = plt.figure(figsize=(10.67, 8), dpi=96)
    ax = Axes3D(fig)
    ax.set_aspect('equal')

    # Plot surface points
    ax.scatter([x for [x,y,z] in evaluated_points], [y for [x,y,z] in evaluated_points], [z for [x,y,z] in evaluated_points])
    ax.plot_trisurf(surfpts[:, 0], surfpts[:, 1], surfpts[:, 2], color='xkcd:gold', alpha=0.5)

    # Add legend to 3D plot, @ref: https://stackoverflow.com/a/20505720
    surface_prx = matplotlib.lines.Line2D([0], [0], linestyle='none', color='xkcd:gold', marker='^')
    ax.legend([surface_prx],
              ['Surface Plot'],
              numpoints=1)
              
    
    set_axes_equal(ax)

    # Rotate the axes and update the plot
    #for angle in range(0, 360, 10):
    #    ax.view_init(30, angle)
    #    plt.draw()
    #    plt.pause(.001)

    # Display the final 3D plot
    #plt.savefig("out.png")
    #plt.show()

  
  #surface.vis = VisMPL.VisSurfWireframe()
  if False:
    surface_full.vis = VisPlotly.VisSurface()
    surface_full.render()
    surface.vis = VisPlotly.VisSurface()
    surface.render()
  #VisMPL.save_figure_as(surface.vis,"output.png")
  #exchange.export_stl(surface0, 'surface0.stl')
  #exchange.export_stl(surface1, 'surface1.stl')
  #exchange.export_stl(surface2, 'surface2.stl')
  exchange.export_stl(surface, output_filename)
  #exchange.export_stl(surface_full, 'surface_full.stl')
  print("File \"{}\" written.".format(output_filename))
  
