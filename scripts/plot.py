#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to visualize python output files.
#

import sys
import numpy as np

import csv
import collections
import copy
from sets import Set
import os
import time
import pickle
import py_reader    # reader utility for opendihu *.py files

files = ""

show_plot = True
if len(sys.argv) > 1:
  try:
    show_plot = int(sys.argv[1])
    files = sys.argv[2:]
  except:
    files = sys.argv[1:]
else:
  # get all input data in current directory
  ls = os.listdir(".")

  # sort files by number in file name
  files = sorted(ls)


# import needed packages from matplotlib
if not show_plot:
  import matplotlib as mpl
  mpl.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import Polygon

# extract the files that are npy files
solution_condition = lambda filename: "solution.npy" in filename
solution_shaped_condition = lambda filename: "solution_shaped.npy" in filename
solution_py_condition = lambda filename: ".py" in filename

solution_files = list(np.extract(map(solution_condition, files), files))
solution_shaped_files = list(np.extract(map(solution_shaped_condition, files), files))
solution_py_files = list(np.extract(map(solution_py_condition, files), files))


# sort files by number in file name
solution_py_files = sorted(solution_py_files)

print "{} files".format(len(solution_py_files))
print solution_py_files[0:min(10,len(solution_py_files))]

data = py_reader.load_data(solution_py_files)

if len(data) == 0:
  print "no data found."
  sys.exit(0)

dimension = data[0]['dimension']

####################
# 1D
if dimension == 1:
  
  min_value, max_value = py_reader.get_min_max(data, "solution", "0")
  min_x, max_x = py_reader.get_min_max(data, "geometry", "x")
  
  print "value range: [{}, {}]".format(min_value, max_value)
  
  # prepare plot
  fig = plt.figure()

  margin = abs(max_value - min_value) * 0.1
  ax = plt.axes(xlim=(min_x, max_x), ylim=(min_value - margin, max_value + margin))
  line, = ax.plot([], [], 'o-', color=(1.0,0.9,0.8), lw=2)
  text = plt.figtext(0.15,0.85,"timestep",size=20)
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  
  def init():
    line.set_data([], [])
    text.set_text("timestep {}".format(0))
    return line,

  def animate(i):
    # display data
    xdata = py_reader.get_values(data[i], "geometry", "x")
    ydata = py_reader.get_values(data[i], "solution", "0")
    line.set_data(xdata, ydata)
    
    # display timestep
    if 'timeStepNo' in data[i]:
      timestep = data[i]['timeStepNo']
    if 'currentTime' in data[i]:
      current_time = data[i]['currentTime']
      
    max_timestep = len(data)-1
      
    text.set_text("timestep {}/{}, t = {}".format(timestep, max_timestep, current_time))
    
    return line,
    
  interval = 5000.0 / len(data)
        
  if len(data) == 1:
    init()
    animate(0)
    plt.savefig("fig.pdf")
    
  else:
    # create animation
    anim = animation.FuncAnimation(fig, animate, init_func=init,
               frames=len(data), interval=interval, blit=False)

    anim.save("anim.mp4")
    
    # create plot with first and last dataset
    # plot first dataset
    line, = ax.plot([], [], 'o-', color=(0.8,1.0,0.9), lw=2, label="t=0")
    line0, = animate(0)
    
    # plot last dataset
    i = len(data)-1
    if 'timeStepNo' in data[i]:
      timestep = data[i]['timeStepNo']
    if 'currentTime' in data[i]:
      current_time = data[i]['currentTime']
      
    line, = ax.plot([], [], 'o-b', lw=2, label="t={}".format(current_time))
    line1, = animate(i)
    
    max_timestep = len(data)-1
    text.set_text("timesteps 0 and {}".format(timestep))
    
    ax = plt.gca()
    #ax.add_line(line0)
    #ax.add_line(line1)
    plt.legend()
      
    plt.savefig("fig.pdf")
    
  if show_plot:
    plt.show()
  
####################
# 2D
if dimension == 2:
  
  field_variable_names = py_reader.get_field_variable_names(data[0])
  
  # classical 2D scalar field variables, like in Laplace eq.
  if "solution" in field_variable_names:
    
    debug = False
    
    min_value, max_value = py_reader.get_min_max(data, "solution", "0")
    min_x, max_x = py_reader.get_min_max(data, "geometry", "x")
    min_y, max_y = py_reader.get_min_max(data, "geometry", "y")
    
    print "value range: [{}, {}]".format(min_value, max_value)
    
    # prepare plot
    fig = plt.figure()

    margin = abs(max_value - min_value) * 0.1
    ax = fig.add_subplot(111, projection='3d', xlim=(min_x, max_x), ylim=(min_y, max_y), zlim=(min_value-margin, max_value+margin))
    
    surface = ax.plot_surface([], [], [], cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1)
    text = plt.figtext(0.15,0.85,"timestep",size=20)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # create mesh
  
    if data[0]["basisFunction"] == "Lagrange":
      n_average_nodes_1D_per_element = data[0]["basisOrder"]
    
    elif data[0]["basisFunction"] == "Hermite":
      n_average_nodes_1D_per_element = 1
    
    if data[0]["meshType"] == "StructuredRegularFixed" or data[0]["meshType"] == "RegularFixed" or data[0]["meshType"] == "StructuredDeformable":
      
      if debug:
        print "basisfunction: [{}], basisOrder: [{}]".format(data[0]["basisFunction"], data[0]["basisOrder"])
      
      nEntries = []
      for i in range(dimension):
        nEntries.append(n_average_nodes_1D_per_element * data[0]["nElements"][i] + 1)

      nEntries = nEntries[::-1]   # reverse list
      
      x_positions = py_reader.get_values(data[0], "geometry", "x")
      y_positions = py_reader.get_values(data[0], "geometry", "y")
      
      X = np.reshape(x_positions, nEntries)
      Y = np.reshape(y_positions, nEntries)
      
      if debug:
        print "nEntries: ", nEntries
        print "x_positions: ", x_positions
        print "X: ",X
        print "y_positions: ", y_positions
        print "Y: ",Y
      
      #print "x_positions shape: {}".format(len(x_positions))
      
    elif data[0]["meshType"] == "UnstructuredDeformable":
      if not data[0]["onlyNodalValues"]:
        print "Error: onlyNodalValues is False, set to True in OutputWriter config!"
      
      x_positions = py_reader.get_values(data[0], "geometry", "x")
      y_positions = py_reader.get_values(data[0], "geometry", "y")
      X = x_positions
      Y = y_positions
      
      triangles = []
      for elemental_dofs in data[0]["elementalDofs"]:
        
        # for linear Lagrange and Hermite
        if n_average_nodes_1D_per_element == 1:
          triangles.append([elemental_dofs[0], elemental_dofs[1], elemental_dofs[3]])
          triangles.append([elemental_dofs[0], elemental_dofs[3], elemental_dofs[2]])
        else:  # for quadratic Lagrange
          triangles.append([elemental_dofs[0], elemental_dofs[1], elemental_dofs[4]])
          triangles.append([elemental_dofs[0], elemental_dofs[4], elemental_dofs[3]])
          triangles.append([elemental_dofs[4], elemental_dofs[1], elemental_dofs[2]])
          triangles.append([elemental_dofs[4], elemental_dofs[2], elemental_dofs[5]])
          triangles.append([elemental_dofs[4], elemental_dofs[5], elemental_dofs[6]])
          triangles.append([elemental_dofs[4], elemental_dofs[6], elemental_dofs[7]])
          triangles.append([elemental_dofs[4], elemental_dofs[6], elemental_dofs[3]])
          triangles.append([elemental_dofs[4], elemental_dofs[7], elemental_dofs[6]])
    
    def animate(i):
      ax.clear()
      
      # display data
      solution_shaped = py_reader.get_values(data[i], "solution", "0")
      
      try:
        Z = np.reshape(solution_shaped, nEntries)
      except:
        Z = solution_shaped
      
      if debug:
        try:
          print "x shape: {}, y shape: {}, z shape: {}".format(X.shape, Y.shape, Z.shape)
        except:
          pass
      
      # for unstructured grid use plot_trisurf
      if data[0]["meshType"] == "UnstructuredDeformable":
        plot = ax.plot_trisurf(X, Y, triangles, Z, cmap=cm.coolwarm, linewidth=1)
        
      # for structured grids use plot_surface
      else:
        plot = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1)
        
      ax.set_zlim(min_value-margin, max_value+margin)
      ax.set_xlabel('X')
      ax.set_ylabel('Y')
      ax.set_zlabel('Z')
      
      # display timestep
      if 'timeStepNo' in data[i]:
        timestep = data[i]['timeStepNo']
        max_timestep = data[-1]['timeStepNo']
      if 'currentTime' in data[i]:
        current_time = data[i]['currentTime']
        
      if timestep == -1:
        text.set_text("t = {}".format(current_time))
      else:
        text.set_text("timestep {}/{}, t = {}".format(timestep, max_timestep, current_time))
      
      return plot,
      
    interval = 5000.0 / len(data)
          
    if len(data) == 1:
      animate(0)
      plt.savefig("fig.pdf")
      
    else:# create animation
      anim = animation.FuncAnimation(fig, animate,
                 frames=len(data), interval=interval, blit=False)

      anim.save("anim.mp4")
      
      # create plot with first and last dataset
      fig = plt.figure(figsize=(5,10))
      ax = fig.add_subplot(211, projection='3d', xlim=(min_x, max_x), ylim=(min_y, max_y), zlim=(min_value-margin, max_value+margin))
      
      # plot first dataset
      plot0, = animate(0)
      
      # plot last dataset
      ax = fig.add_subplot(212, projection='3d', xlim=(min_x, max_x), ylim=(min_y, max_y), zlim=(min_value-margin, max_value+margin))
    
      i = len(data)-1
      if 'timeStepNo' in data[i]:
        timestep = data[i]['timeStepNo']
      if 'currentTime' in data[i]:
        current_time = data[i]['currentTime']
        
      plot1, = animate(i)
      
      text = plt.figtext(0.15,0.85,"timestep",size=20)
      text.set_text("timesteps 0 and {}".format(timestep))
      
      ax = plt.gca()
      #ax.add_line(line0)
      #ax.add_line(line1)
      plt.legend()
        
      plt.tight_layout()
      plt.savefig("fig.pdf")
      
    if show_plot:
      plt.show()
  
  elif "displacements" in field_variable_names:
  # continuum mechanics
    
    debug = False
    if debug:
      # print data
      import pprint 
      pp = pprint.PrettyPrinter()
      pp.pprint(data[0])
    
      print "field_variable_names: ", field_variable_names
    
    min_x_current, max_x_current = py_reader.get_min_max(data, "geometry", "x")
    min_y_current, max_y_current = py_reader.get_min_max(data, "geometry", "y")
    
    min_x_reference, max_x_reference = py_reader.get_min_max(data, "geometryReference", "0")
    min_y_reference, max_y_reference = py_reader.get_min_max(data, "geometryReference", "1")
    
    min_x = min(min_x_current, min_x_reference)
    min_y = min(min_y_current, min_y_reference)
    max_x = max(max_x_current, max_x_reference)
    max_y = max(max_y_current, max_y_reference)
  
    # prepare plot
    fig = plt.figure()

    if debug:
      print "min_x: ", min_x
      print "min_y: ", min_y
      print "max_x: ", max_x
      print "max_y: ", max_y

    # create plot with 30% margins
    margin_x = abs(max_x - min_x) * 0.3
    margin_y = abs(max_y - min_y) * 0.3
    ax = fig.add_subplot(111, xlim=(min_x-margin_x, max_x+margin_x), ylim=(min_y-margin_y, max_y+margin_y))
    text = plt.figtext(0.15,0.85,"timestep",size=20)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
  
    def init():
      pass
        
    def animate(frame_no):
      ax.clear()
    
      dataset = data[frame_no]
      
      # create mesh / polygons
      polygons = []
      if dataset["meshType"] == "StructuredRegularFixed" or dataset["meshType"] == "StructuredDeformable":
        
        if debug:
          print "basisfunction: [{}], basisOrder: [{}]".format(dataset["basisFunction"], dataset["basisOrder"])
        
        if dataset["basisFunction"] == "Lagrange":
          nEntries = dimension * [0]
          for i in range(dimension):
            nEntries[i] = dataset["basisOrder"] * dataset["nElements"][i] + 1
            
        elif dataset["basisFunction"] == "Hermite":
          nEntries = dimension * [0]
          for i in range(dimension):
            nEntries[i] = dataset["nElements"][i] + 1
        
        nEntries = nEntries[::-1]   # reverse list
      
      def create_mesh(x_positions, y_positions, nEntries, **kwargs):
        """ create Polygon objects (quads) that can be plotted """
        polygons = []
        
        X = np.reshape(list(x_positions), nEntries)
        Y = np.reshape(list(y_positions), nEntries)
        
        if debug:
          print "nEntries: ", nEntries
          print "x_positions: ", x_positions
          print "X: ",X
          print "y_positions: ", y_positions
          print "Y: ",Y
          print nEntries
        
        # loop over elements
        for ely in range(nEntries[0]-1):
          for elx in range(nEntries[1]-1):
            point0 = np.array([X[ely][elx], Y[ely][elx]])
            point1 = np.array([X[ely][elx+1], Y[ely][elx+1]])
            point2 = np.array([X[ely+1][elx], Y[ely+1][elx]])
            point3 = np.array([X[ely+1][elx+1], Y[ely+1][elx+1]])
            
            if debug:
              print "polygon (0,1,2,3) = ({},{},{},{})".format(point0, point1, point2, point3)
        
            polygon = Polygon([point0, point1, point3, point2], **kwargs)  # dummy data for xs,ys
            polygons.append(polygon)
            
        return polygons
      
      # parse values of reference and current configuration
      x_positions_current = py_reader.get_values(dataset, "geometry", "x")
      y_positions_current = py_reader.get_values(dataset, "geometry", "y")
      
      x_positions_reference = py_reader.get_values(dataset, "geometryReference", "0")
      y_positions_reference = py_reader.get_values(dataset, "geometryReference", "1")
        
      # create meshes for current and reference configuration
      polygons_current = []
      polygons_current = create_mesh(x_positions_current, y_positions_current, nEntries, fill=False, edgecolor=(0.2, 0.2, 0.8), linewidth=3, label="current configuration")
      polygons_reference = []
      polygons_reference = create_mesh(x_positions_reference, y_positions_reference, nEntries, fill=False, edgecolor=(0.8,0.8,0.8), linewidth=3, label="reference configuration")
      
      
      # add all polygons to axis
      for polygon in polygons_reference + polygons_current:
        ax.add_patch(polygon)
        
      ax.set_aspect('equal','datalim')

      # add legend, every legend entry (current configuration or reference configuration) only once
      handles, labels = ax.get_legend_handles_labels()
      new_handles = []
      new_labels = []
      for (handle,label) in zip(handles,labels):
        if label not in new_labels:
          new_handles.append(handle)
          new_labels.append(label)
      ax.legend(new_handles, new_labels, loc='best')
      
     
      # display timestep
      if 'timeStepNo' in dataset:
        timestep = dataset['timeStepNo']
        max_timestep = data[-1]['timeStepNo']
      if 'currentTime' in dataset:
        current_time = dataset['currentTime']
        
      t = "{}. timestep {}/{}, t = {}".format(frame_no, timestep, max_timestep, current_time)
      text.set_text(t)
      #text = plt.figtext(0.15, 0.85, t, size=20)
    
    interval = 5000.0 / len(data)
          
    if len(data) == 1:
      init()
      animate(0)
      plt.savefig("fig.pdf")
      
      if show_plot:
        plt.show()
      
    else:  # create animation
      
      anim = animation.FuncAnimation(fig, animate, init_func=init,
                 frames=len(data), interval=interval, blit=False)
        
      if show_plot:
        plt.show()
      else:
        anim.save("anim.mp4")
        plt.savefig("fig.pdf")
      
sys.exit(0)
