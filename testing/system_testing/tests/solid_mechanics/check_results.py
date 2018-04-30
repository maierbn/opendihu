#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Compare output data to analytical solution
# Arguments: [<show_plots (0 or 1)>] [<filenames>]
#

import sys
import numpy as np
import csv
import collections
import copy
from sets import Set
import os
import time
import datetime
import scipy.integrate

import py_reader

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

# extract the files that are py files
solution_condition = lambda filename: ".py" in filename
solution_files = list(np.extract(map(solution_condition, files), files))

print("{} files".format(len(solution_files)))

# load data from files
data = py_reader.load_data(solution_files)

if len(data) == 0:
  print("no data found.")
  sys.exit(0)

dimension = data[0]['dimension']

####################
# 2D
if dimension == 2:

  # parse values of reference and current configuration
  x_positions_current = py_reader.get_values(data[0], "geometry", "x")
  y_positions_current = py_reader.get_values(data[0], "geometry", "y")
  
  x_positions_reference = py_reader.get_values(data[0], "geometryReference", "0")
  y_positions_reference = py_reader.get_values(data[0], "geometryReference", "1")
  
  x_displacements = py_reader.get_values(data[0], "displacements", "0")
  y_displacements = py_reader.get_values(data[0], "displacements", "1")
      
  # compare to reference solution
  error_absolute_timestep = []
  error_relative_timestep = []
  
  reference_solution = [
    np.array([0.0, 0.0]),
    np.array([1.0703, 0.0]),
    np.array([0.0, 0.9368]),
    np.array([1.0703, 0.9368]),
  ]
  
  error_list = []
  for (x,y,point_reference_solution) in zip(x_positions_current, y_positions_current, reference_solution):
    
    # get current node positions
    point_solution = np.array([x,y])
    
    # compute difference and error
    difference = point_solution - point_reference_solution
    error_absolute = np.linalg.norm(difference)
    error_list.append(error_absolute)
    
    #print "point_solution: {}, point_reference_solution: {}, difference: {}, error: {}".format(point_solution, point_reference_solution, difference, error_absolute)
    
  error_absolute_mean = np.mean(error_list)
  error_absolute_median = np.median(error_list)
  
  # determine if test passed
  message = ""
  if error_absolute_mean < 1e-1:
    message = "Test passed: "
  else:
    message = "Test failed: "
  
  message += "error mean/median absolute: {:.2e} / {:.2e}".\
    format(error_absolute_mean, error_absolute_median)
  
  # print message and write to log file
  print(message)
  with open("log.txt","a+") as f:
    name = solution_files[0]
    while "/" in name:
      name = name[name.find("/")+1:]
      
    f.write("{}: {}\n".format(datetime.datetime.today().strftime('%d.%m.%Y %H:%M:%S'),message))
  
  
sys.exit(0)
