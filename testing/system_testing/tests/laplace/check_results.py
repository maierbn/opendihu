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
import settings_1d

  
# define analytic solution for testing
def analytic_solution_1d(x):
  # note: this analytic solution would be correct for a non-bounded domain
  bc = settings_1d.bc
  n = settings_1d.n
  physicalExtent = settings_1d.physicalExtent
  b = bc[0]; #a * 0 + b
  a = (bc[n] - b) / physicalExtent # ax + b = bc[n], a = (bc[n] - b) / x[n]
  value = a * x + b
  return value
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

# extract the files that are npy files
solution_condition = lambda filename: ".py" in filename
solution_files = list(np.extract(map(solution_condition, files), files))

print("{} files".format(len(solution_files)))

data = py_reader.load_data(solution_files)

if len(data) == 0:
  print("no data found.")
  sys.exit(0)

dimension = data[0]['dimension']

####################
# 1D
if dimension == 1:
  
  min_value, max_value = py_reader.get_min_max(data, "solution", "0")
  min_x, max_x = py_reader.get_min_max(data, "geometry", "x")
  
  # check if solution diverged
  diverged_tolerance = 1e100
  if min_value < -diverged_tolerance or max_value > diverged_tolerance:
    print("Test failed: solution diverged. Value range: [{}, {}]".format(min_value, max_value))
    sys.exit(0)
  
  # compare to reference solution
  error_absolute_timestep = []
  error_relative_timestep = []
  
  for dataset in data:
    timestep_no = dataset['timeStepNo']
    t = dataset['currentTime']
    xdata = py_reader.get_values(dataset, "geometry", "x")
    ydata = py_reader.get_values(dataset, "solution", "0")
  
    y_analytic = [analytic_solution_1d(x) for x in xdata]
    
    error_absolute = np.array(ydata) - np.array(y_analytic)
    error_absolute_norm = np.linalg.norm(error_absolute) / len(ydata)
    
    error_relative = (np.array(ydata) - np.array(y_analytic)) / ydata
    error_relative_norm = np.linalg.norm(error_relative) / len(ydata)
    
    error_absolute_timestep.append(error_absolute_norm)
    error_relative_timestep.append(error_relative_norm)
    
  # reduce error values for all timesteps
  error_absolute_mean = np.mean(error_absolute_timestep)
  error_absolute_median = np.median(error_absolute_timestep)
  error_relative_mean = np.mean(error_relative_timestep)
  error_relative_median = np.median(error_relative_timestep)
  
  # determine if test passed
  message = ""
  if error_absolute_mean < 1e-1:
    message = "Test passed: "
  else:
    message = "Test failed: "
  
  message += "error mean/median absolute: {:.2e} / {:.2e}, relative: {:.2e} / {:.2e}".\
    format(error_absolute_mean, error_absolute_median, error_relative_mean, error_relative_median)
  
  # print message and write to log file
  print(message)
  with open("log.txt","a+") as f:
    name = solution_files[0]
    while "/" in name:
      name = name[name.find("/")+1:]
      
    f.write("{}: {}\n".format(datetime.datetime.today().strftime('%d.%m.%Y %H:%M:%S'),message))
  
    
####################
# 2D
if dimension == 2:
  
  min_value, max_value = py_reader.get_min_max(data, "solution", "0")
  min_x, max_x = py_reader.get_min_max(data, "geometry", "x")
  min_y, max_y = py_reader.get_min_max(data, "geometry", "y")
  
  print("value range: [{}, {}]".format(min_value, max_value))
  
  # check if solution diverged
  diverged_tolerance = 1e100
  if min_value < -diverged_tolerance or max_value > diverged_tolerance:
    print("Test failed: solution diverged. Value range: [{}, {}]".format(min_value, max_value))
    sys.exit(0)
  
  if True:
    
    # compare to reference solution
    error_absolute_timestep = []
    error_relative_timestep = []
    
    # loop over datasets for timesteps
    for dataset in data:
      # extract time step no and time
      timestep_no = dataset['timeStepNo']
      t = dataset['currentTime']
      
      # extract data
      xdata = np.array(py_reader.get_values(dataset, "geometry", "x"))
      ydata = np.array(py_reader.get_values(dataset, "geometry", "y"))
      fdata = np.array(py_reader.get_values(dataset, "solution", "0"))
      
      # pick a subset of the data to be used for testing
      test_index_subset = []
      for i in range(10):
        test_index_subset.append(int(len(xdata)*i/10.))
            
      xdata_subset = xdata[test_index_subset]
      ydata_subset = ydata[test_index_subset]
      fdata_subset = fdata[test_index_subset]
      
      # compute the analytic solution
      f_analytic = [analytic_solution_2d(x,y,t) for (x,y) in zip(xdata_subset,ydata_subset)]
      
      #print("f_analytic:   ",f_analytic)
      #print("fdata_subset: ",fdata_subset)
      
      error_absolute = np.array(fdata_subset) - np.array(f_analytic)
      error_absolute_norm = np.linalg.norm(error_absolute) / len(fdata_subset)
      
      error_relative = (np.array(fdata_subset) - np.array(f_analytic)) / fdata_subset
      error_relative_norm = np.linalg.norm(error_relative) / len(fdata_subset)
      
      error_absolute_timestep.append(error_absolute_norm)
      error_relative_timestep.append(error_relative_norm)
      
      #print("t={}, error absolute/relative: {}/{}".format(t,error_absolute_norm,error_relative_norm))
      
    # reduce error values for all timesteps
    error_absolute_mean = np.mean(error_absolute_timestep)
    error_absolute_median = np.median(error_absolute_timestep)
    error_relative_mean = np.mean(error_relative_timestep)
    error_relative_median = np.median(error_relative_timestep)
    
    # determine if test passed
    message = ""
    if error_absolute_mean < 1e-1:
      message = "Test passed: "
    else:
      message = "Test failed: "
    
    message += "error mean/median absolute: {:.2e} / {:.2e}, relative: {:.2e} / {:.2e}".\
      format(error_absolute_mean, error_absolute_median, error_relative_mean, error_relative_median)
    
    # print message and write to log file
    print(message)
    with open("log.txt","a+") as f:
      name = solution_files[0]
      while "/" in name:
        name = name[name.find("/")+1:]
        
      f.write("{}: {}\n".format(datetime.datetime.today().strftime('%d.%m.%Y %H:%M:%S'),message))
  
  
  
sys.exit(0)
