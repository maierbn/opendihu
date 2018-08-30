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
import settings_2d

  
# define analytic solution for testing
def analytic_solution_1d(x):
  # u(x) = a + (b-a)*x/l,   a = 
  
  l = settings_1d.physicalExtent
  bc = settings_1d.bc
  a = bc[0]
  b = bc[bc.keys()[-1]]
  value = a + (b-a)*x/l
  return value

def analytic_solution_2d(x1,x2):
  # u(x) = c1*sin(k*pi*x1)*exp(k*pi*x2) + c2*sin(k*pi*x1)*exp(-k*pi*x2)
  
  k = settings_2d.k
  c1 = 1./(2*np.sinh(k*np.pi))
  c2 = -c1
  
  value = c1*np.sin(k*np.pi*x1)*np.exp(k*np.pi*x2) + c2*np.sin(k*np.pi*x1)*np.exp(-k*np.pi*x2)
  #print( "({},{})={}".format(x1,x2,value))
  
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
      test_index_subset = range(len(xdata))
            
      xdata_subset = xdata[test_index_subset]
      ydata_subset = ydata[test_index_subset]
      fdata_subset = fdata[test_index_subset]
      
      # compute the analytic solution
      f_analytic = [analytic_solution_2d(x,y) for (x,y) in zip(xdata_subset,ydata_subset)]
      
  #    if True:  # debugging output and plot
        #print("f_analytic:   ",np.array(f_analytic))
#        print("fdata_subset: ",fdata_subset)


      if "2d_structured_regular_fixed_quadratic" in solution_files[0] and False:

        for (ana,num) in zip(f_analytic,fdata_subset):
          print "ana,num:", ana,num
      if False:  # debugging output and plot
      
        # import needed packages from matplotlib
        if not show_plot:
          import matplotlib as mpl
          mpl.use('Agg')

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import animation
        from matplotlib import cm

        # prepare plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # create mesh
        nEntries = []
        for i in range(dimension):
          nEntries.append(data[0]["nElements"][i] + 1)

        nEntries = nEntries[::-1]   # reverse list
        
        x_positions = py_reader.get_values(data[0], "geometry", "x")
        y_positions = py_reader.get_values(data[0], "geometry", "y")
        
        X = np.reshape(x_positions, nEntries)
        Y = np.reshape(y_positions, nEntries)
      
        print( "nEntries: ", nEntries)
        print( "x_positions: ", len(x_positions))
        print( "X: ",X.size)
        print( "y_positions: ", len(y_positions))
        print( "Y: ",Y.size)
      
        # display data
        #solution_shaped = py_reader.get_values(data[0], "solution", "0")
          
        Z = np.reshape(f_analytic, nEntries)
        print( "x shape: {}, y shape: {}, z shape: {}".format(X.shape, Y.shape, Z.shape))
        
        plt.title("analytical solution")
        ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1)
        plt.show()      
        
      error_absolute = np.array(fdata_subset) - np.array(f_analytic)
      error_absolute_norm = np.linalg.norm(error_absolute) / len(fdata_subset)
      
      error_relative = (np.array(fdata_subset) - np.array(f_analytic)) / fdata_subset
      error_relative_norm = np.linalg.norm(error_relative) / len(fdata_subset)
      
      error_absolute_timestep.append(error_absolute_norm)
      error_relative_timestep.append(error_relative_norm)
      
      #print("error absolute/relative: {}/{}".format(error_absolute_norm,error_relative_norm))
      
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
