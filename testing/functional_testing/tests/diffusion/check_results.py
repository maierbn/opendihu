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
def analytic_solution_2d(x0,x1,t):
  
  def h(x0,x1,t):
    c=settings_2d.c
    return 1./(4*np.pi*c*t)*np.exp(-(x0**2+x1**2) / (4*c*t))
  
  def integrand(y0,y1):
    #print("integrand({},{})={}".format(y0,y1,initial_values_function(y0,y1)*h(x0-y0,x1-y1,t)))
    
    return settings_2d.initial_values_function(y0,y1)*h(x0-y0,x1-y1,t)
  
    
  #(value,error) = scipy.integrate.nquad(integrand, [(0.0, 4.0), (0.0, 3.0)], opts={"limit":5})
  (value,error) = scipy.integrate.nquad(integrand, [(-np.infty, np.infty), (-np.infty, np.infty)])
  
  #print("analytic_solution({},{},{})={} (integration error {})".format(x0,x1,t,value,error))
  return value
  
# define analytic solution for testing
def analytic_solution_1d(x,t):
  # note: this analytic solution would be correct for a non-bounded domain
  
  def h(x,t):
    c=settings_1d.c
    return 1./np.sqrt(4*np.pi*c*t)*np.exp(-x**2 / (4*c*t))
  
  def integrand(y):
    return settings_1d.initial_values_function(y)*h(x-y,t)
    
  (value,error) = scipy.integrate.quad(integrand, -np.infty, np.infty)
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
  
    y_analytic = [analytic_solution_1d(x,t) for x in xdata]
    
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
      
    f.write("{} [{}]: {}\n".format(datetime.datetime.today().strftime('%d.%m.%Y %H:%M:%S'),name,message))
  
  # prepare animation
  fig = plt.figure()

  margin = abs(max_value - min_value) * 0.1
  ax = plt.axes(xlim=(min_x, max_x), ylim=(min_value - margin, max_value + margin))
  line_numerical, = ax.plot([], [], 'o-', lw=2, label='numerical')
  line_analytical, = ax.plot([], [], '+-', lw=1, label='analytical')
  text = plt.figtext(0.2,0.7,"timestep",size=15)
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  plt.legend()
  
  def init():
    line_numerical.set_data([], [])
    line_analytical.set_data([], [])
    text.set_text("timestep {}".format(0))
    return line_numerical,

  def animate(i):
    
    if 'timeStepNo' in data[i]:
      timestep = data[i]['timeStepNo']
    if 'currentTime' in data[i]:
      current_time = data[i]['currentTime']
      t = current_time
      
    # display data
    xdata = py_reader.get_values(data[i], "geometry", "x")
    ydata_numerical = py_reader.get_values(data[i], "solution", "0")
    ydata_analytical = [analytic_solution_1d(x,t) for x in xdata]
    line_numerical.set_data(xdata, ydata_numerical)
    line_analytical.set_data(xdata, ydata_analytical)
    
    # display timestep
    max_timestep = len(data)-1
      
    text.set_text("timestep {}/{}\nt = {}\nabs. error: {}".format(timestep, max_timestep, current_time, error_absolute_timestep[i]))
    
    return line_numerical,
    
  interval = 5000.0 / len(data)
        
  # call the animator.  blit=True means only re-draw the parts that have changed.
  anim = animation.FuncAnimation(fig, animate, init_func=init,
             frames=len(data), interval=interval, blit=False)

  anim.save("numerical_analytical.mp4")
  if show_plot:
    plt.show()
  
    
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
        
      f.write("{} [{}]: {}\n".format(datetime.datetime.today().strftime('%d.%m.%Y %H:%M:%S'),name,message))
    
  # prepare animation
  fig = plt.figure()

  margin = abs(max_value - min_value) * 0.1
  ax = fig.add_subplot(111, projection='3d', xlim=(min_x, max_x), ylim=(min_y, max_y), zlim=(min_value-margin, max_value+margin))
  
  surface_numerical = ax.plot_surface([], [], [], cmap=cm.coolwarm, linewidth=2, label='numerical',rstride=1,cstride=1)
  surface_analytical = ax.plot_surface([], [], [], cmap=cm.coolwarm, linewidth=1, label='analytical',rstride=1,cstride=1)
  
  text = plt.figtext(0.2,0.7,"timestep",size=15)
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  
  # create mesh
  if data[0]["meshType"] == "StructuredRegularFixed" or data[0]["meshType"] == "RegularFixed":
    
    print("basisfunction: [{}], basisOrder: [{}]".format(data[0]["basisFunction"], data[0]["basisOrder"]))
    
    if data[0]["basisFunction"] == "Lagrange":
      nEntries = dimension * [0]
      for i in range(dimension):
        nEntries[i] = data[0]["basisOrder"] * data[0]["nElements"][i] + 1
    
    x_positions = py_reader.get_values(data[0], "geometry", "x")
    y_positions = py_reader.get_values(data[0], "geometry", "y")
    
    nEntries = nEntries[::-1]   # reverse list
    
    X = np.reshape(x_positions, nEntries)
    Y = np.reshape(y_positions, nEntries)
    
    #print("x_positions shape: {}".format(len(x_positions)))
    
  elif data[0]["meshType"] == "StructuredDeformable":
    pass
  
  def init():
    surface_numerical.set_verts([])
    surface_analytical.set_verts([])
    text.set_text("timestep {}".format(0))
    return surface_numerical,

  def animate(i):
    
    # extract timestep and time
    if 'timeStepNo' in data[i]:
      timestep = data[i]['timeStepNo']
    if 'currentTime' in data[i]:
      current_time = data[i]['currentTime']
      t = current_time
    
    # clear old figure
    ax.clear()
    
    # get data for time step
    x_positions = py_reader.get_values(data[i], "geometry", "x")
    y_positions = py_reader.get_values(data[i], "geometry", "y")
    
    xdata = np.reshape(x_positions, nEntries)
    ydata = np.reshape(y_positions, nEntries)
    
    fdata_numerical = py_reader.get_values(data[i], "solution", "0")
    
    # compute analytical solution, this takes very long
    fdata_analytical = [analytic_solution_2d(x,y,t) for y in y_positions for x in x_positions]
    
    # reshape arrays
    fdata_numerical = np.reshape(fdata_numerical, nEntries)
    fdata_analytical = np.reshape(fdata_analytical, nEntries)
    
    #print("xdata: ",xdata)
    #print("ydata: ",ydata)
    #print("fdata_numerical: ",fdata_numerical)
    #print("x shape: {}, y shape: {}, z shape: {}".format(xdata.shape, ydata.shape, fdata_numerical.shape))
    
    # plot 2D functions
    surface_numerical = ax.plot_surface(xdata, ydata, fdata_numerical, cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1)
    surface_analytical = ax.plot_wireframe(xdata, ydata, fdata_analytical, cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1)
    ax.set_zlim(min_value-margin, max_value+margin)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # display timestep
    max_timestep = len(data)-1
    text.set_text("timestep {}/{}\nt = {}\nabs. error: {}".format(timestep, max_timestep, current_time, error_absolute_timestep[i]))
    
    return text,
    
  interval = 5000.0 / len(data)
        
  # call the animator.  blit=True means only re-draw the parts that have changed.
  anim = animation.FuncAnimation(fig, animate, init_func=init,
             frames=len(data),
             interval=interval, blit=False)

  anim.save("numerical_analytical.mp4")
  if show_plot:
    plt.show()
  
  
  
  
sys.exit(0)
