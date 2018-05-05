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
      
    
  if "mooney_rivlin_incompressible_penalty2d_numeric_jacobian" in solution_files[0] \
    or "mooney_rivlin_incompressible_penalty2d_analytic_jacobian" in solution_files[0]:
    import settings_2d
    
    from sympy import *
    from sympy.abc import *
    c0, lambdaValue,tmax,that,lz,ly = symbols('c0,lambdaValue,tmax,that,lz,ly',positive=True)
    
    # (2D) total load constant (traction in reference configuration),  direct solution
    term_direct = solve(2*c0*lambdaValue**4 - tmax/ly*lambdaValue**3 - 2*c0, lambdaValue, simplify=True, positive=True)[3]

    lx_settings = settings_2d.config["FiniteElementMethod"]["physicalExtent"][0]
    ly_settings = settings_2d.config["FiniteElementMethod"]["physicalExtent"][1]
    c0_settings = settings_2d.config["FiniteElementMethod"]["materialParameters"][0]
    c1_settings = settings_2d.config["FiniteElementMethod"]["materialParameters"][1]
    kappa_settings = settings_2d.config["FiniteElementMethod"]["materialParameters"][2]
    tmax_settings = settings_2d.tmax
    dirichlet_bc_settings = settings_2d.config["FiniteElementMethod"]["dirichletBoundaryCondition"]
    
    refconf_volume = lx_settings * ly_settings
    
    analytic_lambda = float(term_direct.subs([(ly,ly_settings), (c0,c0_settings), (tmax,tmax_settings)]).evalf())
    #analytic_lambda = 1.06887593136852
    print "analytic lambda: ",analytic_lambda
    
    # compose reference solution, linear element
    reference_solution = [
      np.array([dirichlet_bc_settings[0], dirichlet_bc_settings[1]]),
      np.array([analytic_lambda*lx_settings, dirichlet_bc_settings[1]]),
      np.array([dirichlet_bc_settings[4], 1./analytic_lambda*ly_settings]),
      np.array([analytic_lambda*lx_settings, 1./analytic_lambda*ly_settings]),
    ]
      
    # calculate volume of reference solution
    p = -reference_solution[0] + reference_solution[3]  # diagonals of quadrilateral
    q = -reference_solution[1] + reference_solution[2]
    reference_volume = 0.5*abs(p[0]*q[1] - p[1]*q[0])   # 2D cross product
    
      
  elif "mooney_rivlin_incompressible_mixed2d_numeric_jacobian_scenario_1" in solution_files[0] \
    or "mooney_rivlin_incompressible_mixed2d_analytic_jacobian_scenario_1" in solution_files[0]:
    import settings_mixed_2d as settings_2d
    
    from sympy import *
    from sympy.abc import *
    c0, lambdaValue,tmax,that,lz,ly = symbols('c0,lambdaValue,tmax,that,lz,ly',positive=True)
    
    # (2D) total load constant (traction in reference configuration),  direct solution
    term_direct = solve(2*c0*lambdaValue**4 - tmax/ly*lambdaValue**3 - 2*c0, lambdaValue, simplify=True, positive=True)[3]

    lx_settings = settings_2d.config["FiniteElementMethod"]["physicalExtent"][0]
    ly_settings = settings_2d.config["FiniteElementMethod"]["physicalExtent"][1]
    c0_settings = settings_2d.config["FiniteElementMethod"]["materialParameters"][0]
    c1_settings = settings_2d.config["FiniteElementMethod"]["materialParameters"][1]
    tmax_settings = settings_2d.tmax
    dirichlet_bc_settings = settings_2d.config["FiniteElementMethod"]["dirichletBoundaryCondition"]
    
    refconf_volume = lx_settings * ly_settings
    
    analytic_lambda = float(term_direct.subs([(ly,ly_settings), (c0,c0_settings), (tmax,tmax_settings)]).evalf())
    #analytic_lambda = 1.06887593136852
    print "analytic lambda: ",analytic_lambda
    
    # compose reference solution, quadratic element
    reference_solution = [
      # bottom row
      np.array([0.0,                             0.0]), 
      np.array([0.5*analytic_lambda*lx_settings, 0.0]),
      np.array([1.0*analytic_lambda*lx_settings, 0.0]),
      # center row
      np.array([0.0,                             0.5/analytic_lambda*ly_settings]), 
      np.array([0.5*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
      np.array([1.0*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
      # top row
      np.array([0.0,                             1.0/analytic_lambda*ly_settings]), 
      np.array([0.5*analytic_lambda*lx_settings, 1.0/analytic_lambda*ly_settings]),
      np.array([1.0*analytic_lambda*lx_settings, 1.0/analytic_lambda*ly_settings]),
    ]
      
    # calculate volume of reference solution
    p = -reference_solution[0] + reference_solution[8]  # diagonals of quadrilateral
    q = -reference_solution[2] + reference_solution[6]
    reference_volume = 0.5*abs(p[0]*q[1] - p[1]*q[0])   # 2D cross product
    
      
  # compare simulation solution to reference solution
  error_absolute_timestep = []
  error_relative_timestep = []
    
  #reference_solution = [
  #  np.array([0.0, 0.0]),
  #  np.array([1.0703, 0.0]),
  #  np.array([0.0, 0.9368]),
  #  np.array([1.0703, 0.9368]),
  #]
  
  simulation_solution = [np.array([x,y]) for (x,y) in zip(x_positions_current, y_positions_current)]
  
  # calculate volume of simulation solution
  if data[0]["basisOrder"] == 1:    # linear geometry
    p = -simulation_solution[0] + simulation_solution[3]  # diagonals of quadrilateral
    q = -simulation_solution[1] + simulation_solution[2]
    solution_volume = 0.5*abs(p[0]*q[1] - p[1]*q[0])   # 2D cross product
  
  elif data[0]["basisOrder"] == 2:   # quadratic geometry
    p = -simulation_solution[0] + simulation_solution[8]  # diagonals of quadrilateral
    q = -simulation_solution[2] + simulation_solution[6]
    solution_volume = 0.5*abs(p[0]*q[1] - p[1]*q[0])   # 2D cross product
  
  print "  volume reference config: ", refconf_volume
  print "  volume reference solution: ", reference_volume
  print "  volume simulated solution: ", solution_volume
  
  error_list = []
  # loop over points
  for (point_simulation_solution,point_reference_solution) in zip(simulation_solution, reference_solution):
        
    # compute difference and error
    difference = point_simulation_solution - point_reference_solution
    error_absolute = np.linalg.norm(difference)
    error_list.append(error_absolute)
    
    print "  point_simulation_solution: {}, point_reference_solution: {}, difference: {}, error: {}".format(point_simulation_solution, point_reference_solution, difference, error_absolute)
    
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
  
  message += ", volume ref.conf., cur.conf.ana., cur.conf.sim: {}, {}, {}".format(refconf_volume, reference_volume, solution_volume)
  
  # print message and write to log file
  print(message)
  with open("log.txt","a+") as f:
    name = solution_files[0]
    while "/" in name:
      name = name[name.find("/")+1:]
      
    f.write("{}: {}\n".format(datetime.datetime.today().strftime('%d.%m.%Y %H:%M:%S'),message))
  
  
sys.exit(0)
