#!../../../../../dependencies/python/install/bin/python3
# -*- coding: utf-8 -*-

import os, sys
import subprocess
import compute_error
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# set global parameters for font sizes
plt.rcParams.update({'font.size': 20})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8

def run(dt_0D, dt_1D, dt_3D, output_name):   
  command = "rm -f out/{output_name}* && mpirun -n 2 ./multiple_fibers ../multiple_fibers_settings.py {dt_0D} {dt_1D} {dt_3D} {output_name}"\
   .format(dt_0D=dt_0D, dt_1D=dt_1D, dt_3D=dt_3D, output_name=output_name)

  try:
    print(command)
    subprocess.check_call(command, shell=True)
  except:
    pass
    
# arguments: <dt_0D> <dt_1D> <dt_3D> <output_name>

if True:
  # reference simulation
  dt_1D = 1e-8                      # timestep width of diffusion
  dt_0D = 1e-8                      # timestep width of ODEs
  dt_3D = 1e-8                      # overall timestep width of splitting
  run(dt_0D, dt_1D, dt_3D, "reference")

  # simulations to evaluate precision
  dt_1D = 1e-6                      # timestep width of diffusion
  dt_0D = 1e-6                      # timestep width of ODEs
  dt_3D = 1e-5                      # overall timestep width of splitting

  # Suggestion by Aaron for 0D: dt_0D=7*1e-5, gives relative error: 1e-4 after t=100

  #####################################
  # evaluate dt_1D
  n = 50
  xlist = np.logspace(-6,-2,n)  # spaced evenly on a log scale.
  ylist = np.zeros(n)

  for i,dt_1D in enumerate(xlist):
    dt_3D = dt_1D
    
    run(dt_0D, dt_1D, dt_3D, "evaluate")
    error = compute_error.compute_error("out/reference", "out/evaluate", "solution")
    ylist[i] = error
    print("dt_1D: {}, error: {}".format(dt_1D,error))

  # print to console and save to file
  with open("dt_1D.csv", "w") as f:
    f.write("dt_1D;error\n")
    for t,error in zip(xlist,ylist):
      print("dt_1D: {}, error: {}".format(t,error))
      f.write("{};{}\n".format(t,error))
      
  # create plot
  plt.figure(figsize=(10,8))
  plt.plot(xlist,ylist,'b+-')
  plt.xlabel('dt_1D')
  plt.ylabel('relative error')
  plt.xscale('log')
  plt.yscale('log')

  plt.savefig("dt_1D.pdf")
  #
  #dt_1D: 1e-06, error: 9.223602548973777e-07
  #dt_1D: 2.782559402207126e-06, error: 6.918122701955344e-07
  #dt_1D: 7.742636826811277e-06, error: 2.3090985353027554e-07
  #dt_1D: 2.1544346900318823e-05, error: 1.3808642879458062e-06
  #dt_1D: 5.994842503189409e-05, error: 5.743946480651069e-06
  #dt_1D: 0.0001668100537200059, error: 1.805607838942555e-05
  #dt_1D: 0.00046415888336127773, error: 5.139642235485197e-05
  #dt_1D: 0.0012915496650148827, error: 0.00013979816121055584
  #dt_1D: 0.003593813663804626, error: 0.0003649962662019944
  #dt_1D: 0.01, error: 0.0009032256016649748
  #
  # -> dt_1D = 5e-4


  #####################################
  # evaluate dt_0D
  dt_1D = 1e-6                      # timestep width of diffusion
  dt_0D = 1e-6                      # timestep width of ODEs
  dt_3D = 1e-5                      # overall timestep width of splitting

  xlist = np.logspace(-6,-2,n)  # spaced evenly on a log scale.
  ylist = np.zeros(n)

  for i,dt_0D in enumerate(xlist):
    dt_3D = dt_0D
    
    run(dt_0D, dt_1D, dt_3D, "evaluate")
    error = compute_error.compute_error("out/reference", "out/evaluate", "solution")
    ylist[i] = error
    print("dt_0D: {}, error: {}".format(dt_0D,error))

  # print to console and save to file
  with open("dt_0D.csv", "w") as f:
    f.write("dt_0D;error\n")
    for t,error in zip(xlist,ylist):
      print("dt_0D: {}, error: {}".format(t,error))
      f.write("{};{}\n".format(t,error))

  # create plot
  plt.figure(figsize=(10,8))
  plt.plot(xlist,ylist,'b+-')
  plt.xlabel('dt_0D')
  plt.ylabel('relative error')
  plt.xscale('log')
  plt.yscale('log')

  plt.savefig("dt_0D.pdf")

  #####################################
  # evaluate dt_3D
  dt_1D = 1e-6                      # timestep width of diffusion
  dt_0D = 1e-6                      # timestep width of ODEs
  dt_3D = 1e-5                      # overall timestep width of splitting

  xlist = np.logspace(-6,-2,n)  # spaced evenly on a log scale.
  ylist = np.zeros(n)

  for i,dt_3D in enumerate(xlist):
    
    run(dt_0D, dt_1D, dt_3D, "evaluate")
    error = compute_error.compute_error("out/reference", "out/evaluate", "solution")
    ylist[i] = error
    print("dt_3D: {}, error: {}".format(dt_3D,error))

  # print to console and save to file
  with open("dt_3D.csv", "w") as f:
    f.write("dt_3D;error\n")
    for t,error in zip(xlist,ylist):
      print("dt_3D: {}, error: {}".format(t,error))
      f.write("{};{}\n".format(t,error))

  # create plot
  plt.figure(figsize=(10,8))
  plt.plot(xlist,ylist,'b+-')
  plt.xlabel('dt_3D')
  plt.ylabel('relative error')
  plt.xscale('log')
  plt.yscale('log')

  plt.savefig("dt_3D.pdf")


  # results: error = 10^-4: dt_1D=10^-3, dt_0D=3*10^-3



#####################################
# evaluate dt_3D for the chosen dt_0D, dt_1D
dt_1D = 1e-3                      # timestep width of diffusion
dt_0D = 3e-3                      # timestep width of ODEs
dt_3D = 1e-3                      # overall timestep width of splitting

xlist = np.logspace(-3,-1,n)  # spaced evenly on a log scale.
ylist = np.zeros(n)

for i,dt_3D in enumerate(xlist):
  
  run(dt_0D, dt_1D, dt_3D, "evaluate")
  error = compute_error.compute_error("out/reference", "out/evaluate", "solution")
  ylist[i] = error
  print("dt_3D: {}, error: {}".format(dt_3D,error))

# print to console and save to file
with open("dt_3Db.csv", "w") as f:
  f.write("dt_3D;error\n")
  for t,error in zip(xlist,ylist):
    print("dt_3D: {}, error: {}".format(t,error))
    f.write("{};{}\n".format(t,error))

# create plot
plt.figure(figsize=(10,8))
plt.plot(xlist,ylist,'b+-')
plt.xlabel('dt_3D')
plt.ylabel('relative error')
plt.xscale('log')
plt.yscale('log')

plt.savefig("dt_3Db.pdf")

