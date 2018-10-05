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

def run(n_processes, n_processes_per_fiber, n_fibers, n_nodes_per_fiber, scenario_name):  
  n_nodes = int(np.ceil(n_processes/24))
# arguments: <n_processes_per_fiber> <n_fibers> <n_nodes_per_fiber> <scenario_name> 
  command = "mpirun -n {} ./cuboid ../cuboid_settings.py {} {} {} {}   # {} nodes"\
   .format(n_processes, n_processes_per_fiber, n_fibers, n_nodes_per_fiber, scenario_name, n_nodes)

  try:
    print(command)
    #subprocess.check_call(command, shell=True)
  except:
    pass
    
##########################
# strong scaling
print("strong scaling")
n_nodes_per_fiber = 1000
n_fibers = 100
n_processes_per_fiber = 100
# n_processes = n_fibers*n_processes_per_fiber
# n_fibers_per_process = n_fibers/n_processes

n = 9
# 1st case: multiple processes per fiber
for n_processes in np.logspace(4,2,n):   # 1e0 to 1e4, spaced evenly on a log scale.
  n_processes = int(np.round(n_processes))
  n_processes_per_fiber = n_processes/n_fibers
  n_processes_per_fiber = int(np.round(n_processes_per_fiber))
  
  run(n_processes, n_processes_per_fiber, n_fibers, n_nodes_per_fiber, "strong_scaling")
        
# 2nd case: multiple fibers per process
print("strong scaling, multiple fibers per process")
for n_processes in np.logspace(2,0,4):   # 1e0 to 1e4, spaced evenly on a log scale.
  n_processes = int(np.round(n_processes))
  if n_processes == 100:
    continue
  
  n_fibers_per_process = int(np.round(n_fibers/n_processes ))
  
  n_processes_per_fiber = -n_fibers_per_process
  run(n_processes, n_processes_per_fiber, n_fibers, n_nodes_per_fiber, "strong_scaling")
    
##########################
# strong scaling on multiple nodes, in factors of 24, one fiber
print("Strong scaling")

n_nodes_per_fiber = 10000*100
n_fibers = 1

# 1st case: multiple processes per fiber
for n_nodes in np.logspace(0,3,n):   # 1e01 to 1e4, spaced evenly on a log scale.
  n_nodes = int(np.round(n_nodes))
  n_processes = n_nodes*24
  n_processes_per_fiber = int(n_processes / n_fibers)
  #n_nodes_per_fiber = int(n_nodes_per_process * n_processes / n_fibers)
  
  run(n_processes, n_processes_per_fiber, n_fibers, n_nodes_per_fiber, "Strong_scaling")
        

##########################
# weak scaling
print("weak scaling")
n_nodes_per_fiber = 1000
n_processes_per_fiber = 10
# n_processes = n_fibers*n_processes_per_fiber

n = 9
for n_processes in np.logspace(0,4,n):   # 1e0 to 1e4, spaced evenly on a log scale.
  n_processes = int(np.round(n_processes))
  n_fibers = int(np.round(n_processes/n_processes_per_fiber))
  run(n_processes, n_processes_per_fiber, n_fibers, n_nodes_per_fiber, "weak_scaling")
    
