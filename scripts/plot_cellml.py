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
# from sets import Set # sets is already included in python3
import os
import time
import pickle
import py_reader    # reader utility for opendihu *.py files

files = ""

plot_name_condition = lambda filename: ".pdf" in filename

show_plot = True
if len(sys.argv) > 1:
  try:
    show_plot = int(sys.argv[1])
    files = sys.argv[2:]
  except:
    files = sys.argv[1:]
  plot_name_lst = list(np.extract(np.array(list(map(plot_name_condition, sys.argv))), sys.argv))
  if len(plot_name_lst)!=1:
    plot_name="fig.pdf"
  else:
    plot_name=plot_name_lst[0]
	
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

solution_files = list(np.extract(np.array(list(map(solution_condition, files))), files))
solution_shaped_files = list(np.extract(np.array(list(map(solution_shaped_condition, files))), files))
solution_py_files = list(np.extract(np.array(list(map(solution_py_condition, files))), files))


# sort files by number in file name (unnessesary?! sorted before in line 34)
solution_py_files = sorted(solution_py_files)

print( "{} files".format(len(solution_py_files)))
print( solution_py_files[0:min(10,len(solution_py_files))])

data = py_reader.load_data(solution_py_files)

if len(data) == 0:
  print( "no data found.")
  sys.exit(0)

dimension = data[0]['dimension']

dataVm=[]
time=[]
for ii in range(0,len(solution_py_files)):
  dataVm.append(data[ii]["data"][1]["components"][0]["values"][0])
  time.append(data[ii]["currentTime"])

min_value = min(dataVm)
max_value = max(dataVm)
min_x = time[0]
max_x = time[-1]

print( "value range: [{}, {}]".format(min_value, max_value))

# prepare plot
fig = plt.figure()

margin = abs(max_value - min_value) * 0.1
ax = plt.axes(xlim=(min_x, max_x), ylim=(min_value - margin, max_value + margin))
line, = ax.plot([], [], 'o-', lw=2)
text = plt.figtext(0.15,0.85,"V_m",size=20)
ax.set_xlabel('time')
ax.set_ylabel('V_m')
line.set_data(time, dataVm)
#print("plotname: {} of type {}".format(plot_name, type(plot_name)))
plt.savefig(plot_name)
    
if show_plot:
  plt.show()
        
sys.exit(0)
