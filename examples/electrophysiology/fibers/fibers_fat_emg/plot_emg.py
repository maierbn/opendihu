#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to visualize EMG output files
# arguments: <path>
#

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

path = "build_release/out"

# parse files
if len(sys.argv) > 1:
  filenames = sys.argv[1:]
else:
  # get all input data in current directory
  files_at_path = os.listdir(path)

  # sort files by number in file name
  filenames = sorted(files_at_path)

filenames = [f for f in filenames if ".csv" in f and "emg" in f]

x_coordinates = set()
y_coordinates = set()

# determine grid of electrodes
for filename in filenames:
  pos3 = filename.rfind(".csv")
  pos2 = filename.rfind("_")
  pos1 = filename.rfind("_",0,pos2)
  base = filename[0:pos1]
  
  y_coordinate = (int)(filename[pos1+1:pos2])
  x_coordinate = (int)(filename[pos2+1:pos3])

  x_coordinates.add(x_coordinate)
  y_coordinates.add(y_coordinate)
  
n_plots_x = len(x_coordinates)
n_plots_y = len(y_coordinates)

print("x:",x_coordinates, n_plots_x)
print("y:",y_coordinates, n_plots_y)

fig = plt.figure()

# labels for entire plot
ax = fig.add_subplot(111)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.xlabel("fiber direction")
plt.ylabel("cross-fiber direction")    

# loop over sub plots
for j,y in enumerate(y_coordinates):
  for i,x in enumerate(x_coordinates):
    n = j*n_plots_x + i + 1
    
    # get data from file
    filename = "{}/{}_{:02}_{:02}.csv".format(path, base,y,x)
    print("load filename: {}".format(filename))
    t,phi = np.loadtxt(filename, delimiter=';', skiprows=1, usecols=(0,1), unpack=True)
    
    ax = fig.add_subplot(n_plots_y,n_plots_x,n)
    ax.plot(t, phi)


plt.show()

